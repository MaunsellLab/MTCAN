'''
RFTuning generates heatmaps and histogram of each unit's RF location.
The script will save the plots for each unit as a PDF in a folder specific
to the day's dataset. 

to do:

Chery March 2022
Modified to save plots as pngs and incorporated changes to track stimCounts
- 04/26/22
Modified to incorporate frame counter, fixed bugs relating to unit indexing
within the stimDesc, added 1-d gaussian filter - 06/10/22
Flipped heatmap and corresponding PSTHs, - 06/12/22
Import using pymat reader - 10/11/22
2D Gauss fit overlay on heatmap (1 sd, 2sd) - 10/18/22
'''

from usefulFns import *
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import seaborn as sns

fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
            'Meetz_221115', 'Meetz_221117',  'Meetz_221124', 'Meetz_221128',
            'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
            'Meetz_230126',  'Akshan_240927', 'Akshan_241016', 'Akshan_241017',
            'Akshan_241025']

unitParams = []
sessionXYMean = []
loc0XY = []
loc1XY = []
loc0XYDiffFromMean = []
loc1XYDiffFromMean = []
loc0DistFromRFCent = []
loc1DistFromRFCent = []
sessionRFGaborSigma = []
t0 = time.time()

for fileIterator in fileList:
    # Load relevant file here with pyMat reader
    monkeyName, seshDate = fileIterator.split('_')
    fileName = f'{monkeyName}_{seshDate}_GRF1_Spikes.mat'

    # ascertain the position of the gabor stimuli for
    # and get the sigma (size) of the Gabor used for MTNC
    fileNameMTNC = f'{monkeyName}_{seshDate}_MTNC_Spikes.mat'
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileNameMTNC)
    corrTrials = correctTrialsMTX(allTrials)
    currTrial = allTrials[corrTrials[0]]
    rfGaborSigma = currTrial['rfGabor']['data']['sigmaDeg']
    sessionRFGaborSigma.append(rfGaborSigma)
    stimDesc = currTrial['stimDesc']['data']
    for count, i in enumerate(stimDesc['stimLoc']):
        if i == 0:
            stim0Azi = stimDesc['azimuthDeg'][count]
            stim0Ele = stimDesc['elevationDeg'][count]
        if i == 1:
            stim1Azi = stimDesc['azimuthDeg'][count]
            stim1Ele = stimDesc['elevationDeg'][count]

    # Load the RF location file now
    os.chdir('../../Python Code')
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

    # create folder and change dir to save PDF's and np.array
    if not os.path.exists('RFLoc Tuning'):
        os.makedirs('RFLoc Tuning')
    os.chdir('RFLoc Tuning/')

    # Tuning code
    correctTrials = correctTrialsGRF(allTrials)
    units = activeUnits('spikeData', allTrials)
    unitsChannel = unitsInfo(units, correctTrials, allTrials)
    seshUnitsID = []
    for unit in units:
        seshUnitsID.append(f'{seshDate}_{unit}')
    seshUnitsID = np.array(seshUnitsID)

    # change stimDesc field to be a list of dictionaries
    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        nStim = len(currTrial['stimDesc']['data']['stimType'])
        currTrial['stimDesc']['data'] = [{k:v[i] for k, v in currTrial['stimDesc']['data'].items()}
                                         for i in range(nStim)]

    frameRateHz = header['displayCalibration']['data']['frameRateHz']
    numAzi = header['map0Settings']['data']['azimuthDeg']['n']
    minAzi = np.float32(header['map0Settings']['data']['azimuthDeg']['minValue'])
    maxAzi = np.float32(header['map0Settings']['data']['azimuthDeg']['maxValue'])
    numEle = header['map0Settings']['data']['elevationDeg']['n']
    minEle = np.float32(header['map0Settings']['data']['elevationDeg']['minValue'])
    maxEle = np.float32(header['map0Settings']['data']['elevationDeg']['maxValue'])
    interstimDurMS = header['mapInterstimDurationMS']['data']
    allTuningMat = np.zeros((len(units), numEle, numAzi))
    allBaselineMat = np.zeros(len(units))
    histPrePostMS = 100
    sponWindowMS = 100
    unitsXY = []

    # assert frame consistency for frame duration
    stimDurFrame = []
    trueStimDurMS = 0
    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        stimDesc = currTrial['stimDesc']['data']
        if 'numMap0Stim' in currTrial:
            map0StimLim = currTrial['numMap0Stim']['data']
            map0Count = 0
            for stim in stimDesc:
                if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                    frameDiff = stim['stimOffFrame'] - stim['stimOnFrame']
                    stimDurFrame.append(frameDiff)
    if len(set(stimDurFrame)) != 1:
        print('stimulus frame duration not consistent across mapping stimuli')
    else:
        trueStimDurMS = np.int32(np.around(1000/frameRateHz * stimDurFrame[0]))

    for uCount, unit in enumerate(units):
        stimCount = np.zeros((numEle, numAzi))
        spikeCountMat = np.zeros((25, numEle, numAzi))
        spikeCountMat[:, :, :] = np.nan
        spikeHists = np.zeros((trueStimDurMS + 2*histPrePostMS+12, numEle, numAzi))
        sponSpikesArr = []

        for corrTrial in correctTrials:
            currTrial = allTrials[corrTrial]
            # stimCount verify
            if 'numMap0Stim' in currTrial:
                map0StimLim = currTrial['numMap0Stim']['data']
                map0Count = 0
                stimDesc = currTrial['stimDesc']['data']
                spikeData = currTrial['spikeData']
                stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
                for stim in stimDesc:
                    # extract spikes from each stimulus aligned to frame counter (on/off)
                    if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                        aziIndex = int(stim['azimuthIndex'])
                        eleIndex = int(stim['elevationIndex'])
                        stCount = int(stimCount[eleIndex][aziIndex])
                        stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'])
                                       / 1000) + stim1TimeS
                        stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'])
                                        / 1000) + stim1TimeS
                        stimCount[eleIndex][aziIndex] += 1
                        map0Count += 1
                        if unit in spikeData['unit']:
                            spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                            unitIndex = np.where(spikeData['unit'] == unit)[0]
                            unitTimeStamps = spikeData['timeStamp'][unitIndex]
                            stimSpikes = np.where((unitTimeStamps >= stimOnTimeS) &
                                                  (unitTimeStamps <= stimOffTimeS))
                            spikeCountMat[stCount][eleIndex][aziIndex] = len(stimSpikes[0])
                            # spontaneous spike count (x ms before stim on)
                            sponSpikes = np.where((unitTimeStamps >= (stimOnTimeS-(sponWindowMS/1000))) &
                                                  (unitTimeStamps <= stimOnTimeS))
                            sponSpikesArr.extend([len(sponSpikes[0])])

                            histStimSpikes = histSpikes(stimOnTimeS, stimOffTimeS,
                                                        histPrePostMS, unitTimeStamps)
                            spikeHists[histStimSpikes, eleIndex, aziIndex] += 1

        spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis=0)
        sponSpikesArr = np.array(sponSpikesArr)
        sponSpikesMean = np.mean(sponSpikesArr)
        allTuningMat[uCount] = spikeCountMean
        allBaselineMat[uCount] = sponSpikesMean

        # ---- 2D Gauss Fit (fail gracefully) ----
        fit_ok = True
        try:
            spikeCountFitArray = spikeCountMean - sponSpikesMean
            com = ndimage.center_of_mass(spikeCountFitArray)
            p_init = models.Gaussian2D(
                amplitude=1, x_mean=com[1], y_mean=com[0],
                x_stddev=1, y_stddev=1, theta=0.0  # theta=None can be finicky; 0.0 is safe
            )
            yi, xi = np.indices(spikeCountFitArray.shape)
            fit_p = fitting.LevMarLSQFitter()
            p = fit_p(p_init, xi, yi, spikeCountFitArray)

            # Extract parameters in a version-agnostic way (.value or float())
            theta = float(getattr(p.theta, "value", p.theta)) * 180 / np.pi
            xStdDev = float(getattr(p.x_stddev, "value", p.x_stddev))
            yStdDev = float(getattr(p.y_stddev, "value", p.y_stddev))
            xMean = float(getattr(p.x_mean, "value", p.x_mean)) + 0.5
            yMean = float(getattr(p.y_mean, "value", p.y_mean)) + 0.5
            amp = float(getattr(p.amplitude, "value", p.amplitude))

            xMeanConv = minAzi + ((maxAzi - minAzi) * (xMean / numAzi))
            yMeanConv = minEle + ((maxEle - minEle) * (yMean / numEle))
            xStdDevConv = xStdDev * ((maxAzi - minAzi) / numAzi)
            yStdDevConv = yStdDev * ((maxEle - minEle) / numEle)

        except Exception:
            # On any failure: mark not OK and fill NaNs
            fit_ok = False
            theta = np.nan
            xStdDev = np.nan
            yStdDev = np.nan
            xMean = np.nan
            yMean = np.nan
            amp = np.nan
            xMeanConv = np.nan
            yMeanConv = np.nan
            xStdDevConv = np.nan
            yStdDevConv = np.nan

        # Always append; NaNs on failure
        RFShapeParams = [theta, xStdDevConv, yStdDevConv, xMeanConv, yMeanConv, amp]
        unitParams.append(RFShapeParams)
        unitsXY.append([xMeanConv, yMeanConv])

        # ---- find position of loc 0 and 1 wrt to center of RF (will propagate NaNs if fit failed) ----
        loc0XY.append((stim0Azi, stim0Ele))
        loc1XY.append((stim1Azi, stim1Ele))
        loc0XYDiffFromMean.append(((stim0Azi - xMeanConv), (stim0Ele - yMeanConv)))
        loc1XYDiffFromMean.append(((stim1Azi - xMeanConv), (stim1Ele - yMeanConv)))
        loc0OffsetFromCenter = np.sqrt(((stim0Azi - xMeanConv) ** 2) + ((stim0Ele - yMeanConv) ** 2))
        loc1OffsetFromCenter = np.sqrt(((stim1Azi - xMeanConv) ** 2) + ((stim1Ele - yMeanConv) ** 2))
        loc0DistFromRFCent.append(loc0OffsetFromCenter)
        loc1DistFromRFCent.append(loc1OffsetFromCenter)

        # # ---- Plot figure ----
        # fig = plt.figure()
        # fig.set_size_inches(6, 8)
        aziLabel = np.around(np.linspace(minAzi, maxAzi, numAzi), 2)
        eleLabel = np.around(np.linspace(minEle, maxEle, numEle), 2)
        #
        # date = header['date']
        #
        # text = fig.text(0.05, 0.75, f'RF tuning for unit {unit}\n\
        # {date}\n\
        # - - - - - - - - - - - - - - - \n\
        # Stimulus Duration: {trueStimDurMS} ms\n\
        # Number of Blocks: {int(stimCount[0][0])}\n\
        # Interstimulus Duration: {interstimDurMS}\n\
        # Channel: {unitsChannel[uCount]}\n\
        # Stim0 Azi: {stim0Azi:.1f}, Ele: {stim0Ele:.1f}\n\
        # Stim1 Azi: {stim1Azi:.1f}, Ele: {stim1Ele:.1f}', size=10, fontweight='bold')
        # text.set_path_effects([path_effects.Normal()])
        #
        # # heatmap (raw data)
        # ax_row1 = []
        # ax_row1.append(plt.subplot2grid((numEle + 3, 6), (0, 3), colspan=3, rowspan=3))
        # # bSmooth = gaussian_filter1d(spikeCountMean, sigma=0.5)
        # bSmooth = gaussian_filter(spikeCountMean, sigma=1)
        # ax_row1[0] = sns.heatmap(spikeCountMean * 1000 / trueStimDurMS, vmin=0)
        # # ax_row1[0].contour(np.arange(.5, bSmooth.shape[1]), np.arange(.5,
        # #                    bSmooth.shape[0]), bSmooth, colors='yellow')
        # ax_row1[0].set_xlabel('azimuth (˚)', fontsize=8)
        # ax_row1[0].set_ylabel('elevation (˚)', fontsize=8)
        # ax_row1[0].set_title('Heatmap of unit RF location', fontsize=9)
        # ax_row1[0].set_yticklabels(eleLabel, fontsize=5)
        # ax_row1[0].set_xticklabels(aziLabel, fontsize=5)
        #
        # # marker lines
        # ax_row1[0].scatter(((stim0Azi + abs(minAzi)) / (maxAzi + abs(minAzi)) * numAzi),
        #                    (stim0Ele + abs(minEle)) / (abs(minEle) + maxEle) * numEle,
        #                    s=100, marker='o')
        # ax_row1[0].scatter(((stim1Azi + abs(minAzi)) / (maxAzi + abs(minAzi)) * numAzi),
        #                    (stim1Ele + abs(minEle)) / (abs(minEle) + maxEle) * numEle,
        #                    s=100, marker='o')
        # ax_row1[0].invert_yaxis()
        #
        # # overlay 1SD, 2SD ellipses ONLY if fit succeeded
        # if fit_ok:
        #     el1SD = Ellipse((xMean, yMean), xStdDev, yStdDev, theta, fill=None, edgecolor='blue')
        #     el1SD.set(linestyle=':')
        #     ax_row1[0].add_artist(el1SD)
        #     el2SD = Ellipse((xMean, yMean), 2 * xStdDev, 2 * yStdDev, theta, fill=None, edgecolor='black')
        #     el2SD.set(linestyle='--')
        #     ax_row1[0].add_artist(el2SD)
        #
        # # ---- PSTHs ----
        # ax_row2 = []
        # for i in np.arange(numEle + 2, 2, -1):
        #     ax = []
        #     for j in np.arange(0, numAzi):
        #         ax.append(plt.subplot2grid((numEle + 3, 6), (i, j)))
        #     ax_row2.append(np.array(ax))
        # ax_row2 = np.array(ax_row2)  # 6 x 6
        #
        # yMax = 0
        # for i in range(numEle):
        #     for j in range(numAzi):
        #         spikeHist = spikeHists[:, i, j] * 1000 / stimCount[i, j]
        #         gaussSmooth = gaussian_filter1d(spikeHist, 5)
        #         if np.max(gaussSmooth) > yMax:
        #             yMax = np.max(gaussSmooth)
        #         ax_row2[i, j].plot(gaussSmooth)
        #         ax_row2[i, j].tick_params(axis='y', which='major', labelsize=4)
        #         ax_row2[i, j].set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
        #                                   2 * histPrePostMS + trueStimDurMS])
        #         ax_row2[i, j].set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
        #                                        trueStimDurMS + histPrePostMS],
        #                                       fontsize=3)
        #         ax_row2[i, j].axvspan(histPrePostMS, histPrePostMS + trueStimDurMS, color='grey', alpha=0.2)
        #         ax_row2[i, j].axhline(y=sponSpikesMean * 1000 / sponWindowMS, linestyle='--', color='grey')
        #         if i == 0 and j == 0:
        #             ax_row2[i, j].set_xlabel('Time (ms)', fontsize=7)
        #             ax_row2[i, j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
        #             ax_row2[i, j].yaxis.set_label_coords(-0.5, 1.75)
        #
        # plt.tight_layout(pad=0.1, w_pad=-0.25, h_pad=-1.45)
        #
        # for i in range(numEle):
        #     for j in range(numAzi):
        #         ax_row2[i, j].set_ylim([0, yMax * 1.1])
        #
        # # ---- save & close ----
        # plt.savefig(f'{unit}.pdf')
        # plt.close('all')

    # calculate average XY RF position for units from session, excluding
    # units with values more than 1.5 SD away from the mean
    unitsXY = np.array(unitsXY)
    unitsXFiltered = rejectOutliers(unitsXY[:, 0], m=2)
    unitsYFiltered = rejectOutliers(unitsXY[:, 1], m=2)
    unitsXYMean = [np.mean(unitsXFiltered), np.mean(unitsYFiltered)]
    sessionXYMean.append(unitsXYMean)

    plt.close('all')
    np.save('unitsRFLocMat.npy', allTuningMat)
    np.save('unitsBaseline.npy', allBaselineMat)
    np.save('seshUnitsID.npy', seshUnitsID)
    seshAziEle = np.array([(azi, ele) for azi in aziLabel for ele in eleLabel])
    np.save("seshAziEle.npy", seshAziEle)
    seshStimLocs = np.array([
        [stim0Azi, stim0Ele],  # stimulus 0: [azi, ele]
        [stim1Azi, stim1Ele]  # stimulus 1: [azi, ele]
    ])
    np.save("seshStimLocs.npy", seshStimLocs)

    os.chdir('../../../Python Code')

print(time.time()-t0)

# # Plot Overall RF Summary from all sessions
# loc0XYDiffFromMean = np.array(loc0XYDiffFromMean)
# loc1XYDiffFromMean = np.array(loc1XYDiffFromMean)
# unitParams = np.array(unitParams)
#
# # figure
# fig = plt.figure()
# fig.set_size_inches(10, 7)
# gs0 = gridspec.GridSpec(2, 2)
# gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0, 0])
# gs01 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0, 1])
# gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1, 0])
# gs03 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1, 1])
#
# ax = fig.add_subplot(gs00[0, 0])
# x = np.linspace(-15, 30, 10)
# y = np.zeros(len(x))
# plt.plot(x, y, color='black')
# y = np.linspace(10, -30, 9)
# x = np.zeros(len(y))
# plt.plot(x, y, color='black')
#
# # overlay 1SD
# badUnits = []
# goodUnits = []
# for count, i in enumerate(unitParams):
#     if -5 < i[3] < 30 and 3 > i[4] > -20:
#         el1SD = Ellipse((i[3], i[4]), i[1], i[2], i[0],
#                         fill=None, edgecolor='blue', alpha=0.2)
#         el1SD.set(linestyle=':')
#         stdMean = ((i[1]*2) + (i[2]*2)) / 2
#         if stdMean < 40:
#             goodUnits.append(count)
#             ax.add_artist(el1SD)
#             ax.scatter(loc0XY[count][0], loc0XY[count][1], color='red')
#             ax.scatter(loc1XY[count][0], loc1XY[count][1], color='blue')
#     else:
#         badUnits.append(count)
# ax.set_xlabel('Azimuth (˚)')
# ax.set_ylabel('Elevation (˚)')
# sessionXYMean = np.array(sessionXYMean)
# ax.scatter(sessionXYMean[:, 0], sessionXYMean[:, 1], alpha=0.6, color='green')
#
# for i in range(len(sessionXYMean)):
#     ax.annotate(i+1, (sessionXYMean[i, 0], sessionXYMean[i, 1]), color='black',
#                 fontsize=6)
#
#
# ax2 = fig.add_subplot(gs01[0, 0])
# diamEccRatio = []
# for i in unitParams:
#     if abs(i[3]) < 30 and abs(i[4]) < 30:
#         unitEcc = np.sqrt((i[3]**2) + (i[4]**2))
#         diam = ((i[1]*2*np.sqrt(2)) + (i[2]*2*np.sqrt(2))) / 2
#         if diam < 40:
#             ax2.scatter(unitEcc, diam, color='grey', alpha=0.8)
#             diamEccRatio.append(diam/unitEcc)
#
# ax2.set_xlim([0, 40])
# ax2.set_ylim([0, 40])
# line = lines.Line2D([0, 1], [0, 1], color='black')
# transform = ax2.transAxes
# line.set_transform(transform)
# ax2.add_line(line)
# ax2.set_xlabel('RF Eccentricity (˚)')
# ax2.set_ylabel('RF Diameter (2DGaussSD * 2 * sqrt(2))')
#
# # ax3 = fig.add_subplot(gs02[0, 0])
#
# ax4 = fig.add_subplot(gs03[0, 0])
# ax4.hist(diamEccRatio, bins=20)
# ax4.set_ylabel('Frequency')
# ax4.set_xlabel('unit Diameter/Eccentricity ratio')
# ax4.axvline(np.median(diamEccRatio), color='black')
# ax4.set_title('Diameter/Eccentricity distribution, vertical line=median',
#               fontsize=7)
#
# plt.tight_layout()
# plt.show()
#
# # plot average RF ellipse with average loc 0 and 1 stim placement
# loc0XAvgOffset = np.median(loc0XYDiffFromMean[goodUnits, :], 0)[0]
# loc0YAvgOffset = np.median(loc0XYDiffFromMean[goodUnits, :], 0)[1]
# loc1XAvgOffset = np.median(loc1XYDiffFromMean[goodUnits, :], 0)[0]
# loc1YAvgOffset = np.median(loc1XYDiffFromMean[goodUnits, :], 0)[1]
# avgEllipseXStd = np.median(unitParams[goodUnits, 1], 0)
# avgEllipseYStd = np.median(unitParams[goodUnits, 2], 0)
#
# fig, ax = plt.subplots()
# x = np.linspace(-5, 5, 11)
# y = np.zeros(len(x))
# plt.plot(x, y, color='white', alpha=0)
# y = np.linspace(5, -5, 21)
# x = np.zeros(len(y))
# plt.plot(x, y, color='white', alpha=0)
# el1SD = Ellipse((0, 0), avgEllipseXStd, avgEllipseYStd,
#                 fill=None, edgecolor='blue', alpha=0.2)
# ax.add_artist(el1SD)
# ax.scatter(0, 0, marker='X')
# ax.scatter(loc0XAvgOffset, loc0YAvgOffset, color='red', label='loc 0')
# ax.scatter(loc1XAvgOffset, loc1YAvgOffset, color='blue', label='loc 1')
# ax.legend()
# plt.show()
#
# # plot RF Gabor sigma as a function of eccentricity
# sessionXYMean = np.array(sessionXYMean)
# sessionRFGaborSigma = np.array(sessionRFGaborSigma)
#
# sessionEccentricity = np.sqrt((sessionXYMean[:, 0] ** 2) +
#                               (sessionXYMean[:, 1] ** 2))
#
# fig, ax = plt.subplots()
# ax.scatter(sessionEccentricity, sessionRFGaborSigma)
# ax.set_xlabel('Session Eccentricity', fontsize=15)
# ax.set_ylabel('Session RF Gabor Sigma', fontsize=15)
# plt.show()
#
#
#
#
# ################################## STOP HERE ###########################################
#
# def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
#     """
#     Create a plot of the covariance confidence ellipse of *x* and *y*.
#
#     Parameters
#     ----------
#     x, y : array-like, shape (n, )
#         Input data.
#
#     ax : matplotlib.axes.Axes
#         The axes object to draw the ellipse into.
#
#     n_std : float
#         The number of standard deviations to determine the ellipse's radiuses.
#
#     **kwargs
#         Forwarded to `~matplotlib.patches.Ellipse`
#
#     Returns
#     -------
#     matplotlib.patches.Ellipse
#     """
#     if x.size != y.size:
#         raise ValueError("x and y must be the same size")
#
#     cov = np.cov(x, y)
#     pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
#     # Using a special case to obtain the eigenvalues of this
#     # two-dimensionl dataset.
#     ell_radius_x = np.sqrt(1 + pearson)
#     ell_radius_y = np.sqrt(1 - pearson)
#     ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
#                       facecolor=facecolor, **kwargs)
#
#     # Calculating the stdandard deviation of x from
#     # the squareroot of the variance and multiplying
#     # with the given number of standard deviations.
#     scale_x = np.sqrt(cov[0, 0]) * n_std
#     mean_x = np.mean(x)
#
#     # calculating the stdandard deviation of y ...
#     scale_y = np.sqrt(cov[1, 1]) * n_std
#     mean_y = np.mean(y)
#
#     transf = transforms.Affine2D() \
#         .rotate_deg(45) \
#         .scale(scale_x, scale_y) \
#         .translate(mean_x, mean_y)
#
#     ellipse.set_transform(transf + ax.transData)
#     return ax.add_patch(ellipse)
#
# ### for testing purposes, to make unit field similar to real data
# for currTrial in allTrials:
#     if 'spikeData' in currTrial:
#         currTrial['spikeData']['unit'] = currTrial['spikeData']['unit'].tolist()
#         for i in range(0,len(currTrial['spikeData']['channel'])):
#             a = str(int(currTrial['spikeData']['channel'][i]))
#             b = str(int(currTrial['spikeData']['unit'][i]))
#             c = a + '_' + b
#             currTrial['spikeData']['unit'][i] = c
#         currTrial['spikeData']['unit'] = np.array(currTrial['spikeData']['unit'])
# ###