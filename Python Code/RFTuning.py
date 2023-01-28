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
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import seaborn as sns


fileList = ['221010', '221013', '221108', '221110', '221115', '221117',
            '221124', '221128', '221206', '221208', '221229', '230123',
            '230126']
unitParams = []
sessionXYMean = []
t0 = time.time()
for fileIterator in fileList:
    # Load relevant file here with pyMat reader
    monkeyName = 'Meetz'
    seshDate = fileIterator
    fileName = f'{monkeyName}_{seshDate}_GRF1_Spikes.mat'

    # ascertain the position of the gabor stimuli for MTNC
    fileNameMTNC = f'{monkeyName}_{seshDate}_MTNC_Spikes.mat'
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileNameMTNC)
    corrTrials = correctTrialsMTX(allTrials)
    currTrial = allTrials[corrTrials[0]]
    stimDesc = currTrial['stimDesc']['data']
    for count, i in enumerate(stimDesc['stimLoc']):
        if i == 0:
            stim0Azi = stimDesc['azimuthDeg'][count]
            stim0Ele = stimDesc['elevationDeg'][count]
        if i == 1:
            stim1Azi = stimDesc['azimuthDeg'][count]
            stim1Ele = stimDesc['elevationDeg'][count]
    os.chdir('../../Python Code')

    # Load the RF location file now
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

    # create folder and change dir to save PDF's and np.array
    if not os.path.exists('RFLoc Tuning'):
        os.makedirs('RFLoc Tuning')
    os.chdir('RFLoc Tuning/')

    # Tuning code
    correctTrials = correctTrialsGRF(allTrials)
    units = activeUnits('spikeData', allTrials)
    unitsChannel = unitsInfo(units, correctTrials, allTrials)

    # change stimDesc field to be a list of dictionaries
    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        nStim = len(currTrial['stimDesc']['data']['stimType'])
        currTrial['stimDesc']['data'] = [{k:v[i] for k,v in currTrial['stimDesc']['data'].items()}
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
    histPrePostMS = 100
    sponWindowMS = 100
    unitsXY = []

    # assert frame consistency for frame duration
    stimDurFrame = []
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
        stimCount = np.zeros((numEle,numAzi))
        spikeCountMat = np.zeros((25,numEle, numAzi))
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
        allTuningMat[uCount] = spikeCountMean - sponSpikesMean

        # 2D Gauss Fit
        spikeCountFitArray = spikeCountMean - sponSpikesMean
        com = ndimage.center_of_mass(spikeCountFitArray)
        p_init = models.Gaussian2D(amplitude=1, x_mean=com[1], y_mean=com[0],
                                   x_stddev=1, y_stddev=1, theta=None,
                                   cov_matrix=None)
        yi, xi = np.indices(spikeCountFitArray.shape)
        fit_p = fitting.LevMarLSQFitter()
        p = fit_p(p_init, xi, yi, spikeCountFitArray)

        theta = p.theta[0] * 180/np.pi
        xStdDev = p.x_stddev[0]
        yStdDev = p.y_stddev[0]
        xMean = p.x_mean[0] + 0.5
        yMean = p.y_mean[0] + 0.5
        amp = p.amplitude[0]
        xMeanConv = minAzi + ((maxAzi-minAzi)*(xMean/numAzi))
        yMeanConv = minEle + ((maxEle-minEle)*(yMean/numEle))
        xStdDevConv = xStdDev * ((maxAzi-minAzi)/numAzi)
        yStdDevConv = yStdDev * ((maxEle-minEle)/numEle)
        RFShapeParams = [theta, xStdDevConv, yStdDevConv,
                         xMeanConv, yMeanConv, amp]
        unitParams.append(RFShapeParams)
        unitsXY.append([xMeanConv, yMeanConv])

        # # Plot figure
        # fig = plt.figure()
        # fig.set_size_inches(6, 8)
        # aziLabel = np.around(np.linspace(minAzi, maxAzi, numAzi), 2)
        # eleLabel = np.around(np.linspace(minEle, maxEle, numEle), 2)
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
        # # heatmap
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
        #                    s=100, marker='o', color='red')
        # ax_row1[0].scatter(((stim1Azi + abs(minAzi)) / (maxAzi + abs(minAzi)) * numAzi),
        #                    (stim1Ele + abs(minEle)) / (abs(minEle) + maxEle) * numEle,
        #                    s=100, marker='o', color='green')
        #
        # ax_row1[0].invert_yaxis()
        #
        # # overlay 1SD, 2SD ellipse
        # el1SD = Ellipse((xMean, yMean), xStdDev, yStdDev, theta,
        #                 fill=None, edgecolor='blue')
        # el1SD.set(linestyle=':')
        # ax_row1[0].add_artist(el1SD)
        # el2SD = Ellipse((xMean, yMean), 2 * xStdDev, 2 * yStdDev, theta,
        #                 fill=None, edgecolor='black')
        # el2SD.set(linestyle='--')
        # ax_row1[0].add_artist(el2SD)
        #
        # ax_row2 = []
        # for i in np.arange(numEle + 2, 2, -1):
        #     ax = []
        #     for j in np.arange(0, numAzi):
        #         ax.append(plt.subplot2grid((numEle + 3, 6), (i, j)))
        #     ax_row2.append(np.array(ax))
        # ax_row2 = np.array(ax_row2)  # 6 x 6
        #
        # # PSTHs
        # yMax = 0
        # for i in range(numEle):
        #     for j in range(numAzi):
        #         spikeHist = spikeHists[:, i, j] * 1000 / stimCount[i, j]
        #         gaussSmooth = gaussian_filter1d(spikeHist, 5)
        #         if max(gaussSmooth) > yMax:
        #             yMax = max(gaussSmooth)
        #         ax_row2[i, j].plot(gaussSmooth)
        #         # ax_row2[i,j].set_title(f"Ele:{eleLabel[i]},Azi:{aziLabel[j]}", fontsize=4)
        #         # ax_row2[i,j].set_ylim(bottom=0)
        #         # ax_row2[i,j].set_ylim([0, 100])
        #         # ax_row2[i,j].set_yticks([0,50,100])
        #         # ax_row2[i,j].set_yticklabels([0,50,100], fontsize=5)
        #         ax_row2[i, j].tick_params(axis='y', which='major', labelsize=4)
        #         ax_row2[i, j].set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
        #                                   2 * histPrePostMS + trueStimDurMS])
        #         ax_row2[i, j].set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
        #                                        trueStimDurMS + histPrePostMS],
        #                                       fontsize=3)
        #         ax_row2[i, j].axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
        #                               color='grey', alpha=0.2)
        #         ax_row2[i, j].axhline(y=sponSpikesMean * 1000 / sponWindowMS,
        #                               linestyle='--', color='grey')
        #         if i == 0 and j == 0:
        #             ax_row2[i, j].set_xlabel('Time (ms)', fontsize=7)
        #             ax_row2[i, j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
        #             ax_row2[i, j].yaxis.set_label_coords(-0.5, 1.75)
        # plt.tight_layout(pad=0.1, w_pad=-0.25, h_pad=-1.45)
        #
        # for i in range(numEle):
        #     for j in range(numAzi):
        #         ax_row2[i, j].set_ylim([0, yMax * 1.1])
        #
        # # saves plot as pdf
        # plt.savefig(f'{unit}.pdf')
        # plt.close('all')
        # continue

    # calculate average XY RF position for units from session, excluding
    # units with values more than 1.5 SD away from the mean
    unitsXY = np.array(unitsXY)
    unitsXFiltered = rejectOutliers(unitsXY[:, 0], m=2)
    unitsYFiltered = rejectOutliers(unitsXY[:, 1], m=2)
    unitsXYMean = [np.mean(unitsXFiltered), np.mean(unitsYFiltered)]
    sessionXYMean.append(unitsXYMean)

    plt.close('all')
    np.save('unitsRFLocMat', allTuningMat)
    os.chdir('../../../Python Code')

print(time.time()-t0)

# Plot Overall RF Summary from all sessions
# figure
fig = plt.figure()
fig.set_size_inches(15, 7)
gs0 = gridspec.GridSpec(1, 2)
gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0])
gs01 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])

ax = fig.add_subplot(gs00[0, 0])
x = np.linspace(-15, 30, 10)
y = np.zeros(len(x))
plt.plot(x, y, color='black')
y = np.linspace(10, -30, 9)
x = np.zeros(len(y))
plt.plot(x, y, color='black')

# overlay 1SD
for i in unitParams:
    if abs(i[3]) < 30 and abs(i[4]) < 30:
        el1SD = Ellipse((i[3], i[4]), i[1], i[2], i[0],
                        fill=None, edgecolor='blue', alpha=0.4)
        el1SD.set(linestyle=':')
        stdMean = ((i[1]*2) + (i[2]*2)) / 2
        if stdMean < 40:
            ax.add_artist(el1SD)
ax.set_xlabel('Azimuth (˚)')
ax.set_ylabel('Elevation (˚)')
sessionXYMean = np.array(sessionXYMean)
ax.scatter(sessionXYMean[:, 0], sessionXYMean[:, 1])

ax2 = fig.add_subplot(gs01[0, 0])
for i in unitParams:
    if abs(i[3]) < 30 and abs(i[4]) < 30:
        unitEcc = np.sqrt((i[3]**2) + (i[4]**2))
        stdMean = ((i[1]*2*np.sqrt(2)) + (i[2]*2*np.sqrt(2))) / 2
        if stdMean < 40:
            ax2.scatter(unitEcc, stdMean, color='grey', alpha=0.8)
ax2.set_xlim([0, 40])
ax2.set_ylim([0, 40])
line = lines.Line2D([0, 1], [0, 1], color='black')
transform = ax2.transAxes
line.set_transform(transform)
ax2.add_line(line)
ax2.set_xlabel('RF Eccentricity (˚)')
ax2.set_ylabel('RF Diameter (2DGaussSD * 2 * sqrt(2))')

plt.tight_layout()
plt.show()


################################## STOP HERE ###########################################

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

### for testing purposes, to make unit field similar to real data
for currTrial in allTrials:
    if 'spikeData' in currTrial:
        currTrial['spikeData']['unit'] = currTrial['spikeData']['unit'].tolist()
        for i in range(0,len(currTrial['spikeData']['channel'])):
            a = str(int(currTrial['spikeData']['channel'][i]))
            b = str(int(currTrial['spikeData']['unit'][i]))
            c = a + '_' + b
            currTrial['spikeData']['unit'][i] = c
        currTrial['spikeData']['unit'] = np.array(currTrial['spikeData']['unit'])
###