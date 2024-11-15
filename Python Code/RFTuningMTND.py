"""
RFTuning generates heatmaps and histogram of each unit's RF location.
The script will save the plots for each unit as a PDF in a folder specific
to the day's dataset.

Chery June 2023
"""

from usefulFns import *
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import seaborn as sns

#### MEETZ
# goodSessions = ['230607', '230608', '230612', '230614']
#
# fileList = ['230720']

#### AKSHAN

# goodSessions = ['240529', '240530', '240603', '240605', '240606', '240607',
#                 '240610', '240611', '240613', '240628', '240703, '240704',
#                 '240705']

fileList = ['241112']

t0 = time.time()

for fileIterator in fileList:
    # Load relevant file here with pyMat reader
    monkeyName = 'Akshan' #'Meetz'
    seshDate = fileIterator
    fileName = f'{monkeyName}_{seshDate}_GRF1_Spikes.mat'
    # fileName = f'{monkeyName}_{seshDate}_GRF1_2_Spikes.mat'

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
        currTrial['stimDesc']['data'] = [{k:v[i] for k, v in currTrial['stimDesc']['data'].items()}
                                         for i in range(nStim)]

    frameRateHz = header['displayCalibration']['data']['frameRateHz']
    numAzi = header['map0Settings']['data']['azimuthDeg']['n']
    minAzi = np.float32(header['map0Settings']['data']['azimuthDeg']['minValue'])
    maxAzi = np.float32(header['map0Settings']['data']['azimuthDeg']['maxValue'])
    numEle = header['map0Settings']['data']['elevationDeg']['n']
    minEle = np.float32(header['map0Settings']['data']['elevationDeg']['minValue'])
    maxEle = np.float32(header['map0Settings']['data']['elevationDeg']['maxValue'])
    aziLabel = np.around(np.linspace(minAzi, maxAzi, numAzi), 2)
    eleLabel = np.around(np.linspace(minEle, maxEle, numEle), 2)
    interstimDurMS = header['mapInterstimDurationMS']['data']
    allTuningMat = np.zeros((len(units), numEle, numAzi))
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
        allTuningMat[uCount] = (spikeCountMean - sponSpikesMean) * 1000 / trueStimDurMS

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

        unitsXY.append([xMeanConv, yMeanConv])

        # Plot figure
        fig = plt.figure()
        fig.set_size_inches(6, 8)
        aziLabel = np.around(np.linspace(minAzi, maxAzi, numAzi), 2)
        eleLabel = np.around(np.linspace(minEle, maxEle, numEle), 2)

        date = header['date']

        text = fig.text(0.05, 0.75, f'RF tuning for unit {unit}\n\
        {date}\n\
        - - - - - - - - - - - - - - - \n\
        Stimulus Duration: {trueStimDurMS} ms\n\
        Number of Blocks: {int(stimCount[0][0])}\n\
        Interstimulus Duration: {interstimDurMS}\n\
        Channel: {unitsChannel[uCount]}', size=10, fontweight='bold')
        text.set_path_effects([path_effects.Normal()])

        # heatmap
        ax_row1 = []
        ax_row1.append(plt.subplot2grid((numEle + 3, 6), (0, 3), colspan=3, rowspan=3))
        # bSmooth = gaussian_filter1d(spikeCountMean, sigma=0.5)
        bSmooth = gaussian_filter(spikeCountMean, sigma=1)
        ax_row1[0] = sns.heatmap(spikeCountMean * 1000 / trueStimDurMS, vmin=0)
        # ax_row1[0].contour(np.arange(.5, bSmooth.shape[1]), np.arange(.5,
        #                    bSmooth.shape[0]), bSmooth, colors='yellow')
        ax_row1[0].set_xlabel('azimuth (˚)', fontsize=8)
        ax_row1[0].set_ylabel('elevation (˚)', fontsize=8)
        ax_row1[0].set_title('Heatmap of unit RF location', fontsize=9)
        ax_row1[0].set_yticklabels(eleLabel, fontsize=5)
        ax_row1[0].set_xticklabels(aziLabel, fontsize=5)

        # marker lines
        ax_row1[0].invert_yaxis()

        # overlay 1SD, 2SD ellipse
        el1SD = Ellipse((xMean, yMean), xStdDev, yStdDev, angle=theta,
                        fill=None, edgecolor='blue')
        el1SD.set(linestyle=':')
        ax_row1[0].add_artist(el1SD)
        el2SD = Ellipse((xMean, yMean), 2 * xStdDev, 2 * yStdDev, angle=theta,
                        fill=None, edgecolor='black')
        el2SD.set(linestyle='--')
        ax_row1[0].add_artist(el2SD)

        ax_row2 = []
        for i in np.arange(numEle + 2, 2, -1):
            ax = []
            for j in np.arange(0, numAzi):
                ax.append(plt.subplot2grid((numEle + 3, 6), (i, j)))
            ax_row2.append(np.array(ax))
        ax_row2 = np.array(ax_row2)  # 6 x 6

        # PSTHs
        yMax = 0
        for i in range(numEle):
            for j in range(numAzi):
                spikeHist = spikeHists[:, i, j] * 1000 / stimCount[i, j]
                gaussSmooth = gaussian_filter1d(spikeHist, 5)
                if max(gaussSmooth) > yMax:
                    yMax = max(gaussSmooth)
                ax_row2[i, j].plot(gaussSmooth)
                # ax_row2[i,j].set_title(f"Ele:{eleLabel[i]},Azi:{aziLabel[j]}", fontsize=4)
                # ax_row2[i,j].set_ylim(bottom=0)
                # ax_row2[i,j].set_ylim([0, 100])
                # ax_row2[i,j].set_yticks([0,50,100])
                # ax_row2[i,j].set_yticklabels([0,50,100], fontsize=5)
                ax_row2[i, j].tick_params(axis='y', which='major', labelsize=4)
                ax_row2[i, j].set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                                          2 * histPrePostMS + trueStimDurMS])
                ax_row2[i, j].set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                                               trueStimDurMS + histPrePostMS],
                                              fontsize=3)
                ax_row2[i, j].axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                                      color='grey', alpha=0.2)
                ax_row2[i, j].axhline(y=sponSpikesMean * 1000 / sponWindowMS,
                                      linestyle='--', color='grey')
                if i == 0 and j == 0:
                    ax_row2[i, j].set_xlabel('Time (ms)', fontsize=7)
                    ax_row2[i, j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                    ax_row2[i, j].yaxis.set_label_coords(-0.5, 1.75)
        plt.tight_layout(pad=0.1, w_pad=-0.25, h_pad=-1.45)

        for i in range(numEle):
            for j in range(numAzi):
                ax_row2[i, j].set_ylim([0, yMax * 1.1])

        # saves plot as pdf
        plt.savefig(f'{unit}.pdf')
        plt.close('all')
        continue

    # calculate average XY RF position for units from session, excluding
    # units with values more than 1.5 SD away from the mean
    unitsXY = np.array(unitsXY)

    plt.close('all')
    np.save('unitsRFLocMat', allTuningMat)
    np.save('eleLabels', eleLabel)
    np.save('aziLabels', aziLabel)
    os.chdir('../../../Python Code')

print(time.time()-t0)
