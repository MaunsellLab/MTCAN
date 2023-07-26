"""
directionTuning generates a polar plot and histograms for each unit.
The script will save the plots for each unit as a PDF in a folder specific
to the day's dataset. This is the version for 12 directions tested

Chery June - 2023
"""

from usefulFns import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

# Start Here:
# Load relevant file here with pyMat reader
monkeyName = 'Meetz'
seshDate = '230720'
fileName = f'{monkeyName}_{seshDate}_GRF2_Spikes.mat'
allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

# create folder and change directory to save PDFs and np.array
if not os.path.exists('Direction Tuning'):
    os.makedirs('Direction Tuning')
os.chdir('Direction Tuning/')

# Tuning code
correctTrials = correctTrialsGRF(allTrials)
units = activeUnits('spikeData', allTrials)
unitsChannel = unitsInfo(units, correctTrials, allTrials)

# change stimDesc field to be a list of dictionaries
for corrTrial in correctTrials:
    currTrial = allTrials[corrTrial]
    nStim = len(currTrial['stimDesc']['data']['stimType'])
    currTrial['stimDesc']['data'] = [{k: v[i] for k, v in currTrial['stimDesc']['data'].items()}
                                     for i in range(nStim)]

frameRateHz = header['displayCalibration']['data']['frameRateHz']
numDir = header['map0Settings']['data']['directionDeg']['n']
interstimDurMS = header['mapInterstimDurationMS']['data']
histPrePostMS = 100
sponWindowMS = 100
allTuningMat = np.zeros((len(units), numDir))
allSEMMat = np.zeros((len(units), numDir))
# numBlocks = header['mappingBlockStatus']['data']['blockLimit']
numBlocks = allTrials[-1]['mappingBlockStatus']['data']['blocksDone']

# assert frame consistency during stimulus duration
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
    print('stimulus frame duration not consistent for mapping stimuli')
else:
    trueStimDurMS = np.int32(np.around(1000 / frameRateHz * stimDurFrame[0]))

# For each unit insert spike counts into matrix during valid stim presentation
for uCount, unit in enumerate(units):
    stimCount = np.zeros((1, numDir))
    spikeCountMat = np.zeros((numBlocks + 1, numDir))
    spikeHists = np.zeros((numDir, trueStimDurMS + 2 * histPrePostMS))
    sponSpikesArr = []

    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        if 'numMap0Stim' in currTrial:
            map0StimLim = currTrial['numMap0Stim']['data']
            map0Count = 0
            stimDesc = currTrial['stimDesc']['data']
            spikeData = currTrial['spikeData']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            for stim in stimDesc:
                if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                    dirIndex = int(stim['directionIndex'])
                    stCount = int(stimCount[0][dirIndex])
                    stimOnTimeS = ((1000 / frameRateHz * stim['stimOnFrame'])
                                   / 1000) + stim1TimeS
                    stimOffTimeS = ((1000 / frameRateHz * stim['stimOffFrame'])
                                    / 1000) + stim1TimeS
                    stimCount[0][dirIndex] += 1
                    map0Count += 1
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)[0]
                        unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        # spike count during stim presentation
                        stimSpikes = np.where((unitTimeStamps >= stimOnTimeS) &
                                              (unitTimeStamps <= stimOffTimeS))
                        spikeCountMat[stCount, dirIndex] = len(stimSpikes[0])
                        # spontaneous spike count (100ms before stim on)
                        sponSpikes = np.where((unitTimeStamps >= (stimOnTimeS - (sponWindowMS / 1000))) &
                                              (unitTimeStamps <= stimOnTimeS))
                        sponSpikesArr.extend([len(sponSpikes[0])])

                        # histograms
                        histStimSpikes = histSpikes(stimOnTimeS, stimOffTimeS,
                                                    histPrePostMS, unitTimeStamps)
                        spikeHists[dirIndex, histStimSpikes] += 1

    # summary stats
    spikeCountMean = np.mean(spikeCountMat[:numBlocks, :], axis=0)
    spikeCountSD = np.std(spikeCountMat[:numBlocks, :], axis=0)
    spikeCountSEM = spikeCountSD / np.sqrt(numBlocks)
    sponSpikesArr = np.array(sponSpikesArr)
    sponSpikesMean = np.mean(sponSpikesArr)
    sponSpikesSEM = np.std(sponSpikesArr) / np.sqrt(len(sponSpikesArr))
    allTuningMat[uCount] = spikeCountMean * 1000 / trueStimDurMS
    allSEMMat[uCount] = spikeCountSEM * 1000 / trueStimDurMS

    # Gaussian Fit
    angleMat = np.arange(180, 900, 30)
    extTunMat = np.concatenate((spikeCountMean[6:], spikeCountMean,
                                spikeCountMean[:6]), axis=0)
    nMax = int(np.argwhere(spikeCountMean == np.max(spikeCountMean))[0]) + 6
    nX = angleMat[nMax - 6:nMax + 7]
    # nY = extTunMat[nMax-3:nMax+4]/max(extTunMat[nMax-3:nMax+3])
    nY = extTunMat[nMax - 6:nMax + 7]
    params = gaussFit(nX, nY)
    nXFull = np.linspace(nX[0], nX[-1], 1000)
    nYFull = gauss(nXFull, *params)
    fitMean = params[2] - 360  # fitted mean
    fitVar = params[3] ** 2  # fitted var

    # Figure
    date = header['date']
    totalSpikes = np.sum(spikeCountMat)
    totalSpikesSec = int(totalSpikes * 1000 / trueStimDurMS)

    fig = plt.figure()
    fig.set_size_inches(6, 8)
    text = fig.text(0.05, 0.80, f'Direction tuning for unit {unit}\n\
    {date}\n\
    - - - - - - - - - - - - - -\n\
    Stimulus Duration: {trueStimDurMS} ms\n\
    Number of Blocks: {numBlocks}\n\
    Interstimulus Duration: {interstimDurMS} ms\n\
    Total Spikes: {totalSpikes}\n\
    Channel: {unitsChannel[uCount]}', size=10, fontweight='bold')
    text.set_path_effects([path_effects.Normal()])
    # Channel: {unitsChannel[uCount]}

    # Polar
    ax_row1 = plt.subplot2grid((12, 6), (0, 3), colspan=3, rowspan=4, polar=True)
    theta = np.radians(np.arange(0, 390, 360/numDir))
    sponTheta = np.radians(np.arange(0, 360, 360 / 100))
    sponTheta = np.append(sponTheta, sponTheta[0])
    r = (np.append(spikeCountMean, spikeCountMean[0])) * 1000 / trueStimDurMS
    err = (np.append(spikeCountSEM, spikeCountSEM[0])) * 1000 / trueStimDurMS
    spon = np.array([sponSpikesMean] * len(sponTheta))
    ax_row1.plot(theta, r, markersize=2)
    ax_row1.errorbar(theta, r, yerr=err, fmt='o', ecolor='black',
                     color='black', markersize=2)
    ax_row1.plot(sponTheta, spon * 1000 / sponWindowMS, linestyle='--')
    ax_row1.plot(np.radians(nXFull % 360), nYFull * 1000 / trueStimDurMS)
    ax_row1.set_theta_zero_location("W")
    ax_row1.set_title('Direction Tuning Polar Plot', fontsize=8)

    # Hists
    spikeHistsRS = np.reshape(spikeHists, (trueStimDurMS + 2 * histPrePostMS, 4, 3))
    stimCountRS = np.reshape(stimCount, (4, 3))
    titleArr = np.reshape(np.arange(0, 360, 30), (4, 3))

    ax_row2 = []
    for i in range(4, 12, 2):
        ax = []
        for j in np.arange(0, 6, 2):
            ax.append(plt.subplot2grid((12, 6), (i, j), colspan=2, rowspan=2))
        ax_row2.append(np.array(ax))
    ax_row2 = np.array(ax_row2)  # 4x3

    yMax = 0
    plotCount = 0
    for i in range(4):
        for j in range(3):
            spikeHist = spikeHists[plotCount, :] * 1000 / numBlocks
            plotCount += 1
            gaussSmooth = gaussian_filter1d(spikeHist, 5)
            if max(gaussSmooth) > yMax:
                yMax = max(gaussSmooth)
            ax_row2[i, j].plot(gaussSmooth)
            histTitle = titleArr[i][j]
            ax_row2[i, j].set_title(f"{histTitle}˚", fontsize=7)
            ax_row2[i, j].set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                                      2 * histPrePostMS + trueStimDurMS])
            ax_row2[i, j].set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                                           trueStimDurMS + histPrePostMS], fontsize=5)
            ax_row2[i, j].axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                                  color='grey', alpha=0.2)
            ax_row2[i, j].axhline(y=sponSpikesMean * 1000 / sponWindowMS,
                                  linestyle='--', color='grey')
            if i == 1 and j == 0:
                ax_row2[i, j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i, j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i, j].yaxis.set_label_coords(-0.2, 0.3)
    plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)

    for i in range(4):
        for j in range(3):
            ax_row2[i, j].set_ylim([0, yMax*1.1])

    # saves plot as pdf
    plt.savefig(f'{unit}.pdf')
    plt.close('all')
    continue

plt.close('all')
np.save('unitsDirTuningMat', allTuningMat)
np.save('unitsDirTuningSEM', allSEMMat)
