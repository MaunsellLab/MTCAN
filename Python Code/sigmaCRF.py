"""
GRF4 (sigma/contrast) main analysis script. This script will extract
the Contrast Response Function for a unit at differently sized sigmas. The script
will extract spike counts during valid trials from active units in a session and put them
in the appropriate matrix. From there we will plot the CRF for ~4 different sizes
of a Gabor at the center of the RF.

Chery - July 2024
"""

from usefulFns import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
import warnings

# Start Here:
# Load relevant file here with pyMat reader
monkeyName = 'Akshan'
seshDate = '240701'
unitOI = 1
fileName = f'{monkeyName}_{seshDate}_GRF4_Spikes.mat'
allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

# create folder and change directory to save PDFs and np.array
if not os.path.exists('SigmaCRF Tuning'):
    os.makedirs('SigmaCRF Tuning')
os.chdir('SigmaCRF Tuning/')

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
numContrasts = header['map0Settings']['data']['contrastPC']['n']
numSigmas = header['map0Settings']['data']['sigmaDeg']['n']
interstimDurMS = header['mapInterstimDurationMS']['data']
sponWindowMS = 100
onLatency = 50 / 1000  # time in MS for counting window latency after stim on
offLatency = 50 / 1000  # time in MS for counting window latency after stim off
histPrePostMS = 100  # 100ms window pre/post stimulus on/off
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

# insert spikes from valid stimulus presentations into spike count matrices
for uCount, unit in enumerate(units):
    if unit == unitOI:
        spikeCountMat = np.zeros((numBlocks+1, numSigmas*numContrasts))
        spikeHists = np.zeros((numSigmas*numContrasts, trueStimDurMS + 2 * histPrePostMS))
        stimCount = np.zeros((numSigmas, numContrasts))
        sponRate = []
        stimCountIndex = np.arange(numSigmas * numContrasts)
        stimCountIndex = stimCountIndex.reshape(numSigmas, numContrasts)

        for corrTrial in correctTrials:
            currTrial = allTrials[corrTrial]
            if 'numMap0Stim' in currTrial:
                map0StimLim = currTrial['numMap0Stim']['data']
                map0Count = 0
                stimDesc = currTrial['stimDesc']['data']
                if 'spikeData' in currTrial:
                    spikeData = currTrial['spikeData']
                    stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
                    for stim in stimDesc:
                        if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                            sigmaIndex = int(stim['sigmaIndex'])
                            contrastIndex = int(stim['contrastIndex'])
                            stCount = int(stimCount[sigmaIndex][contrastIndex])
                            stimOnTimeS = ((1000 / frameRateHz * stim['stimOnFrame'])
                                           / 1000) + stim1TimeS
                            stimOffTimeS = ((1000 / frameRateHz * stim['stimOffFrame'])
                                            / 1000) + stim1TimeS
                            stimCount[sigmaIndex][contrastIndex] += 1
                            stimIndex = stimCountIndex[sigmaIndex, contrastIndex]
                            map0Count += 1

                            if unit in np.array(spikeData['unit']):
                                spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                                unitIndex = np.where(spikeData['unit'] == unit)[0]
                                if type(spikeData['timeStamp']) == np.ndarray:
                                    unitTimeStamps = spikeData['timeStamp'][unitIndex]
                                else:
                                    unitTimeStamps = spikeData['timeStamp']
                                # spike count during stim presentation
                                stimSpikes = np.where((unitTimeStamps >= stimOnTimeS) &
                                                      (unitTimeStamps <= stimOffTimeS))
                                spikeCountMat[stCount, stimIndex] = len(stimSpikes[0])
                                # spontaneous spike count (100ms before stim on)
                                sponSpikes = np.where((unitTimeStamps >= (stimOnTimeS - (sponWindowMS / 1000))) &
                                                      (unitTimeStamps <= stimOnTimeS))
                                sponRate.append(len(sponSpikes[0]))

                                # histograms
                                histStimSpikes = histSpikes(stimOnTimeS, stimOffTimeS,
                                                            histPrePostMS, unitTimeStamps)
                                spikeHists[stimIndex, histStimSpikes] += 1

# mean, SEM, and reshaping of spikeCount matrices
meanSpike = np.mean(spikeCountMat[:numBlocks, :], axis=0)
spikeCountSD = np.std(spikeCountMat[:numBlocks, :], axis=0)
spikeCountSEM = spikeCountSD/np.sqrt(numBlocks)
meanSpikeReshaped = meanSpike.reshape(numSigmas, numContrasts) * 1000/trueStimDurMS
SEMReshaped = spikeCountSEM.reshape(numSigmas, numContrasts) * 1000/trueStimDurMS

# plot CRF for 4 differently sized Gabors
unit = unitOI
baselineResp = np.mean(sponRate) * 1000/trueStimDurMS
warnings.simplefilter('ignore', OptimizeWarning)
contrasts = [0, 3, 6.25, 12.5, 25, 50, 100]
sigmas = [0.1, 0.3333, 0.5667, 0.8]

fig, ax = plt.subplots(2, 2, figsize=(10, 8))
ax = ax.flatten()
for i in range(4):
    response = meanSpikeReshaped[i]
    response = np.insert(response, 0, baselineResp)
    try:
        ax[i].scatter(contrasts, response, color='green')
        initialGuess = [baselineResp, max(response), np.median(contrasts), 2.0]
        pOpt, pCov = curve_fit(contrastFn, contrasts, response,
                               bounds=([baselineResp, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
        xFit = np.logspace(-1, 2, 100)
        yFit = contrastFn(xFit, *pOpt)
        lower, upper = confidenceIntervalCRF(pOpt, pCov, xFit)
        ax[i].plot(xFit, yFit, color='green', alpha=0.5, label=f'{pOpt[2]:.2f}')
    except (RuntimeError, ValueError) as e:
        ax[i].scatter(contrasts, response, color='green')
    ax[i].set_xscale('symlog', linthresh=0.1)
    ax[i].set_xlabel('Contrast (%)')
    ax[i].set_ylabel('Spikes/s')
    ax[i].set_title(f'Contrast Response Function with sigma {sigmas[i]}')
    ax[i].legend()

plt.tight_layout()

# plt.show()

plt.savefig(f'{unit}.pdf')
plt.close('all')
