"""
MTSigma main analysis script. This script will extract spike counts during
valid trials from units active in a session and put them in the appropriate
matrix. From there, this script plots the contrast response function for pref/non-pref
at two points in the RF (center vs peri).

Chery - May 2024
"""

# Imports
from usefulFns import *
import numpy as np
import psignifit as ps
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
import warnings

# fileList = ['240603', '240606', '240610', '240701']
# unitList = ['240603_167', '240606_176', '240610_169', '240701_1']


# Start Here:
# Load relevant file here with pyMat reader
monkeyName = 'Akshan'
seshDate = '240701'
fileName = f'{monkeyName}_{seshDate}_MTSIG_Spikes.mat'
allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

# create folder and change directory to save PDFs and np.array
if not os.path.exists('CRF'):
    os.makedirs('CRF')
os.chdir('CRF/')

# list of indices of correctTrials (non-instruct, valid trialCertify)
corrTrials = correctTrialsMTX(allTrials)

# generate list of unique active units, and their channel
units = activeUnits('spikeData', allTrials)
unitCluster = allTrials[corrTrials[0]]['spikeTempInfo']['cgs']
unitsChannel = unitsInfo(units, corrTrials, allTrials)

# change stimDesc to be list of dictionaries
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    nStim = len(currTrial['stimDesc']['data']['locType'])
    currTrial['stimDesc']['data'] = [{k: v[i] for k, v in currTrial['stimDesc']['data'].items()}
                                     for i in range(nStim)]

# assert: are there correct trials without spikeData
noSpikeData = []
for trialCount, currTrial in enumerate(allTrials):
    trial = currTrial['trial']['data']
    extendedEOT = currTrial['extendedEOT']['data']
    if extendedEOT == 0 and trial['instructTrial'] != 1:
        if 'spikeData' not in currTrial:
            noSpikeData.append(trialCount)

# assert: frame consistency during stimulus duration
frameRateHz = header['frameRateHz']['data']
stimDurFrame = []
trueStimDurMS = 250
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    stimDesc = currTrial['stimDesc']['data']
    for stim in stimDesc:
        if stim['locType'] == 0:
            frameDiff = stim['stimOffFrame'] - stim['stimOnFrame']
            stimDurFrame.append(frameDiff)

if len(set(stimDurFrame)) != 1:
    print('stimulus frame duration not consistent for mapping stimuli')
else:
    trueStimDurMS = np.int32(np.around(1000 / frameRateHz * stimDurFrame[0]))

# initialize lists/arrays/dataframes for counting spikeCounts and for analysis
blocksDone = allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone'][0]
numContrasts = header['blockStatus']['data']['numContrasts'][0]
contrasts = header['blockStatus']['data']['contrasts'][0] * 100
contrasts = np.insert(contrasts, 0, 0)
spikeCountMat = np.zeros((len(units), blocksDone+1, numContrasts*4))
stimCount = np.zeros((2, 2, numContrasts), dtype=int)
sponRate = np.zeros((len(units), numContrasts*4*(blocksDone+1)))
sponIndex = 0
stimCountIndex = np.arange(numContrasts*4)
stimCountIndex = stimCountIndex.reshape(2, 2, numContrasts)
onLatency = 50 / 1000  # time in MS for counting window latency after stim on
offLatency = 50 / 1000  # time in MS for counting window latency after stim off
histPrePostMS = 100  # 100ms window pre/post stimulus on/off
spikeHists = np.zeros((len(units), numContrasts*4, trueStimDurMS + (2*histPrePostMS+1)))

# insert spikes from valid stimulus presentations into spike count matrices
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    if 'spikeData' in currTrial:
        stimDesc = currTrial['stimDesc']['data']
        stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
        for stim in stimDesc:
            if stim['stimType'] == 1 and stim['locType'] != 2:
                stimOnTimeS = ((1000 / frameRateHz * stim['stimOnFrame'])
                               / 1000) + stim1TimeS
                stimOffTimeS = ((1000 / frameRateHz * stim['stimOffFrame'])
                                / 1000) + stim1TimeS
                stimLocation = stim['locType']
                stimDirection = stim['dirType']
                stimContrast = stim['contrastIndex']
                stCount = int(stimCount[stimLocation, stimDirection, stimContrast])
                stimCount[stimLocation, stimDirection, stimContrast] += 1
                stimIndex = stimCountIndex[stimLocation, stimDirection, stimContrast]

                for unitCount, unit in enumerate(units):
                    if unit in currTrial['spikeData']['unit']:
                        unitIndex = np.where(currTrial['spikeData']['unit'] == unit)[0]
                        unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= (stimOnTimeS + onLatency)) &
                                              (unitTimeStamps <= (stimOffTimeS + offLatency)))[0]
                        spikeCountMat[unitCount][stCount][stimIndex] \
                            = len(stimSpikes)
                        sponRate[unitCount][sponIndex] = len(np.where((unitTimeStamps >= (stimOnTimeS - (100/1000))) &
                                                                      (unitTimeStamps <= stimOnTimeS))[0])

                        # PSTHs
                        stimOnPreSNEV = stimOnTimeS - (histPrePostMS / 1000)
                        stimOffPostSNEV = stimOffTimeS + (histPrePostMS / 1000)
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                                         & (unitTimeStamps <= stimOffPostSNEV)
                                                         )] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes * 1000)
                        spikeHists[unitCount, stimIndex, histStimSpikes] += 1

                sponIndex += 1

# mean, SEM, and reshaping of spikeCount matrices (if pref and non pref are different)
meanSpike = np.mean(spikeCountMat[:, :blocksDone, :], axis=1)
spikeCountSD = np.std(spikeCountMat[:, :blocksDone, :], axis=1)
spikeCountSEM = spikeCountSD/np.sqrt(blocksDone)
meanSpikeReshaped = np.zeros((len(units), 2, 2, numContrasts))
SEMReshaped = np.zeros((len(units), 2, 2, numContrasts))
for count, i in enumerate(meanSpikeReshaped):
    meanSpikeReshaped[count] = (meanSpike[count].reshape(2, 2, numContrasts) *
                                1000/trueStimDurMS)
    SEMReshaped[count] = (spikeCountSEM[count].reshape(2, 2, numContrasts) *
                          1000/trueStimDurMS)

# plot CRF
unit = 1
unitID = np.where(units == unit)[0][0]
baselineResp = np.mean(sponRate[unitID][:numContrasts*4*blocksDone]) * 1000/trueStimDurMS

colorID = 0
plotColors = ['green', 'red', 'green', 'red']
plotAlphas = [1, 1, 0.5, 0.5]
plt.figure(figsize=(12, 8))
warnings.simplefilter('ignore', OptimizeWarning)

for i in range(2):
    for j in range(2):
        response = meanSpikeReshaped[unitID][i][j]
        response = np.insert(response, 0, baselineResp)
        try:
            plt.scatter(contrasts, response, color=plotColors[colorID],
                        alpha=plotAlphas[colorID])
            initialGuess = [baselineResp, max(response), np.median(contrasts), 2.0]
            pOpt, pCov = curve_fit(contrastFn, contrasts, response,
                                bounds=([baselineResp, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
            xFit = np.logspace(-1, 2, 100)
            yFit = contrastFn(xFit, *pOpt)
            lower, upper = confidenceIntervalCRF(pOpt, pCov, xFit)
            plt.plot(xFit, yFit, color=plotColors[colorID],
                     alpha=plotAlphas[colorID], label=f'{pOpt[2]:.2f}')
            if j == 2:
                plt.fill_between(xFit, lower, upper, color=plotColors[colorID], alpha=0.2)
        except (RuntimeError, ValueError) as e:
            plt.scatter(contrasts, response, color=plotColors[colorID], alpha=plotAlphas[colorID])
        colorID += 1

plt.xscale('symlog', linthresh=0.1)
plt.xlabel('Contrast (%)')
plt.ylabel('Spikes/s')
plt.title(f'Contrast Response Function, unit {unit}')
plt.legend()
sns.despine(offset=5)

# plt.show()

plt.savefig(f'{unit}.pdf')
plt.close('all')


