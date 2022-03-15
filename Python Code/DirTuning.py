

from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import seaborn as sns

allTrials, header = loadMatFile73('testing_220310_Heatmap_GRF_Spikes.mat')

# for testing purposes, to make unit field similar to real data
for currTrial in allTrials:
    if 'spikeData' in currTrial:
        currTrial['spikeData']['unit'] = currTrial['spikeData']['unit'].tolist()
        for i in range(0,len(currTrial['spikeData']['channel'])):
            a = str(int(currTrial['spikeData']['channel'][i]))
            b = str(int(currTrial['spikeData']['unit'][i]))
            c = a + '_' + b
            currTrial['spikeData']['unit'][i] = c
        currTrial['spikeData']['unit'] = np.array(currTrial['spikeData']['unit'])

def activeUnits(spikeData):

    '''
    function returns the active units across trials for a session as a list

    Inputs: unitData (str): spikeData
    Outputs: units (list): active units for a sessioon

    '''
    units = []
    for currTrial in allTrials:
        if spikeData in currTrial:
            uniqueUnits = np.unique(currTrial[spikeData]['unit'])
            for unique in uniqueUnits:
                if unique not in units:
                    units.append(unique)
    
    return units

units = activeUnits('spikeData')

numDir = int(header['mapSettings']['data']['directionDeg']['n'].tolist())
stimDurMS = int(header['mapStimDurationMS']['data'].tolist())
histPrePostMS = 50

stimCount = np.zeros((1,numDir))
spikeCountMat = np.zeros((50,1, numDir))
spikeCountMat[:,:,:] = np.nan
spikeHists = np.zeros((stimDurMS + 2*histPrePostMS, 1, numDir))

for currTrial in allTrials:
    trial = currTrial['trial']['data']
    trialEnd = currTrial['trialEnd']['data']
    if trial['instructTrial'] != 1 and trialEnd == 0:
        trialStartMS = currTrial['trialStart']['timeMS']
        trialStartSNEV = np.around(currTrial['taskEvents']['trialStart']['timeS'], 3)
        stimDesc = currTrial['stimDesc']
        spikeData = currTrial['spikeData']
        for sCount, stim in enumerate(stimDesc['data']):
            if stim['stimType'] == 2:
                break
            if stim['gaborIndex'] == 1:
                dirIndex = int(stim['directionIndex'])
                stCount = int(stimCount[0][dirIndex])
                stimOnTimeMS = stimDesc['timeMS'][sCount]
                stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                stimOnSNEV = trialStartSNEV + stimDiffS
                if unit in spikeData['unit']:
                    spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                    unitIndex = np.where(spikeData['unit'] == unit)
                    unitTimeStamps = spikeData['timeStamp'][unitIndex]
                    stimSpikes = np.where((unitTimeStamps >= stimOnSNEV) & 
                                    (unitTimeStamps <= stimOnSNEV + stimDurMS/1000))


                    spikeCountMat[stCount][0][dirIndex] = len(stimSpikes[0])
                    stimCount[0][dirIndex] += 1
                    
                    #histograms
                    histSpikes = np.arange(stimOnSNEV - 0.050, stimOnSNEV + \
                                          (stimDurMS+49)/1000, 0.001)
                    for histCount, i in enumerate(range(len(histSpikes))):
                        if np.around(histSpikes[i],3) in unitTimeStamps:
                            spikeHists[histCount][0][directionIndex] += 1


spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)