import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import scipy.io as sp
import math
import os
from collections import defaultdict
from usefulFns import *


'''
# %% for conda/jupyter notebook
    
Fake Data Gen Mat 7.3, with fakeData and spikeData (channel, unit, column) as fields 
'''

def insertStimSpikeData(units, index, stimOnTimeSNEV):
    '''
    this function will generate a poisson spike train and return the normalized
    response for the RF for a given stimulus configuration
    Inputs:
        x: neuron number
        index: the index of the corresponding stimulus configuration
        sigma: semisaturation constant from neuron's contrast response function
    Outputs:
        RFSpikes: normalized response 

    spikes_4rows = np.tile(spikes, (4,1))
    '''


    numNeurons = len(units)

    C0 = stimIndexDict[index][0]['contrast']
    C1 = stimIndexDict[index][1]['contrast']
    L0 = tcDict[1][stimIndexDict[index][0]['direction']]
    L1 = tcDict[1][stimIndexDict[index][1]['direction']]
    sigma = 0.1
    expectedNormSpikeRate = int(((C0*L0) + (C1*L1))/(C0 + C1 + sigma))
    stimDur = 493
    popMean = expectedNormSpikeRate/(1000/stimDur)
    popSTD = math.sqrt(popMean)

    spikes = popMean + np.random.rand(1,numNeurons)     
    # spikes = popMean + popSTD*np.random.rand(1,numNeurons)      
    R = np.zeros((numNeurons,numNeurons))
    for neuronI in range(numNeurons):
        for neuronJ in range(numNeurons):
            if neuronI == neuronJ:
                R[neuronI,neuronJ] = 1
            else:
                R[neuronI, neuronJ] = 0.1

    L = np.linalg.cholesky(R)
    spikes = np.matmul(spikes,L)   
    spikes = np.around(spikes) 

    for count, i in enumerate(spikes[0]):
        if i != 0:
            unit = units[count]
            channelIdentity = int(unit[0:unit.find('_')])
            channel = np.array([channelIdentity] * stimDur)
            spikeTimeMS = (np.sort(np.random.choice(np.arange(stimDur), int(i),
            replace = False)))/1000
            currTrial['spikeData']['timeStamp'] = np.append(currTrial['spikeData'] \
            ['timeStamp'], stimOnTimeSNEV + spikeTimeMS, 0)
            currTrial['spikeData']['unit'] = np.append(currTrial['spikeData'] \
            ['unit'], [unit] * len(spikeTimeMS), 0)
            currTrial['spikeData']['channel'] = np.append(currTrial['spikeData'] \
            ['channel'], [channelIdentity] * len(spikeTimeMS), 0)


def randTuningCurve(numNeurons):
    '''
    Functon will generate random tuning cruves for x number of neurons.
    Function will also iterate through tuning curves to create dictionary of 
    responses to each direction for each neuron 

    Inputs:
        numNueurons (int): number of neurons 
    Outputs:
        tuningMat (2D array): matrix of tuning curve values for each
                              neuron
        tcDictionary (dictionary): maps each neurons response for a direction onto
                                   a dictionary
    '''
    tuningMat = np.zeros((numNeurons + 1, 6))
    tuningMat[0] = np.arange(0,360,60)
    tcDictionary = {}

    for i in range(1, tuningMat.shape[0]):
        np.random.seed(2)
        amp = np.random.randint(15,30)
        y_translate = np.random.randint(30,50)
        x_translate = np.random.randint(60,120)
        tuningMat[i,:] = (amp * np.sin((tuningMat[0,:] * (np.pi / 180)) + x_translate)) + y_translate

    for i, neuron in enumerate(tuningMat[1:,:]):
        tcDictionary[i+1] = {}
        for j, dirResp in enumerate(neuron):
            tcDictionary[i+1][tuningMat[0][j]] = dirResp

    return tuningMat, tcDictionary


# start here 
allTrials, header = loadMatFile73('testing_220304_MTN_Spikes.mat')

# generates a dictionary of stim Index and corresponding directions/contrasts
stimIndexDict = {}
for currTrial in allTrials:
    extendedEOT = currTrial['extendedEOT']['data']
    trial = currTrial['trial']['data']
    if extendedEOT == 0 and trial['instructTrial'] != 1:
        stimDesc = currTrial['stimDesc']['data']
        for stim in stimDesc:
            if stim['stimLoc'] != 2:
                index = int(stim['stimIndex'].tolist())
                if index not in stimIndexDict:
                    stimIndexDict[index] = {}
                    if int(stim['stimLoc'].tolist()) not in stimIndexDict[index]:
                        stimIndexDict[index][int(stim['stimLoc'].tolist())] = \
                        {'direction': int(stim['directionDeg'].tolist()),
                         'contrast': stim['contrast'].tolist()}
                else:
                    if int(stim['stimLoc'].tolist()) not in stimIndexDict[index]:
                        stimIndexDict[index][int(stim['stimLoc'].tolist())] = \
                        {'direction': int(stim['directionDeg'].tolist()),
                         'contrast': stim['contrast'].tolist()}

# trials that should have spikeData but don't
noSpikeData = []
for count, currTrial in enumerate(allTrials):
    extendedEOT = currTrial['extendedEOT']['data']
    trial = currTrial['trial']['data']
    if trial['instructTrial'] != 1 and extendedEOT == 0:
        if 'spikeData' not in currTrial:
            noSpikeData.append(count)

# for testing only, Pre-processing, to generate spikeData field similar 
# to actual spikeData from .nev (will have units as function of channel,
# 1_1, 2_1 etc.)
for currTrial in allTrials:
    if 'spikeData' in currTrial:
        currTrial['spikeData']['unit'] = currTrial['spikeData']['unit'].tolist()
        for i in range(0,len(currTrial['spikeData']['channel'])):
            a = str(int(currTrial['spikeData']['channel'][i]))
            b = str(int(currTrial['spikeData']['unit'][i]))
            c = a + '_' + b
            currTrial['spikeData']['unit'][i] = c
        currTrial['spikeData']['unit'] = np.array(currTrial['spikeData']['unit'])

# code to insert fake spikes during stimulus presentations 
units = activeUnits('spikeData', allTrials)
tuningMat, tcDict = randTuningCurve(len(units))

for currTrial in allTrials:
    if 'spikeData' in currTrial:
        currTrial['spikeData']['channel'] = []
        currTrial['spikeData']['unit'] = []
        currTrial['spikeData']['timeStamp'] = []
        stimDesc = currTrial['stimDesc']['data']
        for count, stim in enumerate(stimDesc):
            if stim['stimLoc'] == 0 and stim['listType'] == 1:
                index = int(stim['stimIndex'].tolist())
                stimOnTimeMS = currTrial['stimDesc']['timeMS'][count] - \
                currTrial['trialStart']['timeMS']
                stimOnTimeSNEV = round(currTrial['taskEvents']['trialStart']\
                ['timeS'].tolist() + (stimOnTimeMS/1000), 3)
                insertStimSpikeData(units, index, stimOnTimeSNEV)

# Fake data for GaborRF heatmap generation 
# Values for fakespikes
fakeAzi = np.random.randint(numAzi)
fakeEle = np.random.randint(numEle)
fakeSpikes = 50 #spikes/sec
baseFR = 10 #spikes/sec
dt = 1/1000

for currTrial in allTrials:
    trial = currTrial['trial']['data']
    trialEnd = currTrial['trialEnd']['data']
    if trial['instructTrial'] != 1 and trialEnd == 0:
        trialStartMS = currTrial['trialStart']['timeMS']
        trialStartSNEV = currTrial['taskEvents']['trialStart']['timeS']
        trialEndSNEV = currTrial['taskEvents']['trialEnd']['timeS']
        trialLen = trialEndSNEV - trialStartSNEV
        stimDesc = currTrial['stimDesc']
        spikeData = currTrial['spikeData']
        spikeData['channel'] = []
        spikeData['unit'] = []
        spikeData['timeStamp'] = []
        for unit in units:
            channelIdentity = int(unit[0:unit.find('_')])
            spikes = np.random.poisson(trialLen*baseFR)
            spikeTimeS = (np.sort(np.random.choice(np.arange(trialLen * 1000),\
                           spikes, replace = False)))/1000
            spikeData['timeStamp'] = np.append(spikeData['timeStamp'], \
                                     trialStartSNEV + spikeTimeS, 0)
            spikeData['unit'] = np.append(spikeData['unit'], \
                                [unit] * len(spikeTimeS), 0)
            spikeData['channel'] = np.append(spikeData['channel'], \
                                   [channelIdentity] * len(spikeTimeS), 0)
        for count, stim in enumerate(stimDesc['data']):
            if stim['stimType'] == 2:
                break
            if stim['gaborIndex'] == 1:
                aziIndex = int(stim['azimuthIndex'])
                eleIndex = int(stim['elevationIndex'])
                stimOnTimeMS = stimDesc['timeMS'][count]
                stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                stimOnSNEV = trialStartSNEV + stimDiffS
                if aziIndex == fakeAzi and eleIndex == fakeEle:
                    for unit in units:
                        channelIdentity = int(unit[0:unit.find('_')])
                        spikeTimeS = (np.sort(np.random.choice(np.arange(stimDurMS),\
                                    int(np.round(fakeSpikes/(1000/stimDurMS))),\
                                    replace = False)))/1000
                        spikeData['timeStamp'] = np.append(spikeData['timeStamp'],\
                                                stimOnSNEV + spikeTimeS, 0)
                        spikeData['unit'] = np.append(spikeData['unit'], \
                                            [unit] * len(spikeTimeS), 0)
                        spikeData['channel'] = np.append(spikeData['channel'], \
                                               [channelIdentity] * len(spikeTimeS), 0)


# fake spikes for dir tuning
tuningMat, tcDict = randTuningCurve(len(units))
baseFR = 10 #spikes/sec
dt = 1/1000

for currTrial in allTrials:
    trial = currTrial['trial']['data']
    trialEnd = currTrial['trialEnd']['data']
    if trial['instructTrial'] != 1 and trialEnd == 0:
        trialStartMS = currTrial['trialStart']['timeMS']
        trialStartSNEV = currTrial['taskEvents']['trialStart']['timeS']
        trialEndSNEV = currTrial['taskEvents']['trialEnd']['timeS']
        trialLen = trialEndSNEV - trialStartSNEV
        stimDesc = currTrial['stimDesc']
        spikeData = currTrial['spikeData']
        spikeData['channel'] = []
        spikeData['unit'] = []
        spikeData['timeStamp'] = []
        for unit in units:
            channelIdentity = int(unit[0:unit.find('_')])
            spikes = np.random.poisson(trialLen*baseFR)
            spikeTimeS = (np.sort(np.random.choice(np.arange(trialLen * 1000),\
                           spikes, replace = False)))/1000
            spikeData['timeStamp'] = np.append(spikeData['timeStamp'], \
                                     trialStartSNEV + spikeTimeS, 0)
            spikeData['unit'] = np.append(spikeData['unit'], \
                                [unit] * len(spikeTimeS), 0)
            spikeData['channel'] = np.append(spikeData['channel'], \
                                   [channelIdentity] * len(spikeTimeS), 0)
        for count, stim in enumerate(stimDesc['data']):
            if stim['stimType'] == 2:
                break
            if stim['gaborIndex'] == 1:
                dirIndex = int(stim['directionIndex'])
                stimOnTimeMS = stimDesc['timeMS'][count]
                stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                stimOnSNEV = trialStartSNEV + stimDiffS
                for unitCount, unit in enumerate(units):
                    expectedSpikes = np.random.poisson(tuningMat[unitCount+1][dirIndex]*stimDurMS/1000)
                    channelIdentity = int(unit[0:unit.find('_')])
                    spikeTimeS = (np.sort(np.random.choice(np.arange(stimDurMS),\
                    expectedSpikes, replace=False)))/1000
                    spikeData['timeStamp'] = np.append(spikeData['timeStamp'],\
                                            stimOnSNEV + spikeTimeS, 0)
                    spikeData['unit'] = np.append(spikeData['unit'], \
                                        [unit] * len(spikeTimeS), 0)
                    spikeData['channel'] = np.append(spikeData['channel'], \
                                            [channelIdentity] * len(spikeTimeS), 0)


# fake spikes for speed tuning
speedMat = [40,80,40,20] #spikes/sec at each speed index
baseFR = 10 #spikes/sec
dt = 1/1000

for currTrial in allTrials:
    trial = currTrial['trial']['data']
    trialEnd = currTrial['trialEnd']['data']
    if trial['instructTrial'] != 1 and trialEnd == 0:
        trialStartMS = currTrial['trialStart']['timeMS']
        trialStartSNEV = currTrial['taskEvents']['trialStart']['timeS']
        trialEndSNEV = currTrial['taskEvents']['trialEnd']['timeS']
        trialLen = trialEndSNEV - trialStartSNEV
        stimDesc = currTrial['stimDesc']
        spikeData = currTrial['spikeData']
        spikeData['channel'] = []
        spikeData['unit'] = []
        spikeData['timeStamp'] = []
        for unit in units:
            channelIdentity = int(unit[0:unit.find('_')])
            spikes = np.random.poisson(trialLen*baseFR)
            spikeTimeS = (np.sort(np.random.choice(np.arange(trialLen * 1000),\
                           spikes, replace = False)))/1000
            spikeData['timeStamp'] = np.append(spikeData['timeStamp'], \
                                     trialStartSNEV + spikeTimeS, 0)
            spikeData['unit'] = np.append(spikeData['unit'], \
                                [unit] * len(spikeTimeS), 0)
            spikeData['channel'] = np.append(spikeData['channel'], \
                                   [channelIdentity] * len(spikeTimeS), 0)
        for count, stim in enumerate(stimDesc['data']):
            if stim['stimType'] == 2:
                break
            if stim['gaborIndex'] == 1:
                speedIndex = int(stim['temporalFreqIndex'])
                stimOnTimeMS = stimDesc['timeMS'][count]
                stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                stimOnSNEV = trialStartSNEV + stimDiffS
                for unitCount, unit in enumerate(units):
                    expectedSpikes = np.random.poisson(speedMat[speedIndex]*stimDurMS/1000)
                    channelIdentity = int(unit[0:unit.find('_')])
                    spikeTimeS = (np.sort(np.random.choice(np.arange(stimDurMS),\
                    expectedSpikes, replace=False)))/1000
                    spikeData['timeStamp'] = np.append(spikeData['timeStamp'],\
                                            stimOnSNEV + spikeTimeS, 0)
                    spikeData['unit'] = np.append(spikeData['unit'], \
                                        [unit] * len(spikeTimeS), 0)
                    spikeData['channel'] = np.append(spikeData['channel'], \
                                            [channelIdentity] * len(spikeTimeS), 0)




for trialCount, corrTrial in enumerate(corrTrials):
    currTrial = allTrials[corrTrial]
    if 'spikeData' in currTrial:
        currTrial['spikeData']['unit'] = []
        currTrial['spikeData']['timeStamp'] = []
        stimDesc = currTrial['stimDesc']['data']
        stim1TimeS = currTrial['taskEvents']['stimulusOn'][0]['time'].tolist()
        for stim in stimDesc:
            if stim['stimLoc'] == 0 and stim['listType'] == 1:
                stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'].tolist())
                              /1000) + stim1TimeS
                stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'].tolist())
                               /1000) + stim1TimeS
                fakeSpikes = np.random.uniform(low=stimOnTimeS, high=stimOffTimeS, size=10).tolist()
                currTrial['spikeData']['timeStamp'].extend(fakeSpikes)
        currTrial['spikeData']['timeStamp'] = np.asarray(currTrial['spikeData']['timeStamp'])
        currTrial['spikeData']['unit'] = np.ones(len(currTrial['spikeData']['timeStamp']))
        currTrial['spikeData']['unit'] = np.asarray(currTrial['spikeData']['unit'])

            



'''
old verison pre v7.3
'''
stimIndexDict = {}
for currTrial in allTrialsData.item()[0]:
    extendedEOT = currTrial['extendedEOT'].item()['data'].item()
    trial = currTrial['trial'].item()['data'].item()
    if extendedEOT == 0 and trial['instructTrial'] != 1:
        stimDesc = currTrial['stimDesc'].item()['data'].item()
        for stim in stimDesc:
            if stim['stimLoc'] != 2:
                index = stim['stimIndex']
                if index not in stimIndexDict:
                    stimIndexDict[index] = {}
                else:
                    if stim['stimLoc'] not in stimIndexDict[index]:
                        stimIndexDict[index][stim['stimLoc']] = \
                        {'direction': stim['directionDeg'],
                         'contrast': stim['contrast']}

'''
backup
'''
for currTrial in allTrials:
    if 'spikeData' in currTrial:
        currTrial['fakeData'] = currTrial['spikeData'].copy()
        currTrial['fakeData']['unit'] = currTrial['fakeData']['unit'].tolist()
        for i in range(0,len(currTrial['spikeData']['channel'])):
            a = str(int(currTrial['fakeData']['channel'][i]))
            b = str(int(currTrial['fakeData']['unit'][i]))
            c = a + '_' + b
            currTrial['fakeData']['unit'][i] = c
        currTrial['fakeData']['unit'] = np.array(currTrial['fakeData']['unit'])              

'''

def poissonSpikeTrain(x, index):
    '''
    this function will generate a poisson spike train and return the normalized
    response for the RF for a given stimulus configuration
    Inputs:
        x: neuron number
        index: the index of the corresponding stimulus configuration
        sigma: semisaturation constant from neuron's contrast response function
    Outputs:
        RFSpikes: normalized response 
    '''
    fRateLocO = tcDict[x][stimIndexDict[index][0]['direction']]
    fRateLoc1 = tcDict[x][stimIndexDict[index][1]['direction']]
    C0 = stimIndexDict[index][0]['contrast']
    C1 = stimIndexDict[index][1]['contrast']
    spikeTrain0 = []
    spikeTrain1 = []
    sigma = 0.3

    dt = 1/1000
    for i in range(500):
        if np.random.uniform(0,1) < fRateLoc0 * dt:
            spikeTrain0.append(1)
        if np.random.uniform(0,1) < fRateLoc1 * dt:
            spikeTrain1.append(1)
    L0 = len(spikeTrain0)
    L1 = len(spikeTrain1)

    RFSpikes = ((C0*L0) + (C1*L1))/(C0+C1+sigma)

    return RFSpikes


allTrialsData, header = loadMatFile('Meetz_2022_0114_MTNAN3.mat')



Run these lines after a dictionary of the stimIndex with corresponding directions
and contrasts at each location has been created

numNeurons = 2
tuningCurves, tcDict = randTuningCurve(numNeurons)

# creates spike count matrix for multiple neurons, each row is a spike response
# to a presentation of the stimulus at a specific index 
spikeCountMat = np.zeros((numNeurons,30,169))
spikeCountMat[:,0,:] = np.arange(0,169)
spikeCountMat[:,1:,:] = np.nan
stimIndexCount = {}

def insertStimSpikeData(x, index, stimOnTimeSNEV):
    '''
    this function will generate a poisson spike train and return the normalized
    response for the RF for a given stimulus configuration
    Inputs:
        x: neuron number
        index: the index of the corresponding stimulus configuration
        sigma: semisaturation constant from neuron's contrast response function
    Outputs:
        RFSpikes: normalized response 
    '''
    
    unitIdentity = units[x]
    unitsField = np.array([unitIdentity] * 493)
    channelIdentity = int(unitIdentity[0:unitIdentity.find('_')])
    channel = np.array([channelIdentity] * 493)


    C0 = stimIndexDict[index][0]['contrast']
    C1 = stimIndexDict[index][1]['contrast']
    L0 = tcDict[x+1][stimIndexDict[index][0]['direction']]
    L1 = tcDict[x+1][stimIndexDict[index][1]['direction']]
    sigma = 0.1

    expectedNormSpikeRate = int(((C0*L0) + (C1*L1))/(C0 + C1 + sigma))
    dt = 1/1000
    baseSpikeTrain = [0] * 493
    for i in range(len(baseSpikeTrain)):
        if np.random.uniform(0,1) < expectedNormSpikeRate * dt:
            baseSpikeTrain[i] = 1
    spikeIndex = np.where(np.array(baseSpikeTrain) == 1)
     
    # add spikes to spikeData
    if len(spikeIndex[0]) != 0:
        spikeIndexS = spikeIndex[0]/1000 #[0] for tuple
        currTrial['spikeData']['timeStamp'] = np.append(currTrial \
        ['spikeData']['timeStamp'], stimOnTimeSNEV + spikeIndexS, 0)
        currTrial['spikeData']['unit'] = np.append(currTrial['spikeData'] \
        ['unit'], [unitIdentity] * len(spikeIndexS), 0)
        currTrial['spikeData']['channel'] = np.append(currTrial['spikeData'] \
        ['channel'], [channelIdentity] * len(spikeIndexS), 0)


# will print trial nunmbers that have targets that appear during RF stimulus presentation
for count, currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data'].item()
    trial = currTrial['trial'].item()['data'].item()
    if extendedEOT == 0 and trial['instructTrial'] != 1 and trial['catchTrial'] != 1:
        stimDesc = currTrial['stimDesc'].item()['data'].item() 
        for n,stim in enumerate(stimDesc):
            if stim['listType'] == 2:
                targetIndex = n
                break
        targetOnFrame = stimDesc['stimOnFrame'][targetIndex]
        for stim in stimDesc:
            if stim['stimLoc'] == 0 and stim['stimOffFrame'] > targetOnFrame:
                print(count, 'this trial has target appear before last RF stimuli is off')



'''
