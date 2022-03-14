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

def activeUnits(unitData):

    '''
    function returns the active units across trials for a session as a list

    Inputs: unitData (str): Are we using fakeData or spikeData
    Outputs: units (list): active units for a sessioon

    '''
    units = []
    for currTrial in allTrials:
        if unitData in currTrial:
            uniqueUnits = np.unique(currTrial[unitData]['unit'])
            for unique in uniqueUnits:
                if unique not in units:
                    units.append(unique)
    
    return units


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
units = activeUnits('spikeData')
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
57.085
58.212
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

# for generating underlying base firing rate
for currTrial in allTrials:
    if 'spikeData' in currTrial:
        trialStartS = currTrial['taskEvents']['trialStart']['timeS']
        trialEndS = currTrial['taskEvents']['trialEnd']['timeS']
        trialLen = int(round(trialEndS - trialStartS, 3) * 1000)
        # generate underlying firing rate 
        for count, unit in enumerate(units):
            channelIdentity = int(unit[0:unit.find('_')])
            baseFR = 5 # same firing rate for each neuron
            dt = 1/1000
            baseSpikeTrain = [0] * trialLen
            for i in range(trialLen):
                if np.random.uniform(0,1) < baseFR * dt:
                    baseSpikeTrain[i] = 1
            spikeIndex = np.where(np.array(baseSpikeTrain) == 1)
            spikeIndexS = spikeIndex[0]/1000 #[0] for tuple 
            if count == 0:
                currTrial['spikeData']['timeStamp'] = trialStartS + spikeIndexS
                currTrial['spikeData']['unit'] = [unit] * len(spikeIndexS)
                currTrial['spikeData']['channel'] = [channelIdentity] * len(spikeIndexS)
            else:
                currTrial['spikeData']['timeStamp'] = np.append(currTrial \
                ['spikeData']['timeStamp'], trialStartS + spikeIndexS, 0)
                currTrial['spikeData']['unit'] = np.append(currTrial['spikeData'] \
                ['unit'], [unit] * len(spikeIndexS), 0)
                currTrial['spikeData']['channel'] = np.append(currTrial['spikeData'] \
                ['channel'], [channelIdentity] * len(spikeIndexS), 0)
                # extend the currTrial['spikeData']['timeStamp'] fields
        stimDesc = currTrial['stimDesc']['data']
        for count, stim in enumerate(stimDesc):
            if stim['stimLoc'] == 0 and stim['listType'] == 1:
                index = int(stim['stimIndex'].tolist())
                stimOnTimeMS = currTrial['stimDesc']['timeMS'][count] - \
                currTrial['trialStart']['timeMS']
                stimOnTimeSNEV = round(currTrial['taskEvents']['trialStart']\
                ['timeS'].tolist() + (stimOnTimeMS/1000), 3)
                for j in range(0, len(units)):
                    insertStimSpikeData(j, index, stimOnTimeSNEV)

# code will add stimuli presentations to index matrix, excluding padding
# stimuli and when target appears before RF stimulus turns off
for n,currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data'].item()
    trial = currTrial['trial'].item()['data'].item()
    targetOnFrame = 20000
    if extendedEOT == 0 and trial['instructTrial'] != 1:
        stimDesc = currTrial['stimDesc'].item()['data'].item()
        for n,stim in enumerate(stimDesc):
            if stim['listType'] == 2:
                targetIndex = n
                break
        targetOnFrame = stimDesc['stimOnFrame'][targetIndex]
        stimCount = 0
        for stim in stimDesc:
            if stim['stimLoc'] == 0 and stim0Count == 0:
                stim0Count += 1
            elif stim['stimLoc'] == 0 and stim0Count != 0 and stim['stimOffFrame'] < targetOnFrame:
                index = stim['stimIndex']
                stimIndexCount[index] = stimIndexCount.get(index, 0) + 1

                for x in range(1,numNeurons+1):
                    spikeCountMat[x-1][stimIndexCount[index]][index] = poissonSpikeTrain(x, index)


meanSpike = np.nanmean(spikeCountMat[:,1:,:], axis = 1)

meanSpikeReshaped = np.zeros((len(units),1,169))
for count,i in enumerate(meanSpikeReshaped):
    i[:,0:6] = meanSpike[count][0:6]
    i[:,6:12] = meanSpike[count][36:42]
    i[:,12] = meanSpike[count][156]
    i[:,13:19] = meanSpike[count][6:12]
    i[:,19:25] = meanSpike[count][42:48]
    i[:,25] = meanSpike[count][157]
    i[:,26:32] = meanSpike[count][12:18]
    i[:,32:38] = meanSpike[count][48:54]
    i[:,38] = meanSpike[count][158]
    i[:,39:45] = meanSpike[count][18:24]
    i[:,45:51] = meanSpike[count][54:60]
    i[:,51] = meanSpike[count][159]
    i[:,52:58] = meanSpike[count][24:30]
    i[:,58:64] = meanSpike[count][60:66]
    i[:,64] = meanSpike[count][160]
    i[:,65:71] = meanSpike[count][30:36]
    i[:,71:77] = meanSpike[count][66:72]
    i[:,77] = meanSpike[count][161]
    i[:,78:84] = meanSpike[count][72:78]
    i[:,84:90] = meanSpike[count][108:114]
    i[:,90] = meanSpike[count][162]
    i[:,91:97] = meanSpike[count][78:84]
    i[:,97:103] = meanSpike[count][114:120]
    i[:,103] = meanSpike[count][163]
    i[:,104:110] = meanSpike[count][84:90]
    i[:,110:116] = meanSpike[count][120:126]
    i[:,116] = meanSpike[count][164]
    i[:,117:123] = meanSpike[count][90:96]
    i[:,123:129] = meanSpike[count][126:132]
    i[:,129] = meanSpike[count][165]
    i[:,130:136] = meanSpike[count][96:102]
    i[:,136:142] = meanSpike[count][132:138]
    i[:,142] = meanSpike[count][166]
    i[:,143:149] = meanSpike[count][102:108]
    i[:,149:155] = meanSpike[count][138:144]
    i[:,155] = meanSpike[count][167]
    i[:,156:168] = meanSpike[count][144:156]
    i[:,168] = meanSpike[count][168]

a = meanSpikeReshaped[0]
b = a.reshape(13,13)
plt.imshow(b, cmap='hot', interpolation='nearest')
    

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
