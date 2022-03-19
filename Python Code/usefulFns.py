import scipy.io as sp
import numpy as np
import matplotlib.pyplot as plt
import math
import os
from collections import defaultdict
import mat73

def loadMatFile73(fileName):
    '''
    Loads the given matfile and assigns variables to access trial data

    Inputs: matfile name, (str)
    Outputs: variables, (nd.array)
    '''
    
    allTrials = mat73.loadmat(f'../Matlab Data/{fileName}', use_attrdict = True)
    allTrialsData = allTrials.trials
    header = allTrials.header

    return allTrialsData, header

def loadMatFile(fileName):
    '''
    Loads the given matfile and assigns variables to access trial data

    Inputs: matfile name, (str)
    Outputs: variables, (nd.array)
    '''
    
    allTrials = sp.loadmat(f'../Matlab Data/{fileName}', squeeze_me = True)
    allTrialsData = allTrials['trials']
    header = allTrials['header']

    return allTrialsData, header

def correctTrialsGRF(allTrials):
    '''
    fn that will filter for correct non-instruct trials (valid) in GRF

    inputs: allTrials (list of trials)
    outputs: correctTrials (list): index of valid trials from allTrials
    '''
    correctTrials = []
    for trialCount, currTrial in enumerate(allTrials):
        trial = currTrial['trial']['data']
        trialEnd = currTrial['trialEnd']['data']
        if trial['instructTrial'] != 1 and trialEnd == 0:
            correctTrials.append(trialCount)
    
    return correctTrials

def correctTrialsNotInstruct(currTrial, catch=True):
    '''
    Function will return a bool if the current trial is correct and not an 
    instruction trial. Also has option to exclude catch trials 

    Inputs: trial (nd.array)
            catch (bool)
    Outputs: bool
    '''

    extendedEOT = currTrial['extendedEOT'].item()['data'].item()
    trial = currTrial['trial'].item()['data'].item()

    if catch == True:
        if extendedEOT == 0 and trial['instructTrial'] == 0:
            return True
        else:
            return False
    else:
        if extendedEOT == 0 and trial['instructTrial'] == 0 and trial['catchTrial'] != 1:
            return True
        else:
            return False


def fieldInTrial(fieldList, trial = None):
    '''
    Function will check whether all fields in a list are in the trial
    
    Inputs: trial (data struct from MATLAB)
            list of fields (list)
    Outputs: bool
    '''
    if trial == None:
        trial = currTrial

    for field in fieldList:
        if field not in trial.dtype.names:
            return False
   
    return True


def targetOnsetIndexStimDesc(stimDesc):
    '''
    fn will identify the index of the target onset stimulus
    fn will put out ListTypes subdictionary within stimDesc 
    
    Inputs: stimDesc (variable)
    Outputs: the index of the target onset stim (int)
    '''

    for count, d in enumerate(stimDesc['listTypes']):
                if 2 in d:
                    targetOnsetStim = count
                    break
    
    return count


def eyePosDurTrial(currTrial):
    '''
    fn will return a defaultdict with the converted x,y deg for each eye

    Inputs: trial (nd.array)
    Outputs: defaultdict
    '''
    eyesXYDeg = defaultdict(list)
    eyeLX = currTrial['eyeLXData'].item()['data'].item() 
    eyeLY = currTrial['eyeLYData'].item()['data'].item()
    eyeRX = currTrial['eyeRXData'].item()['data'].item()
    eyeRY = currTrial['eyeRYData'].item()['data'].item()
    eyeLeftCal = currTrial['eyeLeftCalibrationData'].item()['data'].item()['cal'].item()
    eyeRightCal = currTrial['eyeRightCalibrationData'].item()['data'].item()['cal'].item()
    count = min([len(eyeLX), len(eyeLY), len(eyeRX), len(eyeRY)])

    for s in range(0, count):
        xDegConvert = (eyeLX[s] * eyeLeftCal['m11'].item()) + (eyeLY[s] * eyeLeftCal['m21'].item()) + eyeLeftCal['tX'].item()
        eyesXYDeg['leftX'].append(xDegConvert)
        yDegConvert = (eyeLX[s] * eyeLeftCal['m12'].item()) + (eyeLY[s] * eyeLeftCal['m22'].item()) + eyeLeftCal['tY'].item()
        eyesXYDeg['leftY'].append(yDegConvert)
        xDegConvert = (eyeRX[s] * eyeRightCal['m11'].item()) + (eyeRY[s] * eyeRightCal['m21'].item()) + eyeRightCal['tX'].item()
        eyesXYDeg['rightX'].append(xDegConvert)
        yDegConvert = (eyeRX[s] * eyeRightCal['m12'].item()) + (eyeRY[s] * eyeRightCal['m22'].item()) + eyeRightCal['tY'].item()
        eyesXYDeg['rightY'].append(yDegConvert)

    return eyesXYDeg


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

    spikes = popMean + np.random.rand(1,numNeurons)           
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

    # print(spikes)


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



#testing time to open file
# t0 = time.time(); a = mat73.loadmat('./Documents/Grad School/Maunsell Lab/Analysis/MTCAN/Matlab Data/Meetz_2022_0114_MTNAN_m73.mat', use_attrdict = True); print(time.time() - t0)