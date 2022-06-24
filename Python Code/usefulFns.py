import scipy.io as sp
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy import stats
from itertools import combinations
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import math
import os
from collections import defaultdict
import mat73
import seaborn as sns
import matplotlib.ticker as ticker
import time

def loadMatFile73(NHP, date, fileName):
    '''
    Loads the given matfile and assigns variables to access trial data

    Inputs: matfile name, (str)
    Outputs: variables, (nd.array)
    '''
    os.chdir(f'../{NHP}/{date}/')
    allTrials = mat73.loadmat(f'{fileName}', use_attrdict = True)
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


def correctTrialsMTX(allTrials):
    '''
    Function will filter through allTrials and return a list of 
    correct non-instruct trials in MTN/MTC. This function checks
    for valid trialCertify trials (!=0)

    Inputs: allTrials (list of trials (nd.array))
    Outputs: correctTrials(list): index of valid trials from allTrials
    '''
    correctTrials = []
    for trialCount, currTrial in enumerate(allTrials):
        trial = currTrial['trial']['data']
        extendedEOT = currTrial['extendedEOT']['data']
        trialCertify = currTrial['trialCertify']['data']
        if trial['instructTrial'] != 1 and extendedEOT == 0 and trialCertify == 0:
            correctTrials.append(trialCount)
    
    return correctTrials


def gauss(x, H, A, x0, sigma):
    '''
    equation for gaussian fit
    '''
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def gauss_fit(x, y):
    '''
    apply curve_fit from scipy.optimize 
    '''
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt


def bhattCoef(m1, m2, v1, v2):
    '''
    This function will compute the Bhattacharyya Coefficient of two 
    tuning curves. The BC measures how similar two normal distributions are. 

    Inputs:
        m1 (float) - mean of the first tuning curve
        m2 (float) - mean of the second tuning curve
        v1 (float) - variance of the first tuning curve
        v2 (float) - variance of the second tuning curve
    Outputs: 
        BC (float) - a value b/w 0 and 1 that defines how similar the curves are.
    '''

    BD = (1/4*np.log(1/4*((v1/v2) + (v2/v1) + 2))) + \
         (1/4*(((m1-m2)**2)/(v1+v2)))
    BC = np.exp(-BD)

    return BC


def smooth(y, box_pts):
    '''
    moving point average
    '''
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth


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
    eyeLX = currTrial['eyeLXData']['data']
    eyeLY = currTrial['eyeLYData']['data']
    eyeRX = currTrial['eyeRXData']['data']
    eyeRY = currTrial['eyeRYData']['data']
    eyeLeftCal = currTrial['eyeLeftCalibrationData']['data']['cal']
    eyeRightCal = currTrial['eyeRightCalibrationData']['data']['cal']
    count = min([len(eyeLX), len(eyeLY), len(eyeRX), len(eyeRY)])

    for s in range(0, count):
        xDegConvert = (eyeLX[s] * eyeLeftCal['m11']) + (eyeLY[s] * eyeLeftCal['m21']) + eyeLeftCal['tX']
        eyesXYDeg['leftX'].append(xDegConvert)
        yDegConvert = (eyeLX[s] * eyeLeftCal['m12']) + (eyeLY[s] * eyeLeftCal['m22']) + eyeLeftCal['tY']
        eyesXYDeg['leftY'].append(yDegConvert)
        xDegConvert = (eyeRX[s] * eyeRightCal['m11']) + (eyeRY[s] * eyeRightCal['m21']) + eyeRightCal['tX']
        eyesXYDeg['rightX'].append(xDegConvert)
        yDegConvert = (eyeRX[s] * eyeRightCal['m12']) + (eyeRY[s] * eyeRightCal['m22']) + eyeRightCal['tY']
        eyesXYDeg['rightY'].append(yDegConvert)

    return eyesXYDeg


def activeUnits(unitData, allTrials):

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
    
    return np.sort(units)


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
    tuningMat[0] = np.arange(360,720,60)
    tcDictionary = {}

    for i in range(1, tuningMat.shape[0]):
        np.random.seed(i)
        amp = np.random.randint(15,30)
        y_translate = np.random.randint(30,50)
        x_translate = np.random.randint(60,120)
        tuningMat[i,:] = (amp * np.sin((tuningMat[0,:] * (np.pi / 180)) + x_translate)) + y_translate

    for i, neuron in enumerate(tuningMat[1:,:]):
        tcDictionary[i+1] = {}
        for j, dirResp in enumerate(neuron):
            tcDictionary[i+1][tuningMat[0][j]] = dirResp

    return tuningMat, tcDictionary


# stimCount test
# correctTrials = correctTrialsGRF(allTrials)
# for corrTrial in correctTrials:
#     currTrial = allTrials[corrTrial]
#         if 'numMap0Stim' in currTrial:
    #         map0StimLim = int(currTrial['numMap0Stim']['data'].tolist())
    #         map0Count = 0
    #         stimDesc = currTrial['stimDesc']
    #         for sCount, stim in enumerate(stimDesc['data']):
    #             if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
    #                 aziIndex = int(stim['azimuthIndex'])
    #                 eleIndex = int(stim['elevationIndex'])
    #                 stimCount[eleIndex][aziIndex] = stimCount[eleIndex][aziIndex] + 1
    #                 map0Count += 1


#testing time to open file
# t0 = time.time(); a = mat73.loadmat('./Documents/Grad School/Maunsell Lab/Analysis/MTCAN/Matlab Data/Meetz_2022_0114_MTNAN_m73.mat', use_attrdict = True); print(time.time() - t0)