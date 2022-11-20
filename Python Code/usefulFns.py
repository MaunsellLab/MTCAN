import scipy.io as sp
from scipy.optimize import curve_fit
from scipy import ndimage
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage.filters import gaussian_filter
from scipy import stats
from sklearn.metrics import r2_score
from itertools import combinations
import itertools
import numpy as np
import numpy.ma as ma
from numpy.linalg import inv
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.transforms as transforms
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import math
import os
from collections import defaultdict
import mat73
import seaborn as sns
import pandas as pd
import matplotlib.ticker as ticker
import time
from astropy.modeling import models, fitting
from pymatreader import read_mat
import glob


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


def loadMatFilePyMat(NHP, date, fileName):
    '''
    Loads the given matfile and assigns variables to access trial data

    Inputs: NHP (monkey Name), date of data collection (str), fileName (str)
    Outputs: variables, (dict)
    '''
    os.chdir(f'../{NHP}/{date}/')
    allTrials = read_mat(f'{fileName}')
    allTrialsData = allTrials['trials']
    header = allTrials['header']

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


def unitPrefNullDir(bSmooth):
    '''
    function will return the units preferred and null direction based
    off of the maximum response at either location when there is only one 
    stimulus.

    Inputs: 
        unitCount: unit's index in the units array
        meanSpikeReshaped: array of meanSpike counts for each stimulusIndex
    Outputs:
        prefDirection, nullDirection: the preferred and null direction for the neuron
    '''
    dirArray = np.array([0,60,120,180,240,300])

    maxLoc0 = max(bSmooth[6,:6])
    maxLoc1 = max(bSmooth[:6,6])
    if maxLoc0 > maxLoc1:
        prefDirection = dirArray[np.where(bSmooth[6,:]==maxLoc0)[0][0]]
    else:
        prefDirection = dirArray[np.where(bSmooth[:,6]==maxLoc1)[0][0]]
    nullDirection = (prefDirection + 180)%360

    return prefDirection, nullDirection 


def histSpikes(stimOnTimeS,stimOffTimeS,histPrePostMS,unitTimeStamps):
    '''
    Function will extract the time a unit spikes within the histogram window
    during a stimulus presentation. 

    Inputs:
        stimOnTimeS: when does the stimulus turn on (float)
        stimOffTimeS: when does the stimulus turn off (float)
        histPrePostMS: histogram pre and post window (int)
        unitTimeStamps: all the times the unit fired during that trial (array)
    
    Outputs:
        histStimSpikes: index values for when the unit fired within the hist window 
    '''

    stimOnPreSNEV = stimOnTimeS - histPrePostMS/1000
    stimOffPostSNEV = stimOffTimeS + histPrePostMS/1000
    histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV) &
                     (unitTimeStamps < stimOffPostSNEV))] - stimOnPreSNEV
    histStimSpikes = np.int32(histStimSpikes*1000)

    return histStimSpikes


def vonMises(x,x0,conc,I0):
    '''
    equation for a Von Mises fit
    '''
    return (np.exp(conc * np.cos(x - x0))) / (2*np.pi*I0*conc)


def vonMisesMatt(x,phase, kappa):
    # Von mises distribution
    z = np.exp(kappa * phase) / np.exp(kappa)
    return z / np.mean(z)


def vonMisesFit(x,y):
    '''
    apply curve_fit from scipy.optimize
    '''

    popt, pcov = curve_fit(vonMises, x,y)
    return popt


def logNormal(x,H,A,x0,sigma):
    '''
    equation for log-normal fot
    '''
    return H + A * np.exp(-(x-x0)**2 / (2*sigma**2))


def logNormalFit(x,y):
    '''
    apply curve_fit from scipy.optimize to fit a lognormal
    curve to speed tuning data
    '''
    x = np.log2(x)
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt,pcov = curve_fit(gauss,x,y, p0=[min(y), max(y), mean, sigma])
    return popt


def gauss(x, H, A, x0, sigma):
    '''
    equation for gaussian fit
    '''
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def gaussFit(x, y):
    '''
    apply curve_fit from scipy.optimize 
    '''
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma],
                   bounds=((0,0,0,0),(np.inf,1.50*np.max(y),np.inf,np.inf)))
    return popt


def gauss2dParams(neuronTuningMat):
    '''
    function will return the parameters to fit a 2D gaussian
    onto the neuron RF location heatmap
    returns mean vector (meanVec) and covariance matrix (covMat)
    '''

    com = ndimage.center_of_mass(neuronTuningMat)
    p_init = models.Gaussian2D(amplitude=1, x_mean=com[1], y_mean=com[0], x_stddev=None, 
                            y_stddev=None, theta=None, cov_matrix=None)
    yi, xi = np.indices(neuronTuningMat.shape)
    fit_p = fitting.LevMarLSQFitter()
    p = fit_p(p_init,xi,yi,neuronTuningMat)

    theta = p.theta[0] * 180/np.pi
    xStdDev = p.x_stddev[0]
    yStdDev = p.y_stddev[0]
    xMean = p.x_mean[0]
    yMean = p.y_mean[0]
    amp = p.amplitude[0]

    rho = np.cos(theta)
    covMat = np.array([[xStdDev**2,rho*xStdDev*yStdDev],
                    [rho*xStdDev*yStdDev,yStdDev**2]])
    meanVec = np.array([[xMean],[yMean]])

    return meanVec, covMat, p


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


def bhattCoef2D(m1,m2,cov1,cov2):
    '''
    This function will compute the Bhattacharyya Coefficient for two
    2D-Gaussian distriutions. This is used to measure how much overlap there
    is between two receptive fields.

    Inputs:
        m1 (np.array) - mean vector of the first 2D Gaussian
        m2 (np.array) - mean vector of the second 2D Gaussian
        cov1 (np.array) - covariance matrix of the first 2D Gaussian
        cov2 (np.array) - covariance matrix of the second 2D Gaussian
    Outputs:
        BC (float) - a value b/w 0 and 1 that defines how similar the curves are
    '''

    meanDiff = m1 - m2
    sigma = (cov1 + cov2)/2
    detSigma = np.linalg.det(sigma)
    detCov1 = np.linalg.det(cov1)
    detCov2 = np.linalg.det(cov2)

    X = 1/8 * np.dot(np.dot(np.transpose(meanDiff),inv(sigma)), meanDiff)
    Y = 1/2 * np.log(detSigma/np.sqrt(detCov1 * detCov2))
    BD = X + Y
    BC = np.exp(-BD)

    return BC
    

def smooth(y, box_pts):
    '''
    moving point average
    '''
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth


def spikeHistIndex(loc0Dir,loc0Con, loc1Dir, loc1Con):
    '''
    function will return the index for the each stimulus type PSTH
    '''
    histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == nullDir) & 
    (stimIndexDF['loc0 Contrast'] == lowC) & (stimIndexDF['loc1 Direction'] == nullDir)
    & (stimIndexDF['loc1 Contrast'] == lowC)][0]
    return histIndex


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
    
    return np.sort(units).astype(np.uint16)


def dirClosestToPref(unitPrefDir):
    '''
    function will take the unit's pref direction (based off gaussian fit)
    and return the directed tested that is closest to the preferred direction

    Inputs: unitPrefDir (unit's preferred direction)
    Outputs: prefDir, nullDir (unit's preferred and null direction)
    '''

    extDirArray = np.array([0,60,120,180,240,300,360])
    tempArr = abs(extDirArray-unitPrefDir)
    prefDirLoc = np.where(tempArr == min(tempArr))[0][0]

    if prefDirLoc == 6:
        prefDir = 0
    else:
        prefDir = extDirArray[prefDirLoc]
    nullDir = (prefDir + 180) % 360 

    return prefDir, nullDir


def normFunc0(fixed, BO,A,MU,SIG,S,al,c50,M):
    '''
    curve fit variables for my norm function, when loc0 has a stronger response
    
    BO: gaussian tuning curve baseline offset 
    A: gaussian tuning curve amplitude
    MU: guassian tuning curve mean
    SIG: gaussian tuning curve std dev
    S: scalar for other location gauss function
    al: alpha for normalization
    c50: normalization sigma
    M: baseline resp (blank stimulus)
    '''

    c0,l0,c1,l1 = fixed.T
    num = (c0 * (BO + A*np.exp(-(l0-MU)**2/(2*SIG ** 2)))) + (c1 * S * (BO + A*np.exp(-(l1-MU)**2/(2*SIG ** 2))))
    denom = c0 + (al * c1) + c50

    return  ((num/denom) + M).squeeze()


def normFunc1(fixed, BO, A, MU, SIG, S, al, c50, M):
    '''
    curve fit variables for my norm function, when loc1 has stronger response

    BO: gaussian tuning curve baseline offset 
    A: gaussian tuning curve amplitude
    MU: guassian tuning curve mean
    SIG: gaussian tuning curve std dev
    S: scalar for other location gauss function
    al: alpha for normalization
    c50: normalization sigma
    M: baseline resp (blank stimulus)
    '''

    c0,l0,c1,l1 = fixed.T
    num = (c0 * S * (BO + A*np.exp(-(l0-MU)**2/(2*SIG ** 2)))) + (c1 * (BO + A*np.exp(-(l1-MU)**2/(2*SIG ** 2))))
    denom = (al * c0) + c1 + c50

    return  ((num/denom) + M).squeeze()


def unitsInfo(units, corrTrials, allTrials):
    '''
    function will return the channel the unit was found and whether
    the unit is a multi-unit or single unit.

    Input: 
        units (list): list of active units
        corrTrials (list): list of correc trials
        allTrials (dict): all trials in the dataset
    Outputs: 
        unitChannel (list):  list of which channel the unit was found
    '''

    unitChannels = []
    for unit in units:
        breakNum = 0
        for corrTrial in corrTrials:
            currTrial = allTrials[corrTrial]
            spikeData = currTrial['spikeData']
            for count, sd in enumerate(spikeData['unit']):
                if sd == unit:
                    unitChannels.append(spikeData['channel'][count])
                    breakNum = 10
                    break
            if breakNum == 10:
                break

    return np.array(unitChannels)   


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


''' function fitting working 

def func(fixed, L0_0, L0_60, L0_120, L0_180, L0_240, L0_300, L1_0, L1_60, L1_120,
         L1_180, L1_240, L1_300, aL0, aL1, sig):
    c0,c1,l0,l1 = fixed
    L0 = np.array([L0_0, L0_60, L0_120, L0_180, L0_240, L0_300])
    L1 = np.array([L1_0, L1_60, L1_120, L1_180, L1_240, L1_300])
    return((c0*L0*l0).sum(-1) + (c1*L1*l1).sum(-1))/((aL0*c0[:,0])+(aL1*c1[:,0])+sig)

for unitCount, unit in enumerate(units):
    resp = np.reshape(spikeCountMat[unitCount][1:blocksDone+1,:],(169*blocksDone))
    fixParam = np.tile(np.arange(169), blocksDone)

    c0s, c1s, l0s, l1s = [], [], [], []
    direction_set = np.arange(0, 360, 60)
    for i in fixParam:
        c0 = stimIndexDict[i][0]['contrast']
        l0 = stimIndexDict[i][0]['direction']
        c1 = stimIndexDict[i][1]['contrast']
        l1 = stimIndexDict[i][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(6)
        l1_oh = np.zeros(6)
        l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
        l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1

        c0s.append(np.repeat(c0,6))
        c1s.append(np.repeat(c1,6))
        l0s.append(l0_oh)
        l1s.append(l1_oh)
        #append([c0,l0_oh,c1,l1_oh])
    c0s = np.array(c0s)
    c1s = np.array(c1s)
    l0s = np.array(l0s)
    l1s = np.array(l1s)
    #fix = np.array(fix)
    #print(fix.shape)
    pOpt, pCov = curve_fit(func, (c0s, c1s, l0s, l1s), resp, bounds=(
        (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
        (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        np.inf, np.inf, np.inf, 1,1,1)))
    print(unit,pOpt)

'''

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