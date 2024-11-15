import scipy.io as sp
from scipy.optimize import curve_fit
import scipy.optimize
from scipy import ndimage
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage.filters import gaussian_filter
from scipy import stats
from scipy.stats import mannwhitneyu
from sklearn.metrics import r2_score
from sklearn.metrics import explained_variance_score
from itertools import combinations
import itertools
import numpy as np
import random
import numpy.ma as ma
from numpy.linalg import inv
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.transforms as transforms
import matplotlib.lines as lines
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
from binsreg import *
from scipy.stats import f_oneway
import glob
import matplotlib as mpl
import matplotlib.lines as mlines
from scipy.stats import wilcoxon


# fig saving params
mpl.rcParams['pdf.fonttype'] = 42


def loadMatFile73(NHP, date, fileName):
    """
    Loads the given matfile and assigns variables to access trial data

    Inputs: matfile name, (str)
    Outputs: variables, (nd.array)
    """

    os.chdir(f'../{NHP}/{date}/')
    allTrials = mat73.loadmat(f'{fileName}', use_attrdict=True)
    allTrialsData = allTrials.trials
    header = allTrials.header

    return allTrialsData, header


def loadMatFilePyMat(NHP, date, fileName):
    """
    Loads the given matfile and assigns variables to access trial data

    Inputs: NHP (monkey Name), date of data collection (str), fileName (str)
    Outputs: variables, (dict)
    """

    os.chdir(f'../{NHP}/{date}/')
    allTrials = read_mat(f'{fileName}')
    allTrialsData = allTrials['trials']
    header = allTrials['header']

    return allTrialsData, header


def loadMatFile(fileName):
    """
    Loads the given matfile and assigns variables to access trial data

    Inputs: matfile name, (str)
    Outputs: variables, (nd.array)
    """
    
    allTrials = sp.loadmat(f'../Matlab Data/{fileName}', squeeze_me=True)
    allTrialsData = allTrials['trials']
    header = allTrials['header']

    return allTrialsData, header


def correctTrialsGRF(allTrials):
    """
    fn that will filter for correct non-instruct trials (valid) in GRF

    inputs: allTrials (list of trials)
    outputs: correctTrials (list): index of valid trials from allTrials
    """

    correctTrials = []
    for trialCount, currTrial in enumerate(allTrials):
        trial = currTrial['trial']['data']
        trialEnd = currTrial['trialEnd']['data']
        if trial['instructTrial'] != 1 and trialEnd == 0:
            correctTrials.append(trialCount)
    
    return correctTrials


def correctTrialsMTX(allTrials):
    """
    Function will filter through allTrials and return a list of 
    correct non-instruct trials in MTN/MTC. This function checks
    for valid trialCertify trials (!=0)

    Inputs: allTrials (list of trials (nd.array))
    Outputs: correctTrials(list): index of valid trials from allTrials
    """

    correctTrials = []
    for trialCount, currTrial in enumerate(allTrials):
        trial = currTrial['trial']['data']
        extendedEOT = currTrial['extendedEOT']['data']
        trialCertify = currTrial['trialCertify']['data']
        if trial['instructTrial'] != 1 and extendedEOT == 0 and trialCertify == 0:
            correctTrials.append(trialCount)
    
    return correctTrials


def unitPrefNullDir(bSmooth):
    """
    function will return the units preferred and null direction based
    off of the maximum response at either location when there is only one 
    stimulus.

    Inputs: 
        unitCount: unit's index in the units array
        meanSpikeReshaped: array of meanSpike counts for each stimulusIndex
    Outputs:
        prefDirection, nullDirection: the preferred and null direction for the neuron
    """

    dirArray = np.array([0, 60, 120, 180, 240, 300])

    maxLoc0 = max(bSmooth[6, :6])
    maxLoc1 = max(bSmooth[:6, 6])
    if maxLoc0 > maxLoc1:
        prefDirection = dirArray[np.where(bSmooth[6, :] == maxLoc0)[0][0]]
    else:
        prefDirection = dirArray[np.where(bSmooth[:, 6] == maxLoc1)[0][0]]
    nullDirection = (prefDirection + 180) % 360

    return prefDirection, nullDirection 


def histSpikes(stimOnTimeS, stimOffTimeS, histPrePostMS, unitTimeStamps):
    """
    Function will extract the time a unit spikes within the histogram window
    during a stimulus presentation. 

    Inputs:
        stimOnTimeS: when does the stimulus turn on (float)
        stimOffTimeS: when does the stimulus turn off (float)
        histPrePostMS: histogram pre and post window (int)
        unitTimeStamps: all the times the unit fired during that trial (array)
    
    Outputs:
        histStimSpikes: index values for when the unit fired within the hist window 
    """

    stimOnPreSNEV = stimOnTimeS - histPrePostMS/1000
    stimOffPostSNEV = stimOffTimeS + histPrePostMS/1000
    histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV) &
                                    (unitTimeStamps < stimOffPostSNEV))] - stimOnPreSNEV
    histStimSpikes = np.int32(histStimSpikes*1000)

    return histStimSpikes


def contrastFn(c, r0, rMax, c50, n):
    """
    equation for fitting contrast response functions (Naka-Rushton)
    rMax is maximum response, c50 is the semi-saturation constant
    and n is the exponent that describes the slope of the function
    """

    return r0 + ((rMax * (c ** n)) / ((c ** n) + (c50 ** n)))


def contrastFnNoBaseline(c, rMax, c50, n):
    """
    equation for fitting contrast response functions (Naka-Rushton)
    rMax is maximum response, c50 is the semi-saturation constant
    and n is the exponent that describes the slope of the function
    without a baseline, since responses start at 0 (baseline substracted
    and normalized to 1)
    """

    return (rMax * (c ** n)) / ((c ** n) + (c50 ** n))


def confidenceIntervalCRF(popt, pcov, x, confidence=0.95):
    """
    function to calculate confidence interval for contrast response function fits
    """

    alpha = 1.0 - confidence
    n = len(x)
    p = len(popt)
    dof = max(0, n - p)
    tval = np.sqrt(np.diag(pcov)) * 1.96
    lower = contrastFn(x, *(popt - tval))
    upper = contrastFn(x, *(popt + tval))

    return lower, upper


def gauss(x, H, A, x0, sigma):
    """
    equation for gaussian fit
    """

    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def gaussFit(x, y):
    """
    apply curve_fit from scipy.optimize 
    """
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma],
                           bounds=((-np.inf, 0, -np.inf, 0),
                                   (np.inf, np.inf, np.inf, np.inf)),
                           maxfev=10000000)
    return popt


def gauss2dParams(neuronTuningMat):
    """
    function will return the parameters to fit a 2D gaussian
    onto the neuron RF location heatmap
    returns mean vector (meanVec) and covariance matrix (covMat)
    http://prob140.org/sp18/textbook/notebooks-md/24_01_Bivariate_Normal_Distribution.html
    """

    com = ndimage.center_of_mass(neuronTuningMat)
    p_init = models.Gaussian2D(amplitude=1, x_mean=com[1], y_mean=com[0], x_stddev=None, 
                               y_stddev=None, theta=None, cov_matrix=None)
    yi, xi = np.indices(neuronTuningMat.shape)
    fit_p = fitting.LevMarLSQFitter()
    p = fit_p(p_init, xi, yi, neuronTuningMat)

    theta = p.theta[0] * 180/np.pi
    xStdDev = p.x_stddev[0]
    yStdDev = p.y_stddev[0]
    xMean = p.x_mean[0]
    yMean = p.y_mean[0]
    # amp = p.amplitude[0]

    rho = np.cos(theta)
    covMat = np.array([[xStdDev**2, rho*xStdDev*yStdDev],
                      [rho*xStdDev*yStdDev, yStdDev**2]])
    meanVec = np.array([[xMean], [yMean]])

    return meanVec, covMat, p


def bhattCoef(m1, m2, v1, v2):
    """
    This function will compute the Bhattacharyya Coefficient of two 
    tuning curves. The BC measures how similar two normal distributions are. 

    Inputs:
        m1 (float) - mean of the first tuning curve
        m2 (float) - mean of the second tuning curve
        v1 (float) - variance of the first tuning curve
        v2 (float) - variance of the second tuning curve
    Outputs: 
        BC (float) - a value b/w 0 and 1 that defines how similar the curves are.
    """

    BD = (1/4*np.log(1/4*((v1/v2) + (v2/v1) + 2))) + \
         (1/4*(((m1-m2)**2)/(v1+v2)))
    BC = np.exp(-BD)

    return BC


def bhattCoef2D(m1, m2, cov1, cov2):
    """
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
    """

    meanDiff = m1 - m2
    sigma = (cov1 + cov2)/2
    detSigma = np.linalg.det(sigma)
    detCov1 = np.linalg.det(cov1)
    detCov2 = np.linalg.det(cov2)

    X = 1/8 * np.dot(np.dot(np.transpose(meanDiff), inv(sigma)), meanDiff)
    Y = 1/2 * np.log(detSigma/np.sqrt(detCov1 * detCov2))
    BD = X + Y
    BC = np.exp(-BD)

    return BC
    

def smooth(y, box_pts):
    """
    moving point average
    """
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth


def upperTriMasking(matrix):
    """
    this function will return the indices of the top right triangle
    of the square matrix excluding the diagonal

    Inputs: square matrix
    returns: indices of the top right triangle of the matrix
    """
    m = matrix.shape[0]
    r = np.arange(m)
    mask = r[:, None] < r
    return matrix[mask]


def rejectOutliers(data, m=2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s < m]


def spikeHistIndex(loc0Dir, loc0Con, loc1Dir, loc1Con):
    """
    function will return the index for the each stimulus type PSTH
    """

    histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == nullDir) & 
    (stimIndexDF['loc0 Contrast'] == lowC) & (stimIndexDF['loc1 Direction'] == nullDir)
    & (stimIndexDF['loc1 Contrast'] == lowC)][0]
    return histIndex


def fieldInTrial(fieldList, trial=None):
    """
    Function will check whether all fields in a list are in the trial
    
    Inputs: trial (data struct from MATLAB)
            list of fields (list)
    Outputs: bool
    """
    if trial == None:
        trial = currTrial

    for field in fieldList:
        if field not in trial.dtype.names:
            return False
   
    return True


def eyePosDurTrial(currTrial):
    """
    fn will return a defaultdict with the converted x,y deg for each eye

    Inputs: trial (nd.array)
    Outputs: defaultdict
    """

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
    """
    function returns the active units across trials for a session as a list

    Inputs: unitData (str): Are we using fakeData or spikeData
    Outputs: units (list): active units for a sessioon
    """
    units = []
    for currTrial in allTrials:
        if unitData in currTrial:
            uniqueUnits = np.unique(currTrial[unitData]['unit'])
            for unique in uniqueUnits:
                if unique not in units:
                    units.append(unique)
    
    return np.sort(units).astype(np.uint16)


def dirClosestToPref(unitPrefDir):
    """
    function will take the unit's pref direction (based off gaussian fit)
    and return the directed tested that is closest to the preferred direction

    Inputs: unitPrefDir (unit's preferred direction)
    Outputs: prefDir, nullDir (unit's preferred and null direction from
            directions tested)
    """

    extDirArray = np.array([0, 60, 120, 180, 240, 300, 360])
    tempArr = abs(extDirArray-unitPrefDir)
    prefDirLoc = np.where(tempArr == min(tempArr))[0][0]

    if prefDirLoc == 6:
        prefDir = 0
    else:
        prefDir = extDirArray[prefDirLoc]
    nullDir = (prefDir + 180) % 360 

    return prefDir, nullDir


def fixedValsForGenericNorm(stimMatReIndex, stimIndexDict):
    """
    function will create the fixed values for running the
    generic normalization curve_fit
    """

    c0s, c1s, l0s, l1s = [], [], [], []
    direction_set = np.arange(0, 360, 60)
    if type(stimMatReIndex) == int:
        c0 = stimIndexDict[stimMatReIndex][0]['contrast']
        l0 = stimIndexDict[stimMatReIndex][0]['direction']
        c1 = stimIndexDict[stimMatReIndex][1]['contrast']
        l1 = stimIndexDict[stimMatReIndex][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(6)
        l1_oh = np.zeros(6)
        l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
        l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
        c0s.append(np.repeat(c0, 6))
        c1s.append(np.repeat(c1, 6))
        l0s.append(l0_oh)
        l1s.append(l1_oh)
    else:
        for i in stimMatReIndex:
            c0 = stimIndexDict[i][0]['contrast']
            l0 = stimIndexDict[i][0]['direction']
            c1 = stimIndexDict[i][1]['contrast']
            l1 = stimIndexDict[i][1]['direction']

            # Make one-hot encoding of l0 and l1
            l0_oh = np.zeros(6)
            l1_oh = np.zeros(6)
            l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
            l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
            c0s.append(np.repeat(c0, 6))
            c1s.append(np.repeat(c1, 6))
            l0s.append(l0_oh)
            l1s.append(l1_oh)

    return np.array(c0s), np.array(c1s), np.array(l0s), np.array(l1s)


def fixedValsForEMSGen(stimMatReIndex, stimIndexDict):
    """
    function will create the fixed values for running the
    generic EMS normalization curve_fit
    """

    c0s, c1s, l0s, l1s = [], [], [], []
    direction_set = np.arange(0, 360, 60)
    if type(stimMatReIndex) == int:
        c0 = stimIndexDict[stimMatReIndex][0]['contrast']
        l0 = stimIndexDict[stimMatReIndex][0]['direction']
        c1 = stimIndexDict[stimMatReIndex][1]['contrast']
        l1 = stimIndexDict[stimMatReIndex][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(6)
        l1_oh = np.zeros(6)
        l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
        l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
        c0s.append(np.repeat(c0,6))
        c1s.append(np.repeat(c1,6))
        l0s.append(l0_oh)
        l1s.append(l1_oh)
    else:
        for i in stimMatReIndex:
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

    return np.array(c0s), np.array(c1s), np.array(l0s), np.array(l1s)


def fixedValsForPairedStimL0L6(stimMatReIndex, stimIndexDict, L0L6Resp):
    """
    function will take in the sequence of stimulus indices to be fit
    and return the fixed values for that posistion, including the
    L0L6 parameters for resp at each direction:
    loc 0: contrast
    loc 0: direction
    loc 1: contrast
    loc 1: direction
    """

    c0s, c1s, l0s, l1s, L, S = [], [], [], [], [], []
    direction_set = np.arange(0, 360, 60)
    if type(stimMatReIndex) == int:
        c0 = stimIndexDict[stimMatReIndex][0]['contrast']
        l0 = stimIndexDict[stimMatReIndex][0]['direction']
        c1 = stimIndexDict[stimMatReIndex][1]['contrast']
        l1 = stimIndexDict[stimMatReIndex][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(6)
        l1_oh = np.zeros(6)
        l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
        l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
        c0s.append(np.repeat(c0, 6))
        c1s.append(np.repeat(c1, 6))
        l0s.append(l0_oh)
        l1s.append(l1_oh)
        L.append(L0L6Resp[:-1])
        S.append(np.repeat(L0L6Resp[-1], 6))
    else:
        for i in stimMatReIndex:
            c0 = stimIndexDict[i][0]['contrast']
            l0 = stimIndexDict[i][0]['direction']
            c1 = stimIndexDict[i][1]['contrast']
            l1 = stimIndexDict[i][1]['direction']

            # Make one-hot encoding of l0 and l1
            l0_oh = np.zeros(6)
            l1_oh = np.zeros(6)
            l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
            l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
            c0s.append(np.repeat(c0, 6))
            c1s.append(np.repeat(c1, 6))
            l0s.append(l0_oh)
            l1s.append(l1_oh)
            L.append(L0L6Resp[:-1])
            S.append(np.repeat(L0L6Resp[-1], 6))

    return (np.array(c0s), np.array(c1s), np.array(l0s), np.array(l1s),
            np.array(L), np.array(S))


def fixedValsForEMSGenCondensed(stimMatReIndex, stimIndexDict, nullDir, prefDir):
    """
    function will create the fixed values for running the
    generic EMS normalization curve_fit
    """

    c0s, c1s, l0s, l1s = [], [], [], []
    directionSet = np.array([nullDir, prefDir])
    if type(stimMatReIndex) == int:
        c0 = stimIndexDict[stimMatReIndex][0]['contrast']
        l0 = stimIndexDict[stimMatReIndex][0]['direction']
        c1 = stimIndexDict[stimMatReIndex][1]['contrast']
        l1 = stimIndexDict[stimMatReIndex][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(2)
        l1_oh = np.zeros(2)
        l0_oh[np.argwhere(directionSet == l0).squeeze()] = 1
        l1_oh[np.argwhere(directionSet == l1).squeeze()] = 1
        c0s.append(np.repeat(c0, 2))
        c1s.append(np.repeat(c1, 2))
        l0s.append(l0_oh)
        l1s.append(l1_oh)
    else:
        for i in stimMatReIndex:
            c0 = stimIndexDict[i][0]['contrast']
            l0 = stimIndexDict[i][0]['direction']
            c1 = stimIndexDict[i][1]['contrast']
            l1 = stimIndexDict[i][1]['direction']

            # Make one-hot encoding of l0 and l1
            l0_oh = np.zeros(2)
            l1_oh = np.zeros(2)
            l0_oh[np.argwhere(directionSet == l0).squeeze()] = 1
            l1_oh[np.argwhere(directionSet == l1).squeeze()] = 1
            c0s.append(np.repeat(c0, 2))
            c1s.append(np.repeat(c1, 2))
            l0s.append(l0_oh)
            l1s.append(l1_oh)

    return np.array(c0s), np.array(c1s), np.array(l0s), np.array(l1s)


def fixedValsCurveFitForPairedStim(prefDir, stimMatReIndex, stimIndexDict, gParams):
    """
    function will take in the sequence of stimulus indices to be fit
    and return the fixed values for that posistion, including the
    gaussian parameters for tuning curve:
    loc 0: contrast
    loc 0: direction
    loc 1: contrast
    loc 1: direction
    """

    BO = gParams[0]
    A = gParams[1]
    M = gParams[2]
    Sigma = gParams[3]
    Scal = gParams[4]
    fixedVals = []
    if type(stimMatReIndex) == int:
        c0 = stimIndexDict[stimMatReIndex][0]['contrast']
        l0 = stimIndexDict[stimMatReIndex][0]['direction']
        if abs(l0 - prefDir) > 180:
            if prefDir < 180:
                l0 = prefDir - (360 - (l0 - prefDir))
            else:
                l0 = prefDir + (360 - abs(l0 - prefDir))
        c1 = stimIndexDict[stimMatReIndex][1]['contrast']
        l1 = stimIndexDict[stimMatReIndex][1]['direction']
        if abs(l1 - prefDir) > 180:
            if prefDir < 180:
                l1 = prefDir - (360 - (l1 - prefDir))
            else:
                l1 = prefDir + (360 - abs(l1 - prefDir))
        fixedVals.append((c0, l0, c1, l1, BO, A, M, Sigma, Scal))
    else:
        for i in stimMatReIndex:
            c0 = stimIndexDict[i][0]['contrast']
            l0 = stimIndexDict[i][0]['direction']
            if abs(l0 - prefDir) > 180:
                if prefDir < 180:
                    l0 = prefDir - (360 - (l0 - prefDir))
                else:
                    l0 = prefDir + (360 - abs(l0 - prefDir))
            c1 = stimIndexDict[i][1]['contrast']
            l1 = stimIndexDict[i][1]['direction']
            if abs(l1 - prefDir) > 180:
                if prefDir < 180:
                    l1 = prefDir - (360 - (l1 - prefDir))
                else:
                    l1 = prefDir + (360 - abs(l1 - prefDir))
            fixedVals.append((c0, l0, c1, l1, BO, A, M, Sigma, Scal))
    fixedVals = np.array(fixedVals)

    return fixedVals


def fixedValsForCurveFit(prefDir, stimMatReIndex, stimIndexDict):
    """
    function will take in the sequence of stimulus indexes to be fit
    and return the fixed values for that position:
    loc 0: contrast
    loc 0: direction
    loc 1: contrast
    loc 1: direction
    """

    fixedVals = []
    if type(stimMatReIndex) == int:
        c0 = stimIndexDict[stimMatReIndex][0]['contrast']
        l0 = stimIndexDict[stimMatReIndex][0]['direction']
        if abs(l0 - prefDir) > 180:
            if prefDir < 180:
                l0 = prefDir - (360 - (l0 - prefDir))
            else:
                l0 = prefDir + (360 - abs(l0 - prefDir))
        c1 = stimIndexDict[stimMatReIndex][1]['contrast']
        l1 = stimIndexDict[stimMatReIndex][1]['direction']
        if abs(l1 - prefDir) > 180:
            if prefDir < 180:
                l1 = prefDir - (360 - (l1 - prefDir))
            else:
                l1 = prefDir + (360 - abs(l1 - prefDir))
        fixedVals.append((c0, l0, c1, l1))
    else:
        for i in stimMatReIndex:
            c0 = stimIndexDict[i][0]['contrast']
            l0 = stimIndexDict[i][0]['direction']
            if abs(l0 - prefDir) > 180:
                if prefDir < 180:
                    l0 = prefDir - (360 - (l0 - prefDir))
                else:
                    l0 = prefDir + (360 - abs(l0 - prefDir))
            c1 = stimIndexDict[i][1]['contrast']
            l1 = stimIndexDict[i][1]['direction']
            if abs(l1 - prefDir) > 180:
                if prefDir < 180:
                    l1 = prefDir - (360 - (l1 - prefDir))
                else:
                    l1 = prefDir + (360 - abs(l1 - prefDir))
            fixedVals.append((c0, l0, c1, l1))

    fixedVals = np.array(fixedVals)

    return fixedVals


def lightenColor(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """

    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def unitsInfo(units, corrTrials, allTrials):
    """
    function will return the channel the unit was found

    Input: 
        units (list): list of active units
        corrTrials (list): list of correct trials
        allTrials (dict): all trials in the dataset
    Outputs: 
        unitChannel (list):  list of which channel the unit was found
    """

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
    """
    this function will generate a poisson spike train and return the normalized
    response for the RF for a given stimulus configuration
    Inputs:
        x: neuron number
        index: the index of the corresponding stimulus configuration
        sigma: semisaturation constant from neuron's contrast response function
    Outputs:
        RFSpikes: normalized response 

    spikes_4rows = np.tile(spikes, (4,1))
    """

    numNeurons = len(units)

    C0 = stimIndexDict[index][0]['contrast']
    C1 = stimIndexDict[index][1]['contrast']
    L0 = tcDict[1][stimIndexDict[index][0]['direction']]
    L1 = tcDict[1][stimIndexDict[index][1]['direction']]
    sigma = 0.1
    expectedNormSpikeRate = int(((C0*L0) + (C1*L1))/(C0 + C1 + sigma))
    stimDur = 493
    popMean = expectedNormSpikeRate/(1000/stimDur)

    spikes = popMean + np.random.rand(1, numNeurons)
    R = np.zeros((numNeurons, numNeurons))
    for neuronI in range(numNeurons):
        for neuronJ in range(numNeurons):
            if neuronI == neuronJ:
                R[neuronI, neuronJ] = 1
            else:
                R[neuronI, neuronJ] = 0.1

    L = np.linalg.cholesky(R)
    spikes = np.matmul(spikes, L)
    spikes = np.around(spikes) 

    for count, i in enumerate(spikes[0]):
        if i != 0:
            unit = units[count]
            channelIdentity = int(unit[0:unit.find('_')])
            channel = np.array([channelIdentity] * stimDur)
            spikeTimeMS = (np.sort(np.random.choice(np.arange(stimDur), int(i),
                           replace=False)))/1000
            currTrial['spikeData']['timeStamp'] = np.append(currTrial['spikeData'] \
            ['timeStamp'], stimOnTimeSNEV + spikeTimeMS, 0)
            currTrial['spikeData']['unit'] = np.append(currTrial['spikeData'] \
            ['unit'], [unit] * len(spikeTimeMS), 0)
            currTrial['spikeData']['channel'] = np.append(currTrial['spikeData'] \
            ['channel'], [channelIdentity] * len(spikeTimeMS), 0)


def poissonArrivals(stimOnTimeS, lam, duration, doDecay=False, doRamp=False):

    # Decay lambda over the duration of the spike generation period
    a = -lam / (2 * duration**2)
    if doRamp:
        a *= -1
    t = 0
    spikeTimes = []
    while t < duration:
        dec_lam = (a * t**2 + lam) if doDecay else lam
        t += int(random.expovariate(dec_lam) * 1000)
        if t < duration:
            spikeTimes.append(t)

    return stimOnTimeS + (np.array(spikeTimes)/1000)


def randTuningCurve(numNeurons):
    """
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
    """

    tuningMat = np.zeros((numNeurons + 1, 6))
    tuningMat[0] = np.arange(360, 720, 60)
    tcDictionary = {}

    for i in range(1, tuningMat.shape[0]):
        np.random.seed(i)
        amp = np.random.randint(15, 30)
        y_translate = np.random.randint(30, 50)
        x_translate = np.random.randint(60, 120)
        tuningMat[i, :] = (amp * np.sin((tuningMat[0, :] * (np.pi / 180)) + x_translate)) + y_translate

    for i, neuron in enumerate(tuningMat[1:, :]):
        tcDictionary[i+1] = {}
        for j, dirResp in enumerate(neuron):
            tcDictionary[i+1][tuningMat[0][j]] = dirResp

    return tuningMat, tcDictionary


def vonMises(x, x0, conc, I0):
    """
    equation for a Von Mises fit
    """

    return (np.exp(conc * np.cos(x - x0))) / (2*np.pi*I0*conc)


def vonMisesMatt(x, phase, kappa):
    # Von mises distribution
    z = np.exp(kappa * phase) / np.exp(kappa)
    return z / np.mean(z)


def vonMisesFit(x, y):
    """
    apply curve_fit from scipy.optimize
    """

    popt, pcov = curve_fit(vonMises, x, y)
    return popt


def logNormal(x, H, A, x0, sigma):
    """
    equation for log-normal fot
    """
    return H + A * np.exp(-(x-x0)**2 / (2*sigma**2))


def logNormalFit(x, y):
    """
    apply curve_fit from scipy.optimize to fit a lognormal
    curve to speed tuning data
    """

    x = np.log2(x)
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt


def driverFunc(x, fixedVals, resp):
    """
    this function will aide in estimating my normalization
    equation parameters by minimizing the sum of square error
    """

    yNew = genericNormNoScalar(fixedVals, *x)
    yErr = np.sum((yNew - resp) ** 2)

    return yErr


def driverFuncCondensend(x, fixedVals, resp):
    """
    this function will aide in estimating my normalization
    equation parameters by minimizing the sum of square error
    """

    yNew = genNormCondensed(fixedVals, *x)
    yErr = np.sum((yNew - resp) ** 2)

    return yErr


def driverFuncCondensend1(x, fixedVals, resp):
    """
    this function will aide in estimating my normalization
    equation parameters by minimizing the sum of square error
    """

    yNew = genNormCondensed1(fixedVals, *x)
    yErr = np.sum((yNew - resp) ** 2)

    return yErr


def genericNormNoScalar(fixed, L0_0, L0_60, L0_120, L0_180, L0_240, L0_300,
                        L1_0, L1_60, L1_120, L1_180, L1_240, L1_300,
                        al1):
    """
    this function applies the generic normalization equation without a scalar
    at the weaker location. Instead, it gives that location its own L0-L6 value
    """

    c0, c1, l0, l1 = fixed
    L0 = np.array([L0_0, L0_60, L0_120, L0_180, L0_240, L0_300])
    L1 = np.array([L1_0, L1_60, L1_120, L1_180, L1_240, L1_300])

    # generic norm
    num = ((c0 * (L0 ** 2) * l0).sum(-1) + (c1 * (L1 ** 2) * l1).sum(-1))
    denom = (((1 ** 2) * c0[:, 0]) + ((al1 ** 2) * c1[:, 0]))

    return num/denom

    # # ems
    # loc0 = ((c0 * (L0 ** 2) * l0).sum(-1)) / (c0[:, 0] + ((al1 ** 2) * c1[:, 0]))
    # loc1 = ((c1 * (L1 ** 2) * l1).sum(-1)) / (((1 ** 2) * c0[:, 0]) + c1[:, 0])
    #
    # return loc0 + loc1


def genNormCondensed(fixed, L0_0, L0_1, L1_0, L1_1, al1):
    """
    generic normalization condensed
    """

    c0, c1, l0, l1 = fixed
    L0 = np.array([L0_0, L0_1])
    L1 = np.array([L1_0, L1_1])

    # generic norm
    num = ((c0 * (L0 ** 2) * l0).sum(-1) + (c1 * (L1 ** 2) * l1).sum(-1))
    denom = ((1 * c0[:, 0]) + ((al1 ** 2) * c1[:, 0]))

    return num / denom

    # # ems
    # loc0 = (c0 * (l0 * (L0 ** 2))).sum(-1) / (
    #         c0[:, 0] + (c1[:, 0] * (al1 ** 2)))
    #
    # loc1 = (c1 * (l1 * (L1 ** 2))).sum(-1) / (
    #        (c0[:, 0] * 1) + c1[:, 0])
    #
    # return loc0 + loc1

    # # generic norm
    # num = ((c0 * L0 * l0).sum(-1) + (c1 * L1 * l1).sum(-1))
    # denom = (1 + (al1 * c0[:, 0] * c1[:, 0]) + sig)
    #
    # return num/denom

    # def genNormCondensed(fixed, L0_0, L0_1, L1_0, L1_1, al1, sig):
    #     """
    #     generic normalization condensed
    #     """
    #
    #     c0, c1, l0, l1 = fixed
    #     L0 = np.array([L0_0, L0_1])
    #     L1 = np.array([L1_0, L1_1])
    #
    #     # generic norm
    #     num = ((c0 * L0 * l0).sum(-1) + (c1 * L1 * l1).sum(-1))
    #     denom = ((1 * c0[:, 0]) + (al1 * c1[:, 0]) + sig)
    #
    #     return num / denom


def genNormCondensed1(fixed, L0_0, L0_1, L1_0, L1_1, al1, sig):
    """
    generic normalization condensed
    """

    c0, c1, l0, l1 = fixed
    L0 = np.array([L0_0, L0_1])
    L1 = np.array([L1_0, L1_1])

    # generic norm
    num = ((c0 * L0 * l0).sum(-1) + (c1 * L1 * l1).sum(-1))
    denom = ((al1 * c0[:, 0]) + (1 * c1[:, 0]) + sig)

    return num/denom


def circularHist(ax, x, bins=16, density=True, offset=0, gaps=True):
    """
    Produce a circular histogram of angles on ax.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.PolarAxesSubplot
        axis instance created with subplot_kw=dict(projection='polar').

    x : array
        Angles to plot, expected in units of radians.

    bins : int, optional
        Defines the number of equal-width bins in the range. The default is 16.

    density : bool, optional
        If True plot frequency proportional to area. If False plot frequency
        proportional to radius. The default is True.

    offset : float, optional
        Sets the offset for the location of the 0 direction in units of
        radians. The default is 0.

    gaps : bool, optional
        Whether to allow gaps between bins. When gaps = False the bins are
        forced to partition the entire [-pi, pi] range. The default is True.

    Returns
    -------
    n : array or list of arrays
        The number of values in each bin.

    bins : array
        The edges of the bins.

    patches : `.BarContainer` or list of a single `.Polygon`
        Container of individual artists used to create the histogram
        or list of such containers if there are multiple input datasets.
    """
    # Wrap angles to [-pi, pi)
    x = (x + np.pi) % (2*np.pi) - np.pi

    # Force bins to partition entire circle
    if not gaps:
        bins = np.linspace(-np.pi, np.pi, num=bins+1)

    # Bin data and record counts
    n, bins = np.histogram(x, bins=bins)

    # Compute width of each bin
    widths = np.diff(bins)

    # By default plot frequency proportional to area
    if density:
        # Area to assign each bin
        area = n / x.size
        # Calculate corresponding bin radius
        radius = (area/np.pi) ** .5
    # Otherwise plot frequency proportional to radius
    else:
        radius = n

    # Plot data on ax
    patches = ax.bar(bins[:-1], radius, zorder=1, align='edge', width=widths,
                     edgecolor='C0', fill=False, linewidth=1)

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels for area plots (they are mostly obstructive)
    if density:
        ax.set_yticks([])

    return n, bins, patches


def generateTwoCorrArrays(numSamples, corr):
    """
    this function will generate a two arrays of length numSamples
    that have an inherent correlation

    Inputs: numSamples: the length of the arrays
            corr: the desired correlation b/w the two arrays
    Outputs: arr1, arr2: two correlated arrays
    """

    # Random samples (Uniformly distributed)
    U1 = np.random.rand(numSamples, 1)
    U2 = np.random.rand(numSamples, 1)

    # Random samples (normally distributed uncorrelated)
    S1 = np.sqrt(-2*np.log(U1))*np.cos(2*np.pi*U2)
    S2 = np.sqrt(-2*np.log(U1))*np.sin(2*np.pi*U2)

    # Correlated random samples
    mu_x = 0.5
    mu_y = 0.66
    sigma_x = 0.85
    sigma_y = 1.24
    rho = corr
    x = mu_x + sigma_x * S1
    y = mu_y + sigma_y * (rho*S1 + np.sqrt(1-rho**2)*S2)

    return x.flatten(), y.flatten()


def pairCorrExclude3SD(n1SpikeMat, n2SpikeMat):
    """
    this function will return the pearson's correlation coefficient
    between a pair of neuron's spike counts to a stimulus condition.
    will also return the covariance and SD (numerator/denominator of the
    correlation equation). This function will z-score the spike counts
    for each neuron and exclude trials that are larger than abs(3).

    Inputs: n1SpikeMat/n2SpikeMat: the spike counts for each neuron (array)
    Outputs: pairCorr: pearson's correlation coeff (int)
             pairCov: covariance

    """

    skipTrials = []
    n1Zscore = stats.zscore(n1SpikeMat)
    n2Zscore = stats.zscore(n2SpikeMat)
    n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
    n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
    for x in n1SkipTrials:
        skipTrials.append(x)
    for x in n2SkipTrials:
        skipTrials.append(x)
    goodTrials = [x for x in range(len(n1SpikeMat)) if x not in skipTrials]
    pairStimCorr = stats.pearsonr(n1SpikeMat[goodTrials],
                                  n2SpikeMat[goodTrials])[0]
    pairDCov = np.cov(n1SpikeMat[goodTrials],
                      n2SpikeMat[goodTrials], ddof=1)
    pairDSD = (np.std(n1SpikeMat[goodTrials], ddof=1) *
               np.std(n2SpikeMat[goodTrials], ddof=1))

    return pairStimCorr, pairDCov, pairDSD


def reIndexedRespMat(b, reIndex):
    """
    this function will take the 6x6 (nxn) stimMat matrix for each neuron
    and reIndex it so that the null direction for that neuron is in the top
    left corner
    Inputs: b (2D np.array) - the stimMat matrix for that neuron
            reIndex (np.array) - array that gives reindex operation
    Outputs: bReIndex (2D np.array) - reIndexed stimMat matrix
    """

    bReIndex = np.zeros((7, 7))
    tempMain = b[:6, :6][:, reIndex]
    tempMain = tempMain[:6, :6][reIndex, :]
    temp0Blank = b[:6, 6][reIndex]
    temp1Blank = b[6, :6][reIndex]
    bReIndex[:6, :6] = tempMain
    bReIndex[:6, 6] = temp0Blank
    bReIndex[6, :6] = temp1Blank
    bReIndex[6, 6] = b[6, 6]

    return bReIndex


def reIndexedStimMat(reIndex):
    """
    this function will return the 6x6 grid of stimMat indices as well
    as the reIndexed version of the stimMat indices.
    Inputs: reIndex (np.array) - array that gives reindex operation
    Outputs: stimMat (2D np.array) - indices of stimMat 6x6
             stimMatReIndex (2D np.array) - indices of stimMat reIndex 6x6
    """

    # fixed (independent) variables - matrix of corresponding stim Indexes
    stimMat = np.zeros((7, 7))
    stimMat[:6, :6] = np.arange(36).reshape(6, 6)
    stimMat[6, :6] = np.arange(36, 42)
    stimMat[:, 6] = np.arange(42, 49)

    # reshape fixed variables to match dependent variables (reIndexed)
    stimMatReIndex = np.zeros((7, 7))
    tempMain = stimMat[:6, :6][:, reIndex]
    tempMain = np.squeeze(tempMain[:6, :6][reIndex, :])
    temp0Blank = stimMat[:6, 6][reIndex]
    temp1Blank = stimMat[6, :6][reIndex]
    stimMatReIndex[:6, :6] = tempMain
    stimMatReIndex[:6, 6] = temp0Blank
    stimMatReIndex[6, :6] = temp1Blank
    stimMatReIndex[6, 6] = stimMat[6, 6]

    return stimMat, stimMatReIndex


def getSelAndSuppIndx(loc0Resp, loc1Resp, al0, al1):
    """
    this function will compute the selectivity and suppression index
    for the neuron of interest to different pairings of stimuli
    Inputs: loc0Resp (int): the response to loc0 stimuli from normalization fit
            loc1Resp (int): the response to loc1 stimuli from normalization fit
            al0 (int): the suppression at loc0 (alpha 0 from norm fit)
                            this is usually 1
            al1 (int): the suppression at loc1 (alpha 1 from norm fit)
    Outputs: n1Selectivity (int): selectivity index for neuron
             n1NonPrefSupp (int): suppression index for neuron
    """

    selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
    if selectivity >= 0:
        nonPrefSupp = al1 / (al0 + al1)
    else:
        nonPrefSupp = al0 / (al0 + al1)

    return selectivity, nonPrefSupp


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# END HERE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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

#
# a = meanSpike[unitCount][36:42] * 1000/trueStimDurMS
# b = np.concatenate((a[3:],a[:4]), axis=0)
# degList = np.radians(np.arange(-180,240,60) % 360)
# popt, pcov = curve_fit(vonMises, degList, b)





