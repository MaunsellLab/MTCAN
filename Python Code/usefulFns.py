import scipy.io as sp
import numpy as np
import matplotlib.pyplot as plt
import math
import os
from collections import defaultdict


def loadMatFile(fileName):
    '''
    Loads the given matfile and assigns variables to access trial data

    Inputs: matfile name, (str)
    Outputs: variables, (nd.array)
    '''

    global allTrials, allTrialsData, header
    
    allTrials = sp.loadmat(f'../Matlab Data/{fileName}', squeeze_me = True)
    allTrialsData = allTrials['trials']
    header = allTrials['header']


def fieldInTrial(trial, fieldList):
    '''
    Function will check whether all fields in a list are in the trial
    
    Inputs: trial (data struct from MATLAB)
            list of fields (list)
    Outputs: bool
    '''

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