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
#testing time to open file
# t0 = time.time(); a = mat73.loadmat('./Documents/Grad School/Maunsell Lab/Analysis/MTCAN/Matlab Data/Meetz_2022_0114_MTNAN_m73.mat', use_attrdict = True); print(time.time() - t0)