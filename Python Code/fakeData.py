import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import scipy.io as sp
import math
import os
from collections import defaultdict
from usefulFns import *


def randTuningCurve(numNeurons):
    '''
    functon will generate random tuning cruves for x number of neurons

    Inputs:
        numNueurons (int): number of neurons 
    Outputs:
        tuningMat (2D array): matrix of tuning curve values for each
                              neuron
    '''
    
    tuningMat = np.zeros((numNeurons + 1, 6))
    tuningMat[0] = np.arange(0,360,60)

    for i in range(1, tuningMat.shape[0]):
        y_translate = np.random.randint(10,50)
        x_translate = np.random.randint(60,120)
        tuningMat[i,:] = np.sin((tuningMat[0,:] * np.pi / 180) + x_translate) + y_translate
    return tuningMat


tuningCurves = randTuningCurve(2)

allTrialsData, header = loadMatFile('Meetz_2021_1028_1.mat')

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

# not tested, but code will add stimuli presentations to index matrix, excluding padding
# stimuli and when target appears before RF stimulus turns off

spikeCountMat = np.zeros((30,169))
spikeCountMat[0] = np.arange(0,169)
stimIndexCount = {}

for n,currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data'].item()
    trial = currTrial['trial'].item()['data'].item()
    if extendedEOT == 0 and trial['instructTrial'] != 1:
        stimDesc = currTrial['stimDesc'].item()['data'].item()
        stimCount = 0
        for n,stim in enumerate(stimDesc):
            if stim['listType'] == 2:
                targetIndex = n
                break
        targetOnFrame = stimDesc['stimOnFrame'][targetIndex]
        for stim in stimDesc:
            if stim['stimLoc'] == 0 and stim0Count == 0:
                stim0Count += 1
            elif stim['stimLoc'] == 0 and stim0Count != 0 and stim['stimOffFrame'] < targetOnFrame:
                index = stim['stimIndex']
                stimIndexCount[index] = stimIndexCount.get(index, 0) + 1
                spikeCountMat[stimIndexCount[index]][index] = 1



'''
for trialNum in correctTrials:
    currTrial = allTrialsData.item()[0][trialNum]

# count number of gabors in loc0 and loc1

stimArray = np.zeros((2,1))
for row, i in enumerate(stimDesc['stimLoc']):
    if i == 0:
        
        np.append(stimArray, stimDur ,axis=1)
        
    elif i == 1:
        countLoc1 += 1

# if countLoc0 == countLoc1:
#     stimArray = np.zeros((2,countLoc0))
# else:
#     print('Location 0 and Location 1 have different numbers of stimulu presented')

# for count, i in enumerate(stimDesc['stimLoc']):
#     if i == 0:
#         stimArray[0][]


# stimList = []
# for count, i in enumerate(stimDesc['stimLoc']):
#     if i == 0 or i == 1:
#         stimList.append(count)
'''
