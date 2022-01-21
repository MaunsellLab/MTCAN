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
        amp = np.random.randint(15,30)
        y_translate = np.random.randint(30,50)
        x_translate = np.random.randint(60,120)
        tuningMat[i,:] = (amp * np.sin((tuningMat[0,:] * np.pi / 180) + x_translate)) + y_translate
    return tuningMat


tuningCurves = randTuningCurve(1)
tc_dict = {tuningCurves[0,i]: tuningCurves[1,i] for i in range(tuningCurves.shape[1])}

allTrialsData, header = loadMatFile('Meetz_2021_1028_1.mat')

spikeCountMat = np.zeros((30,169))
spikeCountMat[0] = np.arange(0,169)
stimIndexCount = {}

# not tested, but code will add stimuli presentations to index matrix, excluding padding
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

                fRateLocO = tc_dict[stimIndexDict[index][0]['direction']]
                fRateLoc1 = tc_dict[stimIndexDict[index][1]['direction']]
                C0 = stimIndexDict[index][0]['contrast']
                C1 = stimIndexDict[index][1]['contrast']
                sigma = 0.3

                spikeTrain0 = []
                spikeTrain1 = []

                dt = 1/1000
                for i in range(500):
                    if np.random.uniform(0,1) < fRateLoc0 * dt:
                        spikeTrain0.append(1)
                    if np.random.uniform(0,1) < fRateLoc1 * dt:
                        spikeTrain1.append(1)
                L0 = len(spikeTrain0)
                L1 = len(spikeTrain1)
            
                RFSpikes = ((C0*L0) + (C1*L1))/(C0+C1+sigma)

                spikeCountMat[stimIndexCount[index]][index] = RFSpikes

# generates a dictionary of stim Index and corresponding directions/contrasts
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
                        stimIndexDict[index][stim['stimLoc']] = {'direction': stim['directionDeg'], 'contrast': stim['contrast']}

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

#poisson spike train for stimulus, can turn this into a function 

fRateLocO = tc_dict[stimIndexDict[index][0]['direction']]
fRateLoc1 = tc_dict[stimIndexDict[index][1]['direction']]
C0 = stimIndexDict[index][0]['contrast']
C1 = stimIndexDict[index][1]['contrast']

spikeTrain0 = []
spikeTrain1 = []

dt = 1/1000
for i in range(500):
    if np.random.uniform(0,1) < fRateLoc0 * dt:
        spikeTrain0.append(1)
    if np.random.uniform(0,1) < fRateLoc1 * dt:
        spikeTrain1.append(1)
L0 = len(spikeTrain0)
L1 = len(spikeTrain1)

RFSpikes = ((C0*L0) + (C1*L1))/(C0+C1+sigma)
