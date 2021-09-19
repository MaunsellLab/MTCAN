


import scipy.io as sp
import numpy as np


#importing mat file with multiple trials
allTrials = sp.loadmat('xxx.mat', squeeze_me = True)
allTrialsData = allTrials['trials']

loc0Lst = []
loc1Lst = []
loc2Lst = []
loc3Lst = []

#to run through diff trials 
for i in range(0,len(allTrialsData.item()[0])-1):
    currTrial = allTrialsData.item()[0][i]
    stimDesc = currTrial['stimDesc']
    trial = currTrial['trial']
    trialEnd = currTrial['trialEnd']


    if trialEnd.item()['data'] == 0: 
        if trial.item()['data'].item()['instructTrial'] == 0:
            for count, d in enumerate(stimDesc.item()['data'].item()['listTypes']):
                if 2 in d:
                    targetOnsetStim = count
            attendLoc = trial.item()['data'].item()['attendLoc']
            if attendLoc == 0:
                for d in stimDesc.item()['data'].item()['stimTypes'][1:targetOnsetStim]:
                    if len(loc0Lst) < 9:
                        loc0Lst.append(d[2:4].tolist())
                        lastStimSeq = d[2:4].tolist()
                    else:
                        indexLastStim = [(count,i) for count, i in enumerate(loc0Lst) if i == lastStimSeq]
                        currStimSeq = d[2:4].tolist()
                        indexCurrStim = [(count,i) for count, i in enumerate(loc0Lst) if i == currStimSeq]
                        if indexLastStim[0][0] == len(loc0Lst) - 1:
                            print(indexLastStim[0][0])
                            if indexCurrStim[0][0] != 0:
                                print(indexCurrStim[0][0])
                                print('start sequence not in correct order')
                        else:
                            if indexCurrStim[0][0] != indexLastStim[0][0] + 1:
                                print('not start sequence not in correct order')
                        lastStimSeq = currStimSeq
            elif attendLoc == 1:
                for d in stimDesc.item()['data'].item()['stimTypes'][1:targetOnsetStim]:
                    if len(loc1Lst) < 9:
                        loc1Lst.append(d[2:4].tolist())
                        lastStimSeq = d[2:4].tolist()
                    else:
                        indexLastStim = [(count,i) for count, i in enumerate(loc1Lst) if i == lastStimSeq]
                        currStimSeq = d[2:4].tolist()
                        indexCurrStim = [(count,i) for count, i in enumerate(loc1Lst) if i == currStimSeq]
                        if indexLastStim[0][0] == len(loc1Lst) - 1:
                            if indexCurrStim[0][0] != 0:
                                print('start sequence not in correct order')
                        else:
                            if indexCurrStim[0][0] != indexLastStim[0][0] + 1:
                                print('not start sequence not in correct order')
                        lastStimSeq = currStimSeq
            elif attendLoc == 2:
                for d in stimDesc.item()['data'].item()['stimTypes'][1:targetOnsetStim]:
                    if len(loc2Lst) < 9:
                        loc2Lst.append(d[2:4].tolist())
                        lastStimSeq = d[2:4].tolist()
                    else:
                        indexLastStim = [(count,i) for count, i in enumerate(loc2Lst) if i == lastStimSeq]
                        currStimSeq = d[2:4].tolist()
                        indexCurrStim = [(count,i) for count, i in enumerate(loc2Lst) if i == currStimSeq]
                        if indexLastStim[0][0] == len(loc2Lst) - 1:
                            if indexCurrStim[0][0] != 0:
                                print('start sequence not in correct order')
                        else:
                            if indexCurrStim[0][0] != indexLastStim[0][0] + 1:
                                print('not start sequence not in correct order')
                        lastStimSeq = currStimSeq
            else:
                for d in stimDesc.item()['data'].item()['stimTypes'][1:targetOnsetStim]:
                    if len(loc3Lst) < 9:
                        loc3Lst.append(d[2:4].tolist())
                        lastStimSeq = d[2:4].tolist()
                    else:
                        indexLastStim = [(count,i) for count, i in enumerate(loc3Lst) if i == lastStimSeq]
                        currStimSeq = d[2:4].tolist()
                        indexCurrStim = [(count,i) for count, i in enumerate(loc3Lst) if i == currStimSeq]
                        if indexLastStim[0][0] == len(loc3Lst) - 1:
                            if indexCurrStim[0][0] != 0:
                                print('start sequence not in correct order')
                        else:
                            if indexCurrStim[0][0] != indexLastStim[0][0] + 1:
                                print('not startsequence not in correct order')
                        lastStimSeq = currStimSeq
                        
# Attempt 2

import scipy.io as sp
import numpy as np

#importing mat file with multiple trials
allTrials = sp.loadmat('xxx.mat', squeeze_me = True)
allTrialsData = allTrials['trials']

locLsts = [[],[],[],[]] 
lastStim = [[],[],[],[]] 
currStim = [[],[],[],[]]

#to run through diff trials 
for i in range(0,len(allTrialsData.item()[0])-1):
    print(i)
    currTrial = allTrialsData.item()[0][i]
    trial = currTrial['trial']
    trialEnd = currTrial['trialEnd']

    if trialEnd.item()['data'] == 0: 
        if trial.item()['data'].item()['instructTrial'] == 0:
            stimDesc = currTrial['stimDesc']
            for count, d in enumerate(stimDesc.item()['data'].item()['listTypes']):
                if 2 in d:
                    targetOnsetStim = count
            attendLoc = trial.item()['data'].item()['attendLoc'].tolist()
            stimSeqCheck(attendLoc)

def stimSeqCheck(attendLoc):
    '''
    Function to create a stimulus seq based on attended location.
    Function will also ensure that subsequence trials follow the stimulus
    sequence that was generated. 

    inputs: attendedLoc (integer)
    outputs: appends to locLsts (list of lists), where 0,1,2,3 describe 
             attended loc
             print statement if something is out of sequence
    '''
    for d in stimDesc.item()['data'].item()['stimTypes'][1:targetOnsetStim]:
        if len(locLsts[attendLoc]) < 9:
            locLsts[attendLoc].append(d[2:4].tolist())
            lastStim[attendLoc] = d[2:4].tolist()
        else:
            indexLastStim = [(count,i) for count, i in enumerate(locLsts[attendLoc]) \
                if i == lastStim[attendLoc]]
            currStim[attendLoc] = d[2:4].tolist()
            indexCurrStim = [(count,i) for count, i in enumerate(locLsts[attendLoc]) \
                if i == currStim[attendLoc]]
            if indexLastStim[0][0] == len(locLsts[attendLoc]) - 1:
                if indexCurrStim[0][0] != 0:
                    print('start sequence not in correct order')
            else:
                if indexCurrStim[0][0] != indexLastStim[0][0] + 1:
                    print('middle sequence not in correct order')
            lastStim[attendLoc] = currStim[attendLoc]

# attempt 3

import scipy.io as sp
import numpy as np

# when importing .mat file, make sure it is the combined file with both the 
# header and all trials saved 
# this imports all trials within a .mat file and the accompanying header file
allTrials = sp.loadmat('Meetz_2021_08_25.mat', squeeze_me = True)
allTrialsData = allTrials['trials']
header = allTrials['header']


locDict = {0:{'seq':[],'count':[1]*9,'interstim':[]},
           1:{'seq':[],'count':[1]*9,'interstim':[]},
           2:{'seq':[],'count':[1]*9,'interstim':[]},
           3:{'seq':[],'count':[1]*9,'interstim':[]}}
lastStim = [[],[],[],[]] 
currStim = [[],[],[],[]]


def stimSeqCheck(attendLoc):
    '''
    Function to create a stimulus seq based on attended location and this 
    function will also ensure that subsequence trials follow the stimulus
    sequence that was generated. In addition, the function will generat a list
    of interstim frames in between each stimulus config and check to see that
    the list is frozen during entire dataset.

    inputs: attendedLoc (integer)
    outputs: appends to locDict (dictionary of dictionaries), where 0,1,2,3 
             describe the attended location
             print statement if something is out of sequence (stim sequence or
             interstim sequence)
    '''


    frameOff = stimDesc.item()['stimOffFrame'][0]
    for d in stimDesc.item()[1:targetOnsetStim]:
        if len(locDict[attendLoc]['seq']) < 9:
            locDict[attendLoc]['seq'].append(d['stimTypes'][2:4].tolist())
            lastStim[attendLoc] = d['stimTypes'][2:4].tolist()
            interstim = d['stimOnFrame'] - frameOff
            frameOff = d['stimOffFrame']
            locDict[attendLoc]['interstim'].append(interstim)
        else:
            indexLastStim = [(count,i) for count, i in enumerate(locDict[attendLoc]['seq']) \
                if i == lastStim[attendLoc]]
            currStim[attendLoc] = d['stimTypes'][2:4].tolist()
            indexCurrStim = [(count,i) for count, i in enumerate(locDict[attendLoc]['seq']) \
                if i == currStim[attendLoc]]
            locDict[attendLoc]['count'][indexCurrStim[0][0]] += 1
            if indexLastStim[0][0] == len(locDict[attendLoc]['seq']) - 1:
                if indexCurrStim[0][0] != 0:
                    print('start sequence not in correct order')
            else:
                if indexCurrStim[0][0] != indexLastStim[0][0] + 1:
                    print('middle sequence not in correct order')
            lastStim[attendLoc] = currStim[attendLoc]
            if d['stimOnFrame'] - frameOff != locDict[attendLoc]['interstim'][indexCurrStim[0][0]]:
                print('interstim seq not frozen')
            else:
                frameOff = d['stimOffFrame']


for currTrial in allTrialsData.item()[0]:
    trial = currTrial['trial'].item()['data'].item()
    trialEnd = currTrial['trialEnd'].item()['data']

    if trialEnd == 0: 
        if trial['instructTrial'] == 0:
            stimDesc = currTrial['stimDesc'].item()['data']
            for count, d in enumerate(stimDesc.item()['listTypes']):
                if 2 in d:
                    targetOnsetStim = count
                    break
            attendLoc = trial['attendLoc'].tolist()
            stimSeqCheck(attendLoc)