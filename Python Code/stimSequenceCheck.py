import scipy.io as sp
import numpy as np

# Import .mat file with multiple trials
allTrials = sp.loadmat('alltrials_testing_2021_0912.mat', squeeze_me = True)
allTrialsData = allTrials['trials']
header = allTrials['header']


def stimSeqCheck(attendLoc):
    '''
    Function to create a stimulus seq based on attended location.
    Function will also ensure that subsequence trials follow the stimulus
    sequence that was generated. 

    inputs: attendedLoc (integer)
    outputs: appends to locDict (dictionary of dictionaries), where 0,1,2,3 
             describe the attended location
             print statement if something is out of sequence
    '''
    for d in stimDesc.item()['data'].item()['stimTypes'][1:targetOnsetStim]:
        if len(locDict[attendLoc]['seq']) < 9:
            locDict[attendLoc]['seq'].append(d[2:4].tolist())
            lastStim[attendLoc] = d[2:4].tolist()
        else:
            indexLastStim = [(count,i) for count, i in enumerate(locDict[attendLoc]['seq']) \
                if i == lastStim[attendLoc]]
            currStim[attendLoc] = d[2:4].tolist()
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


locDict = {0:{'seq':[],'count':[1]*9},1:{'seq':[],'count':[1]*9},2:{'seq':[],'count':[1]*9},3:{'seq':[],'count':[1]*9}}
lastStim = [[],[],[],[]] 
currStim = [[],[],[],[]]

for i in range(0,len(allTrialsData.item()[0])-1):
    currTrial = allTrialsData.item()[0][i]
    trial = currTrial['trial']
    trialEnd = currTrial['trialEnd']

    if trialEnd.item()['data'] == 0: 
        if trial.item()['data'].item()['instructTrial'] == 0:
            stimDesc = currTrial['stimDesc']
            for count, d in enumerate(stimDesc.item()['data'].item()['listTypes']):
                if 2 in d:
                    targetOnsetStim = count
                    break
            attendLoc = trial.item()['data'].item()['attendLoc'].tolist()
            stimSeqCheck(attendLoc)
