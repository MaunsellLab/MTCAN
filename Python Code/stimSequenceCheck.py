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
    stimDesc = currTrial['stimDesc']
    trial = currTrial['trial']
    trialEnd = currTrial['trialEnd']

    if trialEnd.item()['data'] == 0: 
        if trial.item()['data'].item()['instructTrial'] == 0:
            for count, d in enumerate(stimDesc.item()['data'].item()['listTypes']):
                if 2 in d:
                    targetOnsetStim = count
            attendLoc = trial.item()['data'].item()['attendLoc']
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
            indexLastStim = [(count,i) for count, i in enumerate(locLsts[attendLoc]) if i == lastStim[attendLoc]]
            currStim[attendLoc] = d[2:4].tolist()
            indexCurrStim = [(count,i) for count, i in enumerate(loc0Lst[attendLoc]) if i == currStim[attendLoc]]
            if indexLastStim[0][0] == len(locLsts[attendLoc]) - 1:
                if indexCurrStim[0][0] != 0:
                    print('start sequence not in correct order')
            else:
                if indexCurrStim[0][0] != indexLastStim[0][0] + 1:
                    print('middle sequence not in correct order')
            lastStim[attendLoc] = currStim[attendLoc]
