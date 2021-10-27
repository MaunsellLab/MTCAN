import scipy.io as sp
import numpy as np
import matplotlib.pyplot as plt

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
            indexLastStim = locDict[attendLoc]['seq'].index(lastStim[attendLoc])
            currStim[attendLoc] = d['stimTypes'][2:4].tolist()
            indexCurrStim = locDict[attendLoc]['seq'].index(currStim[attendLoc])
            locDict[attendLoc]['count'][indexCurrStim] += 1
            if indexLastStim == len(locDict[attendLoc]['seq']) - 1:
                if indexCurrStim != 0:
                    print('start sequence not in correct order')
            else:
                if indexCurrStim != indexLastStim + 1:
                    print('middle sequence not in correct order')
            lastStim[attendLoc] = currStim[attendLoc]
            if d['stimOnFrame'] - frameOff != locDict[attendLoc]['interstim'][indexCurrStim]:
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
