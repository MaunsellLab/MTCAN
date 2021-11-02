'''
This script is intended to make sure stimulus sequences and interstim times 
are frozen. 
The .mat files should have both the header and the trials fields to ensure
that the script runs smoothly. 
'''

from usefulFns import *
import argparse

parser = argparse.ArgumentParser('')
parser.add_argument('fileName', type=str)
args = parser.parse_args()

allTrialsData, header = loadMatFile(args.fileName)
locDict = {0:{'seq':[],'count':[1]*9,'interstim':[]},
           1:{'seq':[],'count':[1]*9,'interstim':[]},
           2:{'seq':[],'count':[1]*9,'interstim':[]},
           3:{'seq':[],'count':[1]*9,'interstim':[]}}
lastStim = [[],[],[],[]] 
currStim = [[],[],[],[]]


def stimSeqCheck(attendLoc, stimSeqLen):
    '''
    Function to create a stimulus seq based on attended location and this 
    function will also ensure that subsequence trials follow the stimulus
    sequence that was generated. In addition, the function will generat a list
    of interstim frames in between each stimulus config and check to see that
    the list is frozen during entire dataset.

    Inputs: attendedLoc (integer)
    Outputs: appends to locDict (dictionary of dictionaries), where 0,1,2,3 
             describe the attended location
             print statement if something is out of sequence (stim sequence or
             interstim sequence)
    '''


    frameOff = stimDesc['stimOffFrame'][0]
    for d in stimDesc[1:stimSeqLen]:
        if len(locDict[attendLoc]['seq']) < 9:
            locDict[attendLoc]['seq'].append(d['stimTypes'][:2].tolist())
            lastStim[attendLoc] = d['stimTypes'][:2].tolist()
            interstim = d['stimOnFrame'] - frameOff
            frameOff = d['stimOffFrame']
            locDict[attendLoc]['interstim'].append(interstim)
        else:
            indexLastStim = locDict[attendLoc]['seq'].index(lastStim[attendLoc])
            currStim[attendLoc] = d['stimTypes'][:2].tolist()
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
            if trial['catchTrial'] != 1:
                stimDesc = currTrial['stimDesc'].item()['data'].item()
                for count, d in enumerate(stimDesc['listTypes']):
                    if 2 in d:
                        stimSeqLen = count
                        break
                attendLoc = trial['attendLoc'].tolist()
                stimSeqCheck(attendLoc, stimSeqLen)
            else:
                stimDesc = currTrial['stimDesc'].item()['data'].item()
                stimSeqLen = len(stimDesc['listTypes']) + 1 
                attendLoc = trial['attendLoc'].tolist()
                stimSeqCheck(attendLoc, stimSeqLen)   

print(locDict)