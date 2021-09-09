import scipy.py as spio
import numpy as np

# import mat file (write a function to input filename)

trial = spio.loadmat('file_name.mat', squeeze_me = True)
trialData = trial['trial']
stimDesc = trialData['stimDesc']
trialDataTrial = trialData['trial']
loc0Lst = []
loc1Lst = []
loc2Lst = []
loc3lst = []

# checks to see if current trial is instruction or not
# if not instruction, will generate targetOnset stimulus
if trialDataTrial.item()['data'].item()['instructTrial'] == 0:
    for count,d in enumerate(stimDesc.item()['data'].item()['listTypes']):
        if 2 in d:
            targetOnsetStim = count
    attendLoc = trialDataTrial.item()['data'].item()['attendLoc']
    if attendLoc == 0:
        for d in stimDesc.item()['data'].item()['stimTypes'][:targetOnsetStim]:
            if len(loc0Lst) < 9:
                loc0Lst.append(d[2:4].tolist())
                lastStimSeq = d[2:4].tolist()
            else:
                # check to see if current stim sequence follows after currStimSeq
                indexLastStim = [(count,i) for count, i in enumerate(loc0Lst) if i == lastStimSeq]
                currStimSeq = d[2:4].tolist()
                indexCurrStim = [(count,i) for count, i in enumerate(loc0Lst) if i == currStimSeq]
                if indexLastStim[0][0] == len(loc0Lst) - 1:
                    if indexCurrStim[0][0] != 0:
                        print('sequence not in correct order')
                else:
                    if indexCurrStim[0][0] != indexLastStim + 1:
                        print('sequence not in correct order')
                lastStimSeq = currStimSeq
    elif attendLoc == 1:
        for d in stimDesc.item()['data'].item()['stimTypes'][:targetOnsetStim]:
            if len(loc1Lst) < 9:
                loc1Lst.append(d[2:4].tolist())
                lastStimSeq = d[2:4].tolist()
            else:
                # check to see if current stim sequence follows after currStimSeq
                indexLastStim = [(count,i) for count, i in enumerate(loc1Lst) if i == lastStimSeq]
                currStimSeq = d[2:4].tolist()
                indexCurrStim = [(count,i) for count, i in enumerate(loc1Lst) if i == currStimSeq]
                if indexLastStim[0][0] == len(loc1Lst) - 1:
                    if indexCurrStim[0][0] != 0:
                        print('sequence not in correct order')
                else:
                    if indexCurrStim[0][0] != indexLastStim + 1:
                        print('sequence not in correct order')
                lastStimSeq = currStimSeq
    elif attendLoc == 2:
        for d in stimDesc.item()['data'].item()['stimTypes'][:targetOnsetStim]:
            if len(loc2Lst) < 9:
                loc2Lst.append(d[2:4].tolist())
                lastStimSeq = d[2:4].tolist()
            else:
                # check to see if current stim sequence follows after currStimSeq
                indexLastStim = [(count,i) for count, i in enumerate(loc2Lst) if i == lastStimSeq]
                currStimSeq = d[2:4].tolist()
                indexCurrStim = [(count,i) for count, i in enumerate(loc2Lst) if i == currStimSeq]
                if indexLastStim[0][0] == len(loc2Lst) - 1:
                    if indexCurrStim[0][0] != 0:
                        print('sequence not in correct order')
                else:
                    if indexCurrStim[0][0] != indexLastStim + 1:
                        print('sequence not in correct order')
                lastStimSeq = currStimSeq
    else
        for d in stimDesc.item()['data'].item()['stimTypes'][:targetOnsetStim]:
            if len(loc3Lst) < 9:
                loc3Lst.append(d[2:4].tolist())
                lastStimSeq = d[2:4].tolist()
            else:
                # check to see if current stim sequence follows after currStimSeq
                indexLastStim = [(count,i) for count, i in enumerate(loc3Lst) if i == lastStimSeq]
                currStimSeq = d[2:4].tolist()
                indexCurrStim = [(count,i) for count, i in enumerate(loc3Lst) if i == currStimSeq]
                if indexLastStim[0][0] == len(loc3Lst) - 1:
                    if indexCurrStim[0][0] != 0:
                        print('sequence not in correct order')
                else:
                    if indexCurrStim[0][0] != indexLastStim + 1:
                        print('sequence not in correct order')
                lastStimSeq = currStimSeq


