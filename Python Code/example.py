


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
