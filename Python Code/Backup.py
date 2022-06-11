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


def poissonSpikeTrain(i, index, sigma):
    '''
    this function will generate a poisson spike train and return the normalized
    response for the RF for a given stimulus configuration
    Inputs:
        i: neuron number
        index: the index of the corresponding stimulus configuration
        sigma: semisaturation constant from neuron's contrast response function
    Outputs:
        RFSpikes: normalized response 
    '''
    # tc_dict = dict[i]
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

    return RFSpikes


allTrialsData, header = loadMatFile('Meetz_2021_1028_1.mat')

tuningCurves = randTuningCurve(1)

#code to iterate through tuning curves to create dict for each neuron 
tcDict = {}
for i, neuron in enumerate(tuningCurves[1:,:]):
    tcDict[i+1] = {}
    for j, dirResp in enumerate(neuron):
        tcDict[i+1][tuningCurves[0][j]] = dirResp


'''
# for multiple neurons
spikeCountMat = np.zeros((numNeurons, 30, 169))
spikeCountMat[:,0,:] = np.arange(0,169)
'''
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

                #input code to change 3D matrix for multiple neurons
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

'''
backup for old multiple neuron for loop
'''

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

                for x in range(1,numNeurons+1):
                    spikeCountMat[x-1][stimIndexCount[index]][index] = poissonSpikeTrain(x, index)
                    # fRateLocO = tcDict[x][stimIndexDict[index][0]['direction']]
                    # fRateLoc1 = tcDict[x][stimIndexDict[index][1]['direction']]
                    # C0 = stimIndexDict[index][0]['contrast']
                    # C1 = stimIndexDict[index][1]['contrast']
                    # sigma = 0.3

                    # spikeTrain0 = []
                    # spikeTrain1 = []

                    # dt = 1/1000
                    # for i in range(500):
                    #     if np.random.uniform(0,1) < fRateLoc0 * dt:
                    #         spikeTrain0.append(1)
                    #     if np.random.uniform(0,1) < fRateLoc1 * dt:
                    #         spikeTrain1.append(1)
                    # L0 = len(spikeTrain0)
                    # L1 = len(spikeTrain1)
                
                    # RFSpikes = ((C0*L0) + (C1*L1))/(C0+C1+sigma)
                    # spikeCountMat[x-1][stimIndexCount[index]][index] = RFSpikes


# direction tuning tbackup, 06/09/22

# Tuning code
units = activeUnits('spikeData', allTrials)
correctTrials = correctTrialsGRF(allTrials)

numDir = np.int32(header['map0Settings']['data']['directionDeg']['n'])
stimDurMS = np.int32(header['mapStimDurationMS']['data'])
histPrePostMS = 50
allTuningMat = np.zeros((len(units),numDir))

os.makedirs('Direction Tuning PDFs')
os.chdir('Direction Tuning PDFs/')

for uCount, unit in enumerate(units):

    stimCount = np.zeros((1,numDir))
    spikeCountMat = np.zeros((50,1, numDir))
    spikeCountMat[:,:,:] = np.nan
    spikeHists = np.zeros((stimDurMS + 2*histPrePostMS, 1, numDir))

    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        if 'numMap0Stim' in currTrial:
            map0StimLim = int(currTrial['numMap0Stim']['data'].tolist())
            map0Count = 0
            trialStartMS = currTrial['trialStart']['timeMS']
            trialStartSNEV = np.around(currTrial['taskEvents']['trialStart']['timeS'], 3)
            stimDesc = currTrial['stimDesc']
            spikeData = currTrial['spikeData']
            for sCount, stim in enumerate(stimDesc['data']):
                if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                    dirIndex = int(stim['directionIndex'])
                    stCount = int(stimCount[0][dirIndex])
                    stimOnTimeMS = stimDesc['timeMS'][sCount]
                    stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                    stimOnSNEV = trialStartSNEV + stimDiffS
                    stimCount[0][dirIndex] += 1
                    map0Count += 1
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)
                        if len(unitIndex) == 1:
                            unitTimeStamps = spikeData['timeStamp']
                        else:
                            unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnSNEV) & 
                                        (unitTimeStamps <= stimOnSNEV + stimDurMS/1000))
                        spikeCountMat[stCount][0][dirIndex] = len(stimSpikes[0])
                        
                        #histograms
                        stimOnPreSNEV = stimOnSNEV - 0.050
                        stimOnPostSNEV = stimOnSNEV + (stimDurMS+49)/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                            & (unitTimeStamps <= stimOnPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        spikeHists[histStimSpikes, 0, dirIndex] += 1


    spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)
    spikeCountSD = ma.std(ma.masked_invalid(spikeCountMat), axis = 0)
    allTuningMat[uCount] = spikeCountMean

    #plot figure
    #polar
    date = header['date']

    fig = plt.figure()
    fig.set_size_inches(6,8)

    text = fig.text(0.05, 0.9, f'Direction tuning for unit {unit}\n{date}',\
                    size=13, fontweight='bold')
    text.set_path_effects([path_effects.Normal()])

    ax_row1 = plt.subplot2grid((10,6), (0,3), colspan = 3, rowspan = 4, polar=True)
    theta = np.radians(np.arange(0,420,360/numDir))
    # theta = theta[::-1]
    r = (np.append(spikeCountMean, spikeCountMean[0][0]))*1000/stimDurMS
    err = (np.append(spikeCountSD, spikeCountSD[0][0]))*1000/stimDurMS
    ax_row1.plot(theta,r)
    ax_row1.errorbar(theta, r, yerr = err,fmt='o', ecolor = 'black', color='black')
    ax_row1.set_theta_zero_location("N")
    ax_row1.set_rmax(100)
    ax_row1.set_title('Direction Tuning Plot', fontsize=8)

    # hists
    spikeHistsRS = np.reshape(spikeHists, (stimDurMS + 2*histPrePostMS,2,3))
    stimCountRS = np.reshape(stimCount, (2,3))
    titleArr = np.reshape(np.arange(0,360,60),(2,3))

    ax_row2 = []
    for countI, i in enumerate(range(4, 10, 3)):
        ax = []
        for countJ, j in enumerate(range(0, 6, 2)):
            ax.append(plt.subplot2grid((10,6), (i,j), colspan = 2, rowspan = 3))
        ax_row2.append(np.array(ax))
    ax_row2 = np.array(ax_row2) # 2 x 3

    for i in range(2):
        for j in range(3):
            spikeHist = spikeHistsRS[:,i,j] * 1000/stimCountRS[i,j]
            histSmooth = smooth(spikeHist,100)
            ax_row2[i,j].plot(histSmooth)
            histTitle = titleArr[i][j]
            ax_row2[i,j].set_title(f"{histTitle}˚", fontsize=7)
            ax_row2[i,j].set_ylim([0, 100])
            ax_row2[i,j].set_yticks([0,50,100])
            ax_row2[i,j].set_yticklabels([0,50,100], fontsize=5)
            ax_row2[i,j].set_xticks([50,50+stimDurMS])
            ax_row2[i,j].set_xticklabels([50,50+stimDurMS], fontsize=5)
            if i == 1 and j == 0:
                ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i,j].yaxis.set_label_coords(-0.2,0.3)
    plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)
    
    # saves plot as pdf 
    plt.savefig(f'{unit}.pdf')

# np.save('unitsDirTuningMat', allTuningMat)