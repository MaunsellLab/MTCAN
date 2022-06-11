'''
to do: 
add hists (depends on number of speeds I want to test for)
incorp frame render to align spikes 
'''

from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

allTrials, header = loadMatFile73('Meetz', '220607', 'Meetz_220607_GRF3_Spikes.mat')

# for testing purposes, to make unit field similar to real data
for currTrial in allTrials:
    if 'spikeData' in currTrial:
        currTrial['spikeData']['unit'] = currTrial['spikeData']['unit'].tolist()
        for i in range(0,len(currTrial['spikeData']['channel'])):
            a = str(int(currTrial['spikeData']['channel'][i]))
            b = str(int(currTrial['spikeData']['unit'][i]))
            c = a + '_' + b
            currTrial['spikeData']['unit'][i] = c
        currTrial['spikeData']['unit'] = np.array(currTrial['spikeData']['unit'])


# Tuning code
units = activeUnits('spikeData', allTrials)
correctTrials = correctTrialsGRF(allTrials)

numSpeeds = np.int32(header['map0Settings']['data']['temporalFreqHz']['n'])
minSpeed = np.int32(header['map0Settings']['data']['temporalFreqHz']['minValue'])
maxSpeed = np.int32(header['map0Settings']['data']['temporalFreqHz']['maxValue'])
stimDurMS = np.int32(header['mapStimDurationMS']['data'])
histPrePostMS = 50
allTuningMat = np.zeros((len(units),numSpeeds))

os.makedirs('Speed Tuning')
os.chdir('Speed Tuning/')

for uCount, unit in enumerate(units):

    stimCount = np.zeros((1,numSpeeds))
    spikeCountMat = np.zeros((50,1, numSpeeds))
    spikeCountMat[:,:,:] = np.nan
    spikeHists = np.zeros((stimDurMS + 2*histPrePostMS, 1, numSpeeds))

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
                    speedIndex = int(stim['temporalFreqIndex'])
                    stCount = int(stimCount[0][speedIndex])
                    stimOnTimeMS = stimDesc['timeMS'][sCount]
                    stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                    stimOnSNEV = trialStartSNEV + stimDiffS
                    stimCount[0][speedIndex] += 1
                    map0Count += 1
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)
                        if len(unitIndex) == 1:
                            unitTimeStamps = spikeData['timeStamp']
                        else:
                            unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnSNEV) & 
                                        (unitTimeStamps <= stimOnSNEV + stimDurMS/1000))
                        spikeCountMat[stCount][0][speedIndex] = len(stimSpikes[0])
                        
                        #histograms
                        stimOnPreSNEV = stimOnSNEV - 0.050
                        stimOnPostSNEV = stimOnSNEV + (stimDurMS+49)/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                        & (unitTimeStamps <= stimOnPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        spikeHists[histStimSpikes, 0, speedIndex] += 1

    spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)
    spikeCountSD = ma.std(ma.masked_invalid(spikeCountMat), axis = 0)
    allTuningMat[uCount] = spikeCountMean

    #plot figure
    date = header['date']

    fig = plt.figure()
    fig.set_size_inches(6,8)

    text = fig.text(0.05, 0.9, f'Speed tuning for unit {unit}\n{date}',\
                    size=13, fontweight='bold')
    text.set_path_effects([path_effects.Normal()])

    ax_row1 = plt.subplot2grid((10,6), (0,3), colspan=3, rowspan=4)
    x = np.linspace(minSpeed,maxSpeed, numSpeeds)
    ax_row1.plot(x,spikeCountMean[0]*1000/stimDurMS)
    ax_row1.errorbar(x, spikeCountMean[0]*1000/stimDurMS, yerr = spikeCountSD[0]*1000/stimDurMS,fmt='o', ecolor = 'black', color='black')
    ax_row1.set_title('Speed Tuning Plot', fontsize=8)
    ax_row1.set_xlabel('Temporal Frequency (Hz)', fontsize = 8)
    ax_row1.set_ylabel('Firing Rate (spikes/sec)', fontsize = 8)

    # hists
    # spikeHistsRS = np.reshape(spikeHists, (stimDurMS + 2*histPrePostMS,2,3))
    # stimCountRS = np.reshape(stimCount, (2,3))

    # ax_row2 = []
    # for countI, i in enumerate(range(4, 10, 3)):
    #     ax = []
    #     for countJ, j in enumerate(range(0, 6, 2)):
    #         ax.append(plt.subplot2grid((10,6), (i,j), colspan = 2, rowspan = 3))
    #     ax_row2.append(np.array(ax))
    # ax_row2 = np.array(ax_row2) # 2 x 3

    # for i in range(2):
    #     for j in range(3):
    #         spikeHist = spikeHistsRS[:,i,j] * 1000/stimCountRS[i,j]
    #         histSmooth = smooth(spikeHist,75)
    #         ax_row2[i,j].plot(histSmooth)
    #         ax_row2[i,j].set_title(f"{i},{j}", fontsize=7)
    #         ax_row2[i,j].set_ylim([0, 70])
    #         ax_row2[i,j].set_yticks([0,35,70])
    #         ax_row2[i,j].set_yticklabels([0,35,70], fontsize=5)
    #         ax_row2[i,j].set_xticks([50,50+stimDurMS])
    #         ax_row2[i,j].set_xticklabels([50,50+stimDurMS], fontsize=5)
    #         if i == 1 and j == 0:
    #             ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
    #             ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
    #             ax_row2[i,j].yaxis.set_label_coords(-0.2,0.3)
    # plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)
    
    # saves plot as pdf 
    plt.savefig(f'{unit}.pdf')

np.save('unitsSpeedTuningMat', allTuningMat)