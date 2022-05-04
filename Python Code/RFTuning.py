'''
RFTuning generates heatmaps and histogram of each unit's RF location.
The script will save the plots for each unit as a PDF in a folder specific
to the day's dataset. 

to do:
incorp frame render to align spikes

Chery March 2022

Modified to save plots as pngs and incorporated changes to track stimCounts
on 04/26/22
'''

from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import seaborn as sns


allTrials, header = loadMatFile73('testing_220310_Heatmap_GRF_Spikes.mat')

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

numAzi = np.int32(header['map0Settings']['data']['azimuthDeg']['n'])
minAzi = np.int32(header['map0Settings']['data']['azimuthDeg']['minValue'])
maxAzi = np.int32(header['map0Settings']['data']['azimuthDeg']['maxValue'])
numEle = np.int32(header['map0Settings']['data']['elevationDeg']['n'])
minEle = np.int32(header['map0Settings']['data']['elevationDeg']['minValue'])
maxEle = np.int32(header['map0Settings']['data']['elevationDeg']['maxValue'])
stimDurMS = np.int32(header['mapStimDurationMS']['data'])
histPrePostMS = 50

for unit in units:
    
    stimCount = np.zeros((numEle,numAzi))
    spikeCountMat = np.zeros((50,numEle, numAzi))
    spikeCountMat[:,:,:] = np.nan
    spikeHists = np.zeros((stimDurMS + 2*histPrePostMS, numEle, numAzi))

    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        
        # stimCount verify 
        if 'numMap0Stim' in currTrial:
            map0StimLim = int(currTrial['numMap0Stim']['data'].tolist())
            map0Count = 0
            trialStartMS = currTrial['trialStart']['timeMS']
            trialStartSNEV = np.around(currTrial['taskEvents']['trialStart']['timeS'], 3)
            stimDesc = currTrial['stimDesc']
            spikeData = currTrial['spikeData']
            for sCount, stim in enumerate(stimDesc['data']):
                if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                    aziIndex = int(stim['azimuthIndex'])
                    eleIndex = int(stim['elevationIndex'])
                    stCount = int(stimCount[eleIndex][aziIndex])
                    stimOnTimeMS = stimDesc['timeMS'][sCount]
                    stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                    stimOnSNEV = trialStartSNEV + stimDiffS
                    map0Count += 1
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)
                        unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnSNEV) & 
                                        (unitTimeStamps <= stimOnSNEV + stimDurMS/1000))
                        spikeCountMat[stCount][eleIndex][aziIndex] = len(stimSpikes[0])
                        stimCount[eleIndex][aziIndex] = stimCount[eleIndex][aziIndex] + 1
                        
                        #histograms
                        # histSpikes = np.arange(stimOnSNEV - 0.050, stimOnSNEV + \
                        #                     (stimDurMS+49)/1000, 0.001)
                        # for histCount, i in enumerate(range(len(histSpikes))):
                        #     if np.around(histSpikes[i],3) in unitTimeStamps:
                        #         spikeHists[histCount][eleIndex][aziIndex] += 1
                        
                        stimOnPreSNEV = stimOnSNEV - 0.050
                        stimOnPostSNEV = stimOnSNEV + (stimDurMS+49)/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)\
                                            & (unitTimeStamps <= stimOnPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        spikeHists[histStimSpikes, eleIndex, aziIndex] += 1



    spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)

    #plot figure
    fig = plt.figure()
    fig.set_size_inches(6,8)
    aziLabel = np.linspace(minAzi,maxAzi, numAzi)
    eleLabel = np.linspace(minEle,maxEle,numEle)

    date = header['date']
    text = fig.text(0.05, 0.9, f'RF tuning for unit {unit}\n{date}', size=13,\
            fontweight='bold')
    text.set_path_effects([path_effects.Normal()])

    ax_row1 = []
    ax_row1.append(plt.subplot2grid((9,6), (0,3), colspan = 3, rowspan = 3)) # ax2
    ax_row1[0] = sns.heatmap(spikeCountMean)
    ax_row1[0].set_xlabel('azimith (deg˚)', fontsize=8)
    ax_row1[0].set_ylabel('elevation (deg˚)', fontsize = 8)
    ax_row1[0].set_title('Heatmap of unit RF location', fontsize=9)
    ax_row1[0].set_yticklabels(eleLabel, fontsize=5)
    ax_row1[0].set_xticklabels(aziLabel, fontsize=5)

    ax_row2 = []
    for i in range(3, 9):
        ax = []
        for j in range(0, 6):
            ax.append(plt.subplot2grid((9,6), (i,j)))
        ax_row2.append(np.array(ax))

    ax_row2 = np.array(ax_row2) # 6 x 6

    for i in range(numEle):
        for j in range(numAzi):
            spikeHist = spikeHists[:,i,j] * 1000/stimCount[i,j]
            histSmooth = smooth(spikeHist,75)
            ax_row2[i,j].plot(histSmooth)
            ax_row2[i,j].set_title(f"{eleLabel[i]},{aziLabel[j]}", fontsize=4)
            ax_row2[i,j].set_ylim([0, 70])
            ax_row2[i,j].set_yticks([0,35,70])
            ax_row2[i,j].set_yticklabels([0,35,70], fontsize=5)
            ax_row2[i,j].set_xticks([50,250])
            # ax_row2[i,j].set_xticklabels([])
            ax_row2[i,j].set_xticklabels([50,250], fontsize=5)
            if i == 5 and j == 0:
                ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i,j].yaxis.set_label_coords(-0.5,1.70)
    plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0)

    # saves plot as png 
    os.makedirs('RFLoc Tuning PNGs')
    os.chdir('RFLoc Tuning PNGs/')
    plt.savefig(f'{unit}.png')
    



'''
smoothHist = savgol_filter(spikeHist, 100,3)
plt.plot(smoothHist)
plt.show()

#moving point average
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

spikeHist = spikeHists[:,1,1] * 1000/ stimCount[1,1]
histSmooth = smooth(spikeHist,75)
plt.plot(histSmooth)
plt.show()
'''

