from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

allTrials, header = loadMatFile73('testing_220317_Dir_GRF_Spikes.mat')

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

def activeUnits(spikeData):

    '''
    function returns the active units across trials for a session as a list

    Inputs: unitData (str): spikeData
    Outputs: units (list): active units for a sessioon

    '''
    units = []
    for currTrial in allTrials:
        if spikeData in currTrial:
            uniqueUnits = np.unique(currTrial[spikeData]['unit'])
            for unique in uniqueUnits:
                if unique not in units:
                    units.append(unique)
    
    return units

units = activeUnits('spikeData')

numDir = int(header['mapSettings']['data']['directionDeg']['n'].tolist())
stimDurMS = int(header['mapStimDurationMS']['data'].tolist())
histPrePostMS = 50

stimCount = np.zeros((1,numDir))
spikeCountMat = np.zeros((50,1, numDir))
spikeCountMat[:,:,:] = np.nan
spikeHists = np.zeros((stimDurMS + 2*histPrePostMS, 1, numDir))

for currTrial in allTrials:
    trial = currTrial['trial']['data']
    trialEnd = currTrial['trialEnd']['data']
    if trial['instructTrial'] != 1 and trialEnd == 0:
        trialStartMS = currTrial['trialStart']['timeMS']
        trialStartSNEV = np.around(currTrial['taskEvents']['trialStart']['timeS'], 3)
        stimDesc = currTrial['stimDesc']
        spikeData = currTrial['spikeData']
        for sCount, stim in enumerate(stimDesc['data']):
            if stim['stimType'] == 2:
                break
            if stim['gaborIndex'] == 1:
                dirIndex = int(stim['directionIndex'])
                stCount = int(stimCount[0][dirIndex])
                stimOnTimeMS = stimDesc['timeMS'][sCount]
                stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                stimOnSNEV = trialStartSNEV + stimDiffS
                if unit in spikeData['unit']:
                    spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                    unitIndex = np.where(spikeData['unit'] == unit)
                    unitTimeStamps = spikeData['timeStamp'][unitIndex]
                    stimSpikes = np.where((unitTimeStamps >= stimOnSNEV) & 
                                    (unitTimeStamps <= stimOnSNEV + stimDurMS/1000))
                    spikeCountMat[stCount][0][dirIndex] = len(stimSpikes[0])
                    stimCount[0][dirIndex] += 1
                    
                    #histograms
                    histSpikes = np.arange(stimOnSNEV - 0.050, stimOnSNEV + \
                                          (stimDurMS+49)/1000, 0.001)
                    for histCount, i in enumerate(range(len(histSpikes))):
                        if np.around(histSpikes[i],3) in unitTimeStamps:
                            spikeHists[histCount][0][dirIndex] += 1


spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)
spikeCountSD = ma.std(ma.masked_invalid(spikeCountMat), axis = 0)

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
r = (np.append(spikeCountMean, spikeCountMean[0][0]))*1000/stimDurMS
err = (np.append(spikeCountSD, spikeCountSD[0][0]))*1000/stimDurMS
ax_row1.plot(theta,r)
ax_row1.errorbar(theta, r, yerr = err,fmt='o', ecolor = 'black', color='black')
ax_row1.set_theta_zero_location("N")
ax_row1.set_rmax(100)
ax_row1.set_title('Direction tuning polar plot', fontsize=8)

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
        ax_row2[i,j].set_xticks([50,250])
        ax_row2[i,j].set_xticklabels([50,250], fontsize=5)
        if i == 1 and j == 0:
            ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
            ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
            ax_row2[i,j].yaxis.set_label_coords(-0.2,0.3)
plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)
plt.show()
