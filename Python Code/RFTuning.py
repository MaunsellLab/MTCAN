'''
things to do
1. change heatmap axis to include min and max azi/ele used
2. ^for PSTHs
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

for unit in units:

    numAzi = int(header['mapSettings']['data']['azimuthDeg']['n'].tolist())
    numEle = int(header['mapSettings']['data']['elevationDeg']['n'].tolist())
    stimDurMS = int(header['mapStimDurationMS']['data'].tolist())
    histPrePostMS = 50

    stimCount = np.zeros((numEle,numAzi))
    spikeCountMat = np.zeros((50,numEle, numAzi))
    spikeCountMat[:,:,:] = np.nan
    spikeHists = np.zeros((stimDurMS + 2*histPrePostMS, numEle, numAzi))


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
                    aziIndex = int(stim['azimuthIndex'])
                    eleIndex = int(stim['elevationIndex'])
                    stCount = int(stimCount[eleIndex][aziIndex])
                    stimOnTimeMS = stimDesc['timeMS'][sCount]
                    stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                    stimOnSNEV = trialStartSNEV + stimDiffS
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)
                        unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnSNEV) & 
                                        (unitTimeStamps <= stimOnSNEV + stimDurMS/1000))
                        spikeCountMat[stCount][eleIndex][aziIndex] = len(stimSpikes[0])
                        stimCount[eleIndex][aziIndex] = stimCount[eleIndex][aziIndex] + 1
                        
                        #histograms
                        histSpikes = np.arange(stimOnSNEV - 0.050, stimOnSNEV + \
                                            (stimDurMS+49)/1000, 0.001)
                        for histCount, i in enumerate(range(len(histSpikes))):
                            if np.around(histSpikes[i],3) in unitTimeStamps:
                                spikeHists[histCount][eleIndex][aziIndex] += 1


    spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)

    #plot figure
    fig = plt.figure()
    fig.set_size_inches(6,8)

    date = header['date']
    text = fig.text(0.05, 0.9, f'RF tuning for unit {unit}\n{date}', size=13)
    text.set_path_effects([path_effects.Normal()])

    ax_row1 = []
    ax_row1.append(plt.subplot2grid((10,6), (0,3), colspan = 3, rowspan = 4)) # ax2
    ax_row1[0] = sns.heatmap(spikeCountMean)
    ax_row1[0].set_xlabel('azimith (deg˚)', fontsize=8)
    ax_row1[0].set_ylabel('elevation (deg˚)', fontsize = 8)
    ax_row1[0].set_title('Heatmap of unit RF location', fontsize=9)

    ax_row2 = []
    for i in range(4, 10):
        ax = []
        for j in range(0, 6):
            ax.append(plt.subplot2grid((10,6), (i,j)))
        ax_row2.append(np.array(ax))

    ax_row2 = np.array(ax_row2) # 6 x 6

    for i in range(numEle):
        for j in range(numAzi):
            spikeHist = spikeHists[:,i,j] * 1000/stimCount[i,j]
            histSmooth = smooth(spikeHist,75)
            ax_row2[i,j].plot(histSmooth)
            ax_row2[i,j].set_title(f"{i},{j}", fontsize=7)
            ax_row2[i,j].set_ylim([0, 70])
            ax_row2[i,j].set_yticks([0,35,70])
            ax_row2[i,j].set_yticklabels([0,35,70], fontsize=5)
            ax_row2[i,j].set_xticks([50,250])
            ax_row2[i,j].set_xticklabels([50,250], fontsize=5)
            if i == 5 and j == 0:
                ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i,j].yaxis.set_label_coords(-0.5,1.70)
    plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)
    plt.show() 
    #plt save as pdf


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

