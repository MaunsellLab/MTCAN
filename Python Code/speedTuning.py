'''
to do: 
add hists (depends on number of speeds I want to test for)
 
'''

from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

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


# load data
allTrials, header = loadMatFile73('Meetz', '220607', 'Meetz_220607_GRF3_Spikes.mat')

# create folder and change dir to save PDF's and np.array
if not os.path.exists('Speed Tuning'):
    os.makedirs('Speed Tuning')
os.chdir('Speed Tuning/')


# Tuning code
units = activeUnits('spikeData', allTrials)
correctTrials = correctTrialsGRF(allTrials)

frameRateHz = header['displayCalibration']['data']['frameRateHz'].tolist()
numSpeeds = np.int32(header['map0Settings']['data']['temporalFreqHz']['n'])
minSpeed = np.int32(header['map0Settings']['data']['temporalFreqHz']['minValue'])
maxSpeed = np.int32(header['map0Settings']['data']['temporalFreqHz']['maxValue'])
stimDurMS = np.int32(header['mapStimDurationMS']['data'])
histPrePostMS = 100
allTuningMat = np.zeros((len(units),numSpeeds))

# assert frame consistency for frame duration
stimDurFrame = []
for corrTrial in correctTrials:
    currTrial = allTrials[corrTrial]
    stimDesc = currTrial['stimDesc']['data']
    if 'numMap0Stim' in currTrial:
        map0StimLim = int(currTrial['numMap0Stim']['data'].tolist())
        map0Count = 0
        for stim in stimDesc:
            if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                frameDiff = stim['stimOffFrame'].tolist() - stim['stimOnFrame'].tolist()
                stimDurFrame.append(frameDiff)
if len(set(stimDurFrame)) != 1:
    print('stimulus frame duration not consistent for mapping stimuli')
else: 
    trueStimDurMS = np.around(1000/frameRateHz * stimDurFrame[0])


for uCount, unit in enumerate(units):

    stimCount = np.zeros((1,numSpeeds+1))
    spikeCountMat = np.zeros((50,1, numSpeeds))
    spikeCountMat[:,:,:] = np.nan
    spikeHists = np.zeros((numSpeeds+1, stimDurMS + 2*histPrePostMS+12))

    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        if 'numMap0Stim' in currTrial:
            map0StimLim = int(currTrial['numMap0Stim']['data'].tolist())
            map0Count = 0
            stimDesc = currTrial['stimDesc']['data']
            spikeData = currTrial['spikeData']
            stim1TimeS = currTrial['taskEvents']['stimulusOn'][0]['time'].tolist()
            for stim in stimDesc:
                if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                    speedIndex = int(stim['temporalFreqIndex'])
                    stCount = int(stimCount[0][speedIndex])
                    stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'].tolist())
                                   /1000) + stim1TimeS
                    stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'].tolist())
                                   /1000) + stim1TimeS
                    stimCount[0][speedIndex] += 1
                    map0Count += 1
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)[0]
                        unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnTimeS) & 
                                        (unitTimeStamps <= stimOffTimeS))
                        spikeCountMat[stCount][0][speedIndex] = len(stimSpikes[0])
                        
                        #histograms
                        stimOnPreSNEV = stimOnTimeS - histPrePostMS/1000
                        stimOffPostSNEV = stimOffTimeS + histPrePostMS/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                        & (unitTimeStamps <= stimOffPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        spikeHists[speedIndex, histStimSpikes] += 1

    spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)
    spikeCountSD = ma.std(ma.masked_invalid(spikeCountMat), axis = 0)
    allTuningMat[uCount] = spikeCountMean

    #plot figure
    date = header['date']

    fig = plt.figure()
    fig.set_size_inches(6,8)

    text = fig.text(0.05, 0.85, f'Speed tuning for unit {unit}\n{date}\n- - - - -\n\
    Stimulus Duration = {stimDurMS} ms\nNumber of Blocks = {int(stimCount[0][0])}',\
                    size=10, fontweight='bold')
    text.set_path_effects([path_effects.Normal()])

    ax_row1 = plt.subplot2grid((10,6), (0,3), colspan=3, rowspan=4)
    x = np.linspace(minSpeed,maxSpeed, numSpeeds)
    # x = np.logspace(1,5,num=6,base=2)
    ax_row1.plot(x,spikeCountMean[0]*1000/stimDurMS)
    ax_row1.errorbar(x, spikeCountMean[0]*1000/stimDurMS, yerr = spikeCountSD[0]*1000/stimDurMS,fmt='o', ecolor = 'black', color='black')
    ax_row1.set_title('Speed Tuning Plot', fontsize=8)
    ax_row1.set_xlabel('Temporal Frequency (Hz)', fontsize = 8)
    ax_row1.set_ylabel('Firing Rate (spikes/sec)', fontsize = 8)

    # hists


    ax_row2 = []
    for countI, i in enumerate(range(4, 10, 3)):
        ax = []
        for countJ, j in enumerate(range(0, 6, 2)):
            ax.append(plt.subplot2grid((10,6), (i,j), colspan = 2, rowspan = 3))
        ax_row2.append(np.array(ax))
    ax_row2 = np.array(ax_row2) # 2 x 3

    yMax = 0
    plotCount = 0
    for i in range(2):
        for j in range(3):
            spikeHist = spikeHists[plotCount,:] * 1000/stimCount[0][plotCount]
            plotCount += 1
            gaussSmooth = gaussian_filter1d(spikeHist, 15)
            if max(gaussSmooth) > yMax:
                yMax = max(gaussSmooth)
            ax_row2[i,j].plot(gaussSmooth)
            ax_row2[i,j].set_title(f"{i},{j}", fontsize=7)
            ax_row2[i,j].set_ylim(bottom=0)
            ax_row2[i,j].set_yticks([0,50,100])
            ax_row2[i,j].set_yticklabels([0,50,100], fontsize=5)
            ax_row2[i,j].set_xticks([0,histPrePostMS,histPrePostMS+stimDurMS,2*histPrePostMS+stimDurMS])
            ax_row2[i,j].set_xticklabels([-(histPrePostMS), 0, 0+stimDurMS, stimDurMS+histPrePostMS], fontsize=5)
            ax_row2[i,j].axvspan(histPrePostMS, histPrePostMS+stimDurMS, color='grey', alpha=0.2)
            if i == 1 and j == 0:
                ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i,j].yaxis.set_label_coords(-0.2,0.3)
    plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)

    for i in range(2):
        for j in range(3):
            ax_row2[i,j].set_ylim([0, yMax*1.1])
    # saves plot as pdf 
    plt.savefig(f'{unit}.pdf')

np.save('unitsSpeedTuningMat', allTuningMat)
plt.close('all')