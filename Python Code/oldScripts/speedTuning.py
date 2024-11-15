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
        # currTrial['spikeData']['unit'] = currTrial['spikeData']['unit'].tolist()
        for i in range(0,len(currTrial['spikeData']['channel'])):
            a = str(int(currTrial['spikeData']['channel'][i]))
            b = str(int(currTrial['spikeData']['unit'][i]))
            c = a + '_' + b
            currTrial['spikeData']['unitChannel'][i] = c
        currTrial['spikeData']['unit'] = np.array(currTrial['spikeData']['unit'])

##### Start Here
# Load relevant file here with pyMat reader 

monkeyName = 'Meetz'
seshDate = '221108'
fileName = f'{monkeyName}_{seshDate}_GRF3_Spikes.mat'

allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)
# allTrials, header = loadMatFilePyMat('Meetz', '221010', 'Meetz_221010_GRF3_Spikes.mat')


# create folder and change dir to save PDF's and np.array
if not os.path.exists('Speed Tuning'):
    os.makedirs('Speed Tuning')
os.chdir('Speed Tuning/')


## Tuning code
correctTrials = correctTrialsGRF(allTrials)
units = activeUnits('spikeData', allTrials)
unitChannel = unitsInfo(units, correctTrials, allTrials)


## Assert: stimDesc field is in right format for easy access
for corrTrial in correctTrials:
    currTrial = allTrials[corrTrial]
    nStim = len(currTrial['stimDesc']['data']['stimType'])
    currTrial['stimDesc']['data'] = [{k:v[i] for k,v in currTrial['stimDesc']['data'].items()} 
                                    for i in range(nStim)]


## Declare Variables 
frameRateHz = header['displayCalibration']['data']['frameRateHz']
numTempFq = header['map0Settings']['data']['temporalFreqHz']['n']
minTempFq = np.float32(header['map0Settings']['data']['temporalFreqHz']['minValue'])
maxTempFq = np.float32(header['map0Settings']['data']['temporalFreqHz']['maxValue'])
spatialFreq = np.float32(header['map0Settings']['data']['spatialFreqCPD']['minValue'])
numDirs = header['map0Settings']['data']['directionDeg']['n']
minDir = np.int32(header['map0Settings']['data']['directionDeg']['minValue'])
maxDir = np.int32(header['map0Settings']['data']['directionDeg']['maxValue'])
interstimDurMS = header['mapInterstimDurationMS']['data']
histPrePostMS = 100
sponWindowMS = 100
numBlocks = allTrials[-1]['mappingBlockStatus']['data']['blocksDone']
allTuningMat = np.zeros((len(units),numTempFq))

# assert frame consistency for frame duration
stimDurFrame = []
for corrTrial in correctTrials:
    currTrial = allTrials[corrTrial]
    stimDesc = currTrial['stimDesc']['data']
    if 'numMap0Stim' in currTrial:
        map0StimLim = currTrial['numMap0Stim']['data']
        map0Count = 0
        for stim in stimDesc:
            if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                frameDiff = stim['stimOffFrame'] - stim['stimOnFrame']
                stimDurFrame.append(frameDiff)
if len(set(stimDurFrame)) != 1:
    print('stimulus frame duration not consistent for mapping stimuli')
else: 
    trueStimDurMS = np.int32(np.around(1000/frameRateHz * stimDurFrame[0]))


for uCount, unit in enumerate(units):

    stimCount = np.zeros((numDirs,numTempFq))
    spikeCountMat = np.zeros((numDirs,numBlocks+1, numTempFq))
    spikeHists = np.zeros((numDirs, numTempFq, trueStimDurMS+2*histPrePostMS+1))
    sponSpikesArr = []

    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        if 'numMap0Stim' in currTrial:
            map0StimLim = currTrial['numMap0Stim']['data']
            map0Count = 0
            stimDesc = currTrial['stimDesc']['data']
            spikeData = currTrial['spikeData']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            for stim in stimDesc:
                if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                    speedIndex = int(stim['temporalFreqIndex'])
                    dirIndex = int(stim['directionIndex'])
                    stCount = int(stimCount[dirIndex][speedIndex])
                    stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'])
                                   /1000) + stim1TimeS
                    stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'])
                                   /1000) + stim1TimeS
                    stimCount[dirIndex][speedIndex] += 1
                    map0Count += 1
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)[0]
                        unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        # spike count during stim presentation
                        stimSpikes = np.where((unitTimeStamps >= stimOnTimeS+0.05) & 
                                        (unitTimeStamps <= stimOffTimeS+0.075))
                        spikeCountMat[dirIndex][stCount][speedIndex] = len(stimSpikes[0])
                        # spike count during spontaneous (-100ms -> 0ms stim on) 
                        sponSpikes = np.where((unitTimeStamps >= (stimOnTimeS-(sponWindowMS/1000))) & 
                                        (unitTimeStamps <= stimOnTimeS))
                        sponSpikesArr.extend([len(sponSpikes[0])])

                        #histograms
                        stimOnPreSNEV = stimOnTimeS - histPrePostMS/1000
                        stimOffPostSNEV = stimOffTimeS + histPrePostMS/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                        & (unitTimeStamps <= stimOffPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        spikeHists[dirIndex, speedIndex, histStimSpikes] += 1

    spikeCountMean = np.mean(spikeCountMat[:,:numBlocks+1,:], axis=1)
    spikeCountSD = np.std(spikeCountMat[:,:numBlocks+1,:], axis=1)
    spikeCountSEM = spikeCountSD/np.sqrt(numBlocks)
    sponSpikesArr = np.array(sponSpikesArr)
    sponSpikesMean = np.mean(sponSpikesArr)
    sponSpikesSEM = np.std(sponSpikesArr)/np.sqrt(len(sponSpikesArr))

    # allTuningMat[uCount] = spikeCountMean

    ############# 
    # Figure
    ############# 

    date = header['date']
    tempFreqList = [minTempFq]
    a = minTempFq
    for i in range(numTempFq-1):
        a = a * 2
        tempFreqList.append(a)
    tempFreq = np.array(tempFreqList)
    speed = np.around(tempFreq/spatialFreq,2)
    dirs = np.arange(0,360,360/numDirs)
    sponList = np.array([sponSpikesMean] * len(speed))
    sponSEM = np.array([sponSpikesSEM] * len(speed))

    fig = plt.figure()
    fig.set_size_inches(6,8)
    text = fig.text(0.05, 0.85, f'Speed tuning for unit {unit}\n\
    {date}\n\
    - - - - - - - - - -\n\
    Stimulus Duration: {trueStimDurMS} ms\n\
    Number of Blocks: {numBlocks}\n\
    Interstimulus Duration: {interstimDurMS}\n\
    Channel: {unitChannel[uCount]}',size=8, fontweight='bold')
    text.set_path_effects([path_effects.Normal()])

    #### Line Graph
    ax_row1 = plt.subplot2grid((10,6), (0,3), colspan=3, rowspan=4)
    for i,dir in enumerate(dirs):
        ax_row1.plot(speed,spikeCountMean[i]*1000/trueStimDurMS, label=f'{dir}˚')
        ax_row1.errorbar(speed,spikeCountMean[i]*1000/trueStimDurMS,
        yerr = spikeCountSEM[i]*1000/trueStimDurMS, fmt='o',ecolor='black',color='black')
    ax_row1.plot(speed,sponList*1000/sponWindowMS, '--', label='spon')
    ax_row1.errorbar(speed,sponList*1000/sponWindowMS,
            yerr = sponSEM*1000/sponWindowMS, fmt='o',ecolor='black',color='black')
    
    ax_row1.set_title('Speed Tuning Plot', fontsize=8)
    ax_row1.set_xlabel('Speed (˚/sec)', fontsize = 8)
    ax_row1.set_ylim(bottom=0)
    ax_row1.set_xticks(speed)
    ax_row1.set_xticklabels(np.around(speed), fontsize=3)
    ax_row1.set_ylabel('Firing Rate (spikes/sec)', fontsize = 8)
    ax_row1.legend(prop={'size': 6})


    #### Hists
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
            for k,dir in enumerate(dirs):
                spikeHist = spikeHists[k,plotCount,:] * 1000/stimCount[0][plotCount]
                gaussSmooth = gaussian_filter1d(spikeHist, 5)
                if max(gaussSmooth) > yMax:
                    yMax = max(gaussSmooth)
                ax_row2[i,j].plot(gaussSmooth)
            ax_row2[i,j].set_title(f"{speed[plotCount]} ˚/sec", fontsize=7)
            plotCount += 1
            # ax_row2[i,j].set_ylim(bottom=0)
            # ax_row2[i,j].set_yticks([0,50,100])
            # ax_row2[i,j].set_yticklabels([0,50,100], fontsize=5)
            ax_row2[i,j].set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
            ax_row2[i,j].set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
            ax_row2[i,j].axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
            ax_row2[i,j].axhline(y=sponSpikesMean*1000/sponWindowMS, linestyle='--', color='grey')
            if i == 1 and j == 0:
                ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i,j].yaxis.set_label_coords(-0.3,0.3)
    plt.tight_layout(pad=0.5, w_pad=0.2, h_pad=0.2)

    for i in range(2):
        for j in range(3):
            ax_row2[i,j].set_ylim([0, yMax*1.1])
    # saves plot as pdf
    plt.savefig(f'{unit}.pdf')

# np.save('unitsSpeedTuningMat', allTuningMat)
plt.close('all')

'''
log normal fit
'''
