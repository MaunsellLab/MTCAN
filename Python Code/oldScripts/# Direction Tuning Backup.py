# Direction Tuning Backup

'''
to do:
incorp frame render to align spikes
'''

from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import time
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


### for testing purposes, to make unit field similar to real data
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
allTrials, header = loadMatFile73('Meetz', '220607', 'Meetz_220607_GRF2_Spikes.mat')


# Tuning code
units = activeUnits('spikeData', allTrials)
correctTrials = correctTrialsGRF(allTrials)

#made edit here
frameRateHz = header['displayCalibration']['data']['frameRateHz'].tolist()
numDir = np.int32(header['map0Settings']['data']['directionDeg']['n'])
stimDurMS = np.int32(header['mapStimDurationMS']['data'])
histPrePostMS = 100
allTuningMat = np.zeros((len(units),numDir))

if not os.path.exists('Direction Tuning PDFs'):
    os.makedirs('Direction Tuning PDFs')
os.chdir('Direction Tuning PDFs/')


for uCount, unit in enumerate(units):
    all_data = []
    all_unit_ids = []
    all_stim_dirs = []

    stimCount = np.zeros((1,numDir))
    spikeCountMat = np.zeros((50,1, numDir))
    spikeCountMat[:,:,:] = np.nan
    #made edit here
    spikeHists = np.zeros((numDir, stimDurMS + 2*histPrePostMS+12))

    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        if 'numMap0Stim' in currTrial:
            map0StimLim = int(currTrial['numMap0Stim']['data'].tolist())
            map0Count = 0
            stimDesc = currTrial['stimDesc']['data']
            spikeData = currTrial['spikeData']
            #made edit here
            stim1TimeS = currTrial['taskEvents']['stimulusOn'][0]['time'].tolist()
            for stim in stimDesc:
                if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                    dirIndex = int(stim['directionIndex'])
                    stCount = int(stimCount[0][dirIndex])
                    #made edit here
                    stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'].tolist())
                                   /1000) + stim1TimeS
                    
                    stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'].tolist())
                                   /1000) + stim1TimeS
                    stimCount[0][dirIndex] += 1
                    map0Count += 1
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)[0]
                        #print(unitIndex)
                        if len(unitIndex) == 1: # maybe double check on this
                            unitTimeStamps = spikeData['timeStamp']
                            #continue
                        else:
                            unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        #made edit here
                        stimSpikes = np.where((unitTimeStamps >= stimOnTimeS) & 
                                        (unitTimeStamps <= stimOffTimeS))
                        spikeCountMat[stCount, 0, dirIndex] = len(stimSpikes[0])
                        
                        #check my histograms
                        #histograms
                        stimOnPreSNEV = stimOnTimeS - histPrePostMS/1000
                        stimOffPostSNEV = stimOffTimeS + histPrePostMS/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                            & (unitTimeStamps < stimOffPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)

                        histStimSpikes_binarized = np.zeros(spikeHists.shape[-1])
                        histStimSpikes_binarized[histStimSpikes] = 1
                        spikeHists[dirIndex, histStimSpikes] += 1
                        all_data.append(histStimSpikes_binarized)
                        #all_unit_ids.append(unit)
                        all_stim_dirs.append(dirIndex)


    all_data = np.array(all_data)
    #all_unit_ids = np.array(all_unit_ids)
    all_stim_dirs = np.array(all_stim_dirs)

    date = header['date']
    fig = plt.figure()
    fig.set_size_inches(6,8)

    text = fig.text(0.05, 0.85, f'Direction tuning for unit {unit}\n{date}\n- - - - -\n \
    Stimulus Duration = {stimDurMS} ms\nNumber of Blocks = {int(stimCount[0][0])}',\
                    size=10, fontweight='bold')
    text.set_path_effects([path_effects.Normal()])

    '''
    ax_row1 = plt.subplot2grid((10,6), (0,3), colspan = 3, rowspan = 4, polar=True)
    theta = np.radians(np.arange(0,420,360/numDir))
    r = (np.append(spikeCountMean, spikeCountMean[0][0]))*1000/stimDurMS
    err = (np.append(spikeCountSD, spikeCountSD[0][0]))*1000/stimDurMS
    ax_row1.plot(theta,r)
    ax_row1.errorbar(theta, r, yerr = err,fmt='o', ecolor = 'black', color='black')
    ax_row1.set_theta_zero_location("N")
    ax_row1.set_rmax(100)
    ax_row1.set_title('Direction Tuning Plot', fontsize=8)
    '''

    # hists

    #made edit here
    # spikeHistsRS = np.reshape(spikeHists, (stimDurMS + 162,2,3))
    stimCountRS = np.reshape(stimCount, (2,3))
    titleArr = np.reshape(np.arange(0,360,60),(2,3))

    ax_row2 = []
    for countI, i in enumerate(range(4, 10, 3)):
        ax = []
        for countJ, j in enumerate(range(0, 6, 2)):
            ax.append(plt.subplot2grid((10,6), (i,j), colspan = 2, rowspan = 3))
        ax_row2.append(np.array(ax))
    ax_row2 = np.array(ax_row2) # 2 x 3


    plotCount = 0
    for i in range(2):
        for j in range(3):
            spikeHist = spikeHists[plotCount,:] * 1000/stimCountRS[i,j] 
            #rel_rows = np.where(all_stim_dirs == plotCount)[0]

            #spikeHist = all_data[rel_rows,:].mean(0)

            plotCount += 1
            histSmooth = smooth(spikeHist,50)#*1000
            ax_row2[i,j].plot(histSmooth)
            histTitle = titleArr[i][j]
            ax_row2[i,j].set_title(f"{histTitle}˚", fontsize=7)
            ax_row2[i,j].set_ylim([0, 100])
            ax_row2[i,j].set_yticks([0,50,100])
            ax_row2[i,j].set_yticklabels([0,50,100], fontsize=5)
            # ax_row2[i,j].set_xticks([histPrePostMS,histPrePostMS+stimDurMS])
            ax_row2[i,j].set_xticks([0,50,100,150,200,250,300,350,400,450])
            ax_row2[i,j].set_xticklabels([0,50,100,150,200,250,300,350,400,450], fontsize=5)
            if i == 1 and j == 0:
                ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i,j].yaxis.set_label_coords(-0.2,0.3)
    plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)
    
    # saves plot as pdf 
    plt.savefig(f'{unit}.pdf')
    continue
plt.close('all')


    ######################
    spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)
    spikeCountSD = ma.std(ma.masked_invalid(spikeCountMat), axis = 0)
    allTuningMat[uCount] = spikeCountMean

    #plot figure
    #polar
    date = header['date']

    fig = plt.figure()
    fig.set_size_inches(6,8)

    text = fig.text(0.05, 0.85, f'Direction tuning for unit {unit}\n{date}\n- - - - -\n \
    Stimulus Duration = {stimDurMS} ms\nNumber of Blocks = {int(stimCount[0][0])}',\
                    size=10, fontweight='bold')
    text.set_path_effects([path_effects.Normal()])


    ax_row1 = plt.subplot2grid((10,6), (0,3), colspan = 3, rowspan = 4, polar=True)
    theta = np.radians(np.arange(0,420,360/numDir))
    r = (np.append(spikeCountMean, spikeCountMean[0][0]))*1000/stimDurMS
    err = (np.append(spikeCountSD, spikeCountSD[0][0]))*1000/stimDurMS
    ax_row1.plot(theta,r)
    ax_row1.errorbar(theta, r, yerr = err,fmt='o', ecolor = 'black', color='black')
    ax_row1.set_theta_zero_location("N")
    ax_row1.set_rmax(100)
    ax_row1.set_title('Direction Tuning Plot', fontsize=8)

    # hists

    #made edit here
    # spikeHistsRS = np.reshape(spikeHists, (stimDurMS + 162,2,3))
    stimCountRS = np.reshape(stimCount, (2,3))
    titleArr = np.reshape(np.arange(0,360,60),(2,3))

    ax_row2 = []
    for countI, i in enumerate(range(4, 10, 3)):
        ax = []
        for countJ, j in enumerate(range(0, 6, 2)):
            ax.append(plt.subplot2grid((10,6), (i,j), colspan = 2, rowspan = 3))
        ax_row2.append(np.array(ax))
    ax_row2 = np.array(ax_row2) # 2 x 3

    plotCount = 0
    for i in range(2):
        for j in range(3):
            spikeHist = spikeHists[plotCount,:] * 1000/stimCountRS[i,j] 
            plotCount += 1
            histSmooth = smooth(spikeHist,50)
            ax_row2[i,j].plot(histSmooth)
            histTitle = titleArr[i][j]
            ax_row2[i,j].set_title(f"{histTitle}˚", fontsize=7)
            ax_row2[i,j].set_ylim([0, 100])
            ax_row2[i,j].set_yticks([0,50,100])
            ax_row2[i,j].set_yticklabels([0,50,100], fontsize=5)
            ax_row2[i,j].set_xticks([75,75+stimDurMS])
            ax_row2[i,j].set_xticklabels([75,75+stimDurMS], fontsize=5)
            if i == 1 and j == 0:
                ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i,j].yaxis.set_label_coords(-0.2,0.3)
    plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)
    
    # saves plot as pdf 
    plt.savefig(f'{unit}.pdf')

plt.close('all')
# np.save('unitsDirTuningMat', allTuningMat)

