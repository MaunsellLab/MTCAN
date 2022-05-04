'''
to do:
incorp frame render to align spikes
'''

from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import time
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


# load data
allTrials, header = loadMatFile73('Testing', 'testing_220317', 'testing_220317_Dir_GRF_Spikes.mat')


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


# Tuning code
units = activeUnits('spikeData', allTrials)
correctTrials = correctTrialsGRF(allTrials)

numDir = np.int32(header['map0Settings']['data']['directionDeg']['n'])
stimDurMS = np.int32(header['mapStimDurationMS']['data'])
histPrePostMS = 50
allTuningMat = np.zeros((len(units),numDir))

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
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)
                        unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnSNEV) & 
                                        (unitTimeStamps <= stimOnSNEV + stimDurMS/1000))
                        spikeCountMat[stCount][0][dirIndex] = len(stimSpikes[0])
                        stimCount[0][dirIndex] += 1
                        
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
            histSmooth = smooth(spikeHist,75)
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
    
    # saves plot as png 
    os.makedirs('Direction Tuning')
    os.chdir('Direction Tuning/')
    plt.savefig(f'{unit}.png')

np.save('unitsDirTuningMat', allTuningMat)




#gauss fit

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gauss_fit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt

%matplotlib notebook

tc = np.array([10,23,45,80,37,16,12])
x = np.array([0,60,120,180,240,300,360])
x_full = np.linspace(0, 360, 1000)
params = gauss_fit(x, tcNorm)
y_full_fit = gauss(x_full, *params)

plt.plot(x_full, y_full_fit, '--r', label='fit')
plt.scatter(x, tcNorm, label= 'not fit')
plt.legend()

