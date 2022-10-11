'''
MTN2 Analysis Script


CHECK LOW/HIGH Contrast conditions 

to do:
convert heatmap to spikes/sec, it's at spikes/stimDurMS
'''

import seaborn as sns
import numpy as np
import numpy.ma as ma
from usefulFns import *
from itertools import combinations
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time

'''
### testing
allTrials, header = loadMatFile73('Testing', 'Meetz_220621', 'Meetz_220621_2_MTN.mat')
###
'''

####### START HERE ######

# load relevant file 
allTrials, header = loadMatFile73('Meetz', '221010', 'Meetz_221010_MTNC_Spikes.mat')


if not os.path.exists('Normalization'):
    os.makedirs('Normalization')
os.chdir('Normalization/')


## generate list of unique active units
units = activeUnits('spikeData', allTrials)


## list of indices of correctTrials (non-instruct, valid trialCertify)
corrTrials = correctTrialsMTX(allTrials)

## assert: stim sequence list is frozen
seqList = []
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    stimDesc = currTrial['stimDesc']['data']
    for stim in stimDesc:
        if stim['stimLoc'] == 0 and stim['listType'] == 1:
            if len(seqList) < 49:
                seqList.append(int(stim['stimIndex'].tolist()))
                seqArr = np.array(seqList)
                lastIndex = int(stim['stimIndex'].tolist())
            else:
                posLastIndex = np.where(seqArr==lastIndex)[0][0]
                if posLastIndex == len(seqArr)-1:
                    if int(stim['stimIndex'].tolist()) != seqArr[0]:
                        print('out of sequence')
                    else:
                        lastIndex = int(stim['stimIndex'].tolist())
                else:
                    if int(stim['stimIndex'].tolist()) != seqArr[posLastIndex+1]:
                        print('out of sequence')
                    else:
                        lastIndex = int(stim['stimIndex'].tolist())


## assert: are there correct trials without spikeData
noSpikeData = []
for trialCount, currTrial in enumerate(allTrials):
    trial = currTrial['trial']['data']
    extendedEOT = currTrial['extendedEOT']['data']
    if extendedEOT == 0 and trial['instructTrial'] != 1:
        if 'spikeData' not in currTrial:
            noSpikeData.append(trialCount)


## assert: frame consistency during stimlus duration
frameRateHz = header['frameRateHz']['data'].tolist()
stimDurFrame = []
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    stimDesc = currTrial['stimDesc']['data']
    for stim in stimDesc:
        if stim['stimLoc'] == 0:
            frameDiff = stim['stimOffFrame'].tolist() - stim['stimOnFrame'].tolist()
            stimDurFrame.append(frameDiff)
if len(set(stimDurFrame)) != 1:
    print('stimulus frame duration not consistent for mapping stimuli')
else: 
    trueStimDurMS = np.int32(np.around(1000/frameRateHz * stimDurFrame[0]))


# generates a dictionary, numpy array, and Pandas Dataframe of stim Index 
# and corresponding directions/contrasts
stimIndexDict = {}
for currTrial in allTrials:
    extendedEOT = currTrial['extendedEOT']['data']
    trial = currTrial['trial']['data']
    if extendedEOT == 0 and trial['instructTrial'] != 1:
        stimDesc = currTrial['stimDesc']['data']
        for stim in stimDesc:
            if stim['stimLoc'] != 2:
                index = int(stim['stimIndex'].tolist())
                if index not in stimIndexDict:
                    stimIndexDict[index] = {}
                    if int(stim['stimLoc'].tolist()) not in stimIndexDict[index]:
                        stimIndexDict[index][int(stim['stimLoc'].tolist())] = \
                        {'direction': int(stim['directionDeg'].tolist()),
                         'contrast': np.around(stim['contrast'],2).tolist()}
                else:
                    if int(stim['stimLoc'].tolist()) not in stimIndexDict[index]:
                        stimIndexDict[index][int(stim['stimLoc'].tolist())] = \
                        {'direction': int(stim['directionDeg'].tolist()),
                         'contrast': np.around(stim['contrast'],2).tolist()}
stimIndexArray = np.zeros((49,4))
for i in range(len(stimIndexDict)):
    stimIndexArray[i][0] = stimIndexDict[i][0]['direction']
    stimIndexArray[i][1] = stimIndexDict[i][0]['contrast']
    stimIndexArray[i][2] = stimIndexDict[i][1]['direction']
    stimIndexArray[i][3] = stimIndexDict[i][1]['contrast']
stimIndexDF = pd.DataFrame(stimIndexArray, columns=['loc0 Direction', 'loc0 Contrast',
                                            'loc1 Direction', 'loc1 Contrast'])


#initialize lists/arrays/dataframes for counting spikeCounts and for analysis
blocksDone = int(allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone']
                 .tolist()) 
highContrast, zeroContrast = max(stimIndexDF['loc0 Contrast'].unique()), \
                             min(stimIndexDF['loc0 Contrast'].unique())
zeroDir = 0
dirArray = np.array([0,60,120,180,240,300])
spikeCountMat = np.zeros((len(units),blocksDone+1,49))
spikeCountLong = []
sponSpikeCountLong = []
histPrePostMS = 100 #100ms window pre/post stimlus on/off
sponWindowMS = 50 #50ms window before stimulus onset
spikeHists = np.zeros((len(units),49, trueStimDurMS+(2*histPrePostMS+1)))
stimIndexCount = np.zeros(49) 


# insert spike counts into matrix of unique stimulus sets
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    if 'spikeData' in currTrial:
        stimDesc = currTrial['stimDesc']['data']
        stim1TimeS = currTrial['taskEvents']['stimulusOn'][0]['time'].tolist()
        fixateTimeS = currTrial['taskEvents']['fixate']['time'].tolist()
        for stim in stimDesc:
            if stim['stimLoc'] == 0 and stim['listType'] == 1:
                stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'].tolist())
                                /1000) + stim1TimeS
                stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'].tolist())
                                /1000) + stim1TimeS
                stimIndex = np.int32(stim['stimIndex'])
                stCount = int(stimIndexCount[stimIndex])
                stimIndexCount[stimIndex] += 1
                for unitCount, unit in enumerate(units):
                    if unit in currTrial['spikeData']['unit']:
                        unitIndex = np.where(currTrial['spikeData']['unit'] == unit)[0]
                        #added 50ms onset latency for spike counts (100 for offset)
                        unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= (stimOnTimeS+0.05)) & 
                                    (unitTimeStamps <= (stimOffTimeS+0.05)))[0]
                        spikeCountMat[unitCount][stCount][stimIndex] \
                        = len(stimSpikes)
                        spikeCountLong.append([unit, stimIndex, stimIndexCount[stimIndex], len(stimSpikes)])

                        #Spontaneous Spikes
                        sponSpikes = np.where((unitTimeStamps >= (stimOnTimeS-(sponWindowMS/1000))) 
                                            & (unitTimeStamps <= stimOnTimeS))[0]
                        sponSpikeCountLong.append([unit,len(sponSpikes)])

                        #PSTHs
                        stimOnPreSNEV = stimOnTimeS - (histPrePostMS/1000)
                        stimOffPostSNEV = stimOffTimeS + (histPrePostMS/1000)
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)\
                                    & (unitTimeStamps <= stimOffPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        # spikeHists[unitCount, stimIndex, histStimSpikes] += 1
                        spikeHists[unitCount][stimIndex][histStimSpikes] += 1


# mean, SEM, and reshaping of spikeCount matrices 
# create pandas dataframe of spikeCount with corresponding unit, stimIndex
spikeCountDF = pd.DataFrame(spikeCountLong, columns=['unit', 'stimIndex', 
                                                     'stimCount', 'stimSpikes'])
sponSpikeCountDF = pd.DataFrame(sponSpikeCountLong, columns =['unit', 'sponSpikes'])

sponSpikesMean = np.zeros(len(units)) # spikes in 50ms window
sponSpikesSEM = np.zeros(len(units)) 
for unitCount, unit in enumerate(units):
    sponSpikesMean[unitCount] = sponSpikeCountDF.loc[sponSpikeCountDF['unit']==unit].mean()[1]
    sponSpikesSEM[unitCount] = sponSpikeCountDF.loc[sponSpikeCountDF['unit']==unit].mean()[1]
meanSpike = np.mean(spikeCountMat[:,:blocksDone,:], axis = 1)
spikeCountSD = np.std(spikeCountMat[:,:blocksDone,:], axis = 1)
spikeCountSEM = spikeCountSD/np.sqrt(blocksDone)
meanSpikeReshaped = np.zeros((len(units),1,49))
for count,i in enumerate(meanSpikeReshaped):
    #row1 of 7x7 grid when reshaped
    i[:,0:6] = meanSpike[count][0:6]
    i[:,6] = meanSpike[count][42]
    #row2
    i[:,7:13] = meanSpike[count][6:12]
    i[:,13] = meanSpike[count][43]
    #row3
    i[:,14:20] = meanSpike[count][12:18]
    i[:,20] = meanSpike[count][44]
    #row4
    i[:,21:27] = meanSpike[count][18:24]
    i[:,27] = meanSpike[count][45]
    #row5
    i[:,28:34] = meanSpike[count][24:30]
    i[:,34] = meanSpike[count][46]
    #row6
    i[:,35:41] = meanSpike[count][30:36]
    i[:,41] = meanSpike[count][47]
    #row7
    i[:,42:48] = meanSpike[count][36:42]
    i[:,48] = meanSpike[count][48]


## heatmap of normalization
for unit in range(len(units)):

    b = meanSpikeReshaped[unit].reshape(7,7)
    bSmooth = gaussian_filter(b, sigma=1) * 1000/trueStimDurMS
    prefDir, nullDir = unitPrefNullDir(bSmooth)

    #using seaborn
    ax = sns.heatmap(bSmooth, square=True, linewidths=0.2, vmin=0)
    ax.set_xticks(np.arange(7)+0.5)
    ax.set_title(f'heatmap of normalization for {units[unit]}')
    ax.set_xticklabels(['0','60','120','180','240','300','blank'], rotation = 45)
    ax.xaxis.set_ticks_position("top")

    ax.set_yticks(np.arange(7)+0.5)
    ax.set_yticklabels(['0','60','120','180','240','300','blank'], rotation = 0)
    plt.tight_layout()
    plt.savefig(f'{units[unit]}NormHeatmapNonGaussianFilter.pdf')
    plt.close('all')

for unit in range(len(units)):

    b = meanSpikeReshaped[unit].reshape(7,7)
    bSmooth = gaussian_filter(b, sigma=1) * 1000/trueStimDurMS
    prefDir, nullDir = unitPrefNullDir(meanSpikeReshaped)

    # a = meanSpikeReshaped[unit]
    # b = a.reshape(13,13)
    # bSmooth = gaussian_filter(b, sigma=1)

    # maxLoc0 = max(bSmooth[12,:])
    # maxLoc1 = max(bSmooth[:,12])
    # if maxLoc0 > maxLoc1:
    #     prefDir = dirArray[np.where(bSmooth[12,:]==maxLoc0)[0][0]]
    # else:
    #     prefDir = dirArray[np.where(bSmooth[:,12]==maxLoc1)[0][0]]
    # nullDir = (prefDir + 180)%360
    nullIndex = np.where(dirArray==nullDir)[0][0]
    reIndex = (np.array([0,1,2,3,4,5])+nullIndex) % 6

    bSmoothReIndex = bSmooth[:6,:6]
    bSmoothReIndex = bSmoothReIndex[:,reIndex]
    bSmoothReIndex = bSmoothReIndex[reIndex,:]
    b0Blank = bSmooth[:6,6]
    b0Blank = b0Blank[reIndex]
    b1Blank = bSmooth[6,:6]
    b1Blank = b1Blank[reIndex]
    
    gridspec_kw = {"height_ratios":[6,1], "width_ratios" : [6,1]}
    heatmapkws = dict(square=False, cbar=False, linewidths=1.0, vmin=0, vmax=np.max(bSmooth))
    asp = bSmooth.shape[0]/float(bSmooth.shape[1])
    figW = 8
    figH = figW*asp
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(figW, figH), gridspec_kw=gridspec_kw)
    # left = 0.07; right=0.7
    # bottom = 0.1; top = 0.1
    # plt.subplots_adjust(left=left,right=right,bottom=bottom,top=top,wspace=0.1,hspace=0.1*asp)
    plt.subplots_adjust(wspace=0.1,hspace=0.1*asp)
    sns.heatmap(bSmoothReIndex, ax=axes[0,0], xticklabels=True, yticklabels=True, **heatmapkws)
    sns.heatmap(b0Blank, ax=axes[0,1], xticklabels=True, yticklabels=False, **heatmapkws)
    sns.heatmap(b1Blank, ax=axes[1,0], xticklabels=False, yticklabels=True, **heatmapkws)
    sns.heatmap(bSmooth[6,6], ax=axes[1,1], xticklabels=False, yticklabels=False, **heatmapkws)

    axes[0,0].xaxis.set_ticks_position("top")
    axes[0,1].xaxis.set_ticks_position("top")
    axes[0,0].set_xticklabels(dirArray[:6][reIndex[:6]], rotation = 45)
    axes[0,1].set_xticklabels(dirArray[:6][reIndex[:6]], rotation = 45)
    axes[0,0].set_yticklabels(dirArray[:6][reIndex[:6]], rotation = 0)
    axes[1,0].set_yticklabels(dirArray[:6][reIndex[:6]], rotation = 0)
    axes[0,0].set_xlabel('loc 0 contrast')
    axes[0,0].xaxis.set_label_position('top') 
    axes[0,0].set_ylabel('loc 1 contrast')
    plt.tight_layout()
    plt.show()


## scipy curveFit Normalization parameters
def func(fixed, L_0, L_60, L_120, L_180, L_240, L_300, sL0, sL1, aL, sig):
    c0,c1,l0,l1 = fixed
    L = np.array([L_0, L_60, L_120, L_180, L_240, L_300])
    # return((c0*sL0*L*l0).sum(-1) + (c1*sL1*L*l1).sum(-1))/((aL0*c0[:,0])+(aL1*c1[:,0])+sig)
    # return (((c0*sL0*(l0*L)).sum(-1) + (c1*sL1*(l1*L)).sum(-1))/((aL0*c0[:,0])+(aL1*c1[:,0])+sig))
    return (((c0*sL0*(l0*L)).sum(-1) + (c1*sL1*(l1*L)).sum(-1))/((aL*c0[:,0])+(aL*c1[:,0])+sig))
for unitCount, unit in enumerate(units):
    resp = np.reshape(spikeCountMat[unitCount][:blocksDone+1,:],(49*blocksDone))
    fixParam = np.tile(np.arange(49), blocksDone)

    c0s, c1s, l0s, l1s = [], [], [], []
    direction_set = np.arange(0, 360, 60)
    for i in fixParam:
        c0 = stimIndexDict[i][0]['contrast']
        l0 = stimIndexDict[i][0]['direction']
        c1 = stimIndexDict[i][1]['contrast']
        l1 = stimIndexDict[i][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(6)
        l1_oh = np.zeros(6)
        l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
        l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
        c0s.append(np.repeat(c0,6))
        c1s.append(np.repeat(c1,6))
        l0s.append(l0_oh)
        l1s.append(l1_oh)

    c0s = np.array(c0s)
    c1s = np.array(c1s)
    l0s = np.array(l0s)
    l1s = np.array(l1s)

    pOpt, pCov = curve_fit(func, (c0s, c1s, l0s, l1s), resp, bounds=(
        (0,0,0,0,0,0,0,0,0,0),(np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,1,1,1,1)))
    print(unit,pOpt)


## PSTHs for P+P, N+N, P+N, and converse for other location
for unitCount, unit in enumerate(units):
    print(unit)
    b = meanSpikeReshaped[unitCount].reshape(7,7)
    bSmooth = gaussian_filter(b, sigma=1) * 1000/trueStimDurMS
    prefDir, nullDir = unitPrefNullDir(b)

    # figure 
    fig = plt.figure()
    fig.set_size_inches(4,2)
    ax = []
    for col in range(2):
        ax.append(plt.subplot2grid((1,2), (0,col)))
    ax = np.array(ax)

    locDirList = [(prefDir,prefDir),(nullDir,nullDir),(prefDir,nullDir),
                  (prefDir,prefDir),(nullDir,nullDir),(nullDir,prefDir)]
    yMax = 0
    plotCount = 0
    subCount = 0
    for locDir in locDirList:
        loc0Dir, loc1Dir = locDir
        histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) & 
                                      (stimIndexDF['loc0 Contrast'] == highContrast) & 
                                      (stimIndexDF['loc1 Direction'] == loc1Dir) & 
                                      (stimIndexDF['loc1 Contrast'] == highContrast)][0]
        dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
        smoothPlot = gaussian_filter1d(dirPlot,5)
        if max(smoothPlot) > yMax:
            yMax = max(smoothPlot)
        ax[plotCount].plot(smoothPlot, label=f'{loc0Dir}+{loc1Dir}')
        ax[plotCount].set_title(f'loc0 direction {loc0Dir} loc1 direction {loc1Dir}', fontsize= 5)
        ax[plotCount].set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
        ax[plotCount].set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
        ax[plotCount].set_xlim([0,trueStimDurMS+(2*histPrePostMS+1)])
        ax[plotCount].axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.1)
        ax[plotCount].axhline(y=sponSpikesMean[unitCount]*1000/sponWindowMS, linestyle='--', color='grey')
        ax[plotCount].set_xlabel('Stimulus Duration (ms)', fontsize=7)
        if plotCount == 0 or plotCount == 1:
            ax[plotCount].legend(loc='upper left', prop={'size': 4})
        if plotCount == 0:
            ax[plotCount].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
        subCount += 1
        if subCount % 3 == 0:
            plotCount += 1
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    for i in range(2):
        ax[i].set_ylim([0,yMax*1.1])
    plt.savefig(f'{unit}PrefNull.pdf')
    plt.close('all')


## mean response of PP, PN, and NN to low/high contrast
# create pd.dataframe for plotting
pnContrastPlotDF = pd.DataFrame(columns=['unit', 'stimIndex', 'stimCount', 
                                'stimSpikes', 'contrast', 'prefNullStr'])
for unitCount, unit in enumerate(units):
    b = meanSpikeReshaped[unitCount].reshape(7,7)
    bSmooth = gaussian_filter(b, sigma=1) * 1000/trueStimDurMS
    prefDir, nullDir = unitPrefNullDir(bSmooth)
    print(unit, prefDir, nullDir)
    orientCount = 0
    # orientList = ['pref+pref', 'pref+null', 'null+pref', 'null+null', 
    #               'prefDir+blank', 'nullDir+blank','blank+prefDir', 'blank+nullDir']
    orientList = ['pref', 'pref+null', 'null+pref', 'null', 
                  'pref', 'null','pref', 'null']
    for j in [(prefDir,prefDir),(prefDir,nullDir),(nullDir,prefDir),(nullDir,nullDir),]:
        loc0Dir, loc1Dir = j
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == highContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == highContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit) 
                    & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = highContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1

    #Loc 0 Pref/Null only
    for i in [(prefDir, zeroDir), (nullDir, zeroDir)]:
        loc0Dir, loc1Dir = i  
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == highContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == zeroContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit) 
                    & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = zeroContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1

    #loc 1 Pref/Null only
    for x in [(zeroDir,prefDir), (zeroDir,nullDir)]:
        loc0Dir, loc1Dir = x  
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == zeroContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == highContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit) 
                    & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = zeroContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1


#figure
for unit in units:
    unitDF = pnContrastPlotDF.loc[pnContrastPlotDF['unit'] == unit]
    unitDF['stimSpikes'] = unitDF['stimSpikes'] * 1000/trueStimDurMS #spikeCounts in spikes/sec
    ax = sns.catplot(data=unitDF, x='contrast', y='stimSpikes', hue='prefNullStr', 
                kind='point', errorbar="se")
    ax.set(ylim=(0, None))
    ax.refline(y=(sponSpikesMean[unitCount]*1000/sponWindowMS), color = 'grey')
    ax.set(xlabel='Contrast', ylabel='Firing Rate (spikes/sec)')
    plt.savefig(f'{unit}BarPlot.pdf')
    plt.close('all')


## PSTHs for P+P, N+N, P+I, and converse for other location
for unitCount, unit in enumerate(units):
    if unit != 61: ####### REMOVE THIS
        print(unit)
        b = meanSpikeReshaped[unitCount].reshape(13,13)
        bSmooth = gaussian_filter(b, sigma=1)
        maxLoc0 = max(bSmooth[12,:])
        maxLoc1 = max(bSmooth[:,12])
        if maxLoc0 > maxLoc1:
            prefDir = dirArray[np.where(bSmooth[12,:]==maxLoc0)[0][0]]
        else:
            prefDir = dirArray[np.where(bSmooth[:,12]==maxLoc1)[0][0]]
        nullDir = (prefDir + 180)%360
        otherDir = [i for i in dirArray[:6] if i != prefDir and i != nullDir]
        contrastList = [[lowContrast, lowContrast],[lowContrast,highContrast],
                        [highContrast,lowContrast],[highContrast,highContrast]]
        locDirList = [(prefDir,prefDir),(nullDir,nullDir),(prefDir,nullDir),
            (prefDir,otherDir[0]),(prefDir,otherDir[1]),(prefDir,otherDir[2]),
            (prefDir,otherDir[3]),(prefDir,prefDir),(nullDir,nullDir),(nullDir,prefDir),
            (otherDir[0],prefDir),(otherDir[1],prefDir),(otherDir[2],prefDir),
            (otherDir[3],prefDir)]
        
        fig = plt.figure()
        fig.set_size_inches(8,4)
        ax = []
        for row in range(2):
            for col in range(4):
                ax.append(plt.subplot2grid((2,4), (row,col)))
        ax = np.array(ax)

        yMax = 0
        plotCount = 0
        for x in contrastList:
            subCount = 0
            loc0Con, loc1Con = x
            for locDir in locDirList:
                loc0Dir, loc1Dir = locDir
                histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) & 
                (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == loc1Dir)
                & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
                dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
                smoothPlot = gaussian_filter1d(dirPlot,5)
                if max(smoothPlot) > yMax:
                    yMax = max(smoothPlot)
                ax[plotCount].plot(smoothPlot, label=f'{loc0Dir}+{loc1Dir}')
                ax[plotCount].set_title(f'loc0 contrast {loc0Con} loc1 contrast {loc1Con}', fontsize= 5)
                ax[plotCount].set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
                ax[plotCount].set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
                ax[plotCount].set_xlim([0,trueStimDurMS+(2*histPrePostMS+1)])
                ax[plotCount].axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.1)
                if plotCount == 0 or plotCount == 1:
                    ax[plotCount].legend(loc='upper left', prop={'size': 4})
                subCount += 1
                if subCount % 7 == 0:
                    plotCount += 1
        plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
        for i in range(8):
            ax[i].set_ylim([0,yMax*1.1])
        plt.savefig(f'{unit}.pdf')
        plt.close('all')


## correlations incomplete 
#filtered units test
combs = [i for i in combinations(filterUnits, 2)]
#normal units 
combs = [i for i in combinations(units, 2)]
corrMat = np.zeros(len(combs))

# z-scored spikeCountMat
zSpikeCountMat = stats.zscore(spikeCountMat[:,:blocksDone,:], axis=1, nan_policy='omit') 
zSpikeCountMat = np.nan_to_num(zSpikeCountMat)
zSpikeCountMat = np.reshape(zSpikeCountMat,(len(units),blocksDone*49))
for count, i in enumerate(combs):
    n1 = np.where(units == i[0])[0][0]
    n2 = np.where(units == i[1])[0][0]
    pairCorr = stats.pearsonr(zSpikeCountMat[n1],zSpikeCountMat[n2])
    corrMat[count] = pairCorr[0]
    # for j in range(np.shape(spikeCountMat)[2]):
    #     stimCorr = stats.pearsonr(zSpikeCountMat[n1,1:blocksDone+1,j],
    #                              zSpikeCountMat[n2,1:blocksDone+1,j])
    #     # stimCorr = stats.pearsonr(spikeCountMat[n1,:blocksDone+1,j],
    #     #                          spikeCountMat[n2,:blocksDone+1,j])

popCorr = np.mean(corrMat)


## RF location tuning similarity b/w neurons Bhattacharyya Distance 2D
RFLocMat = np.load('../RFLoc Tuning/unitsRFLocMat.npy')
combs = [i for i in combinations(units, 2)]
pairLocSimScore = np.zeros((len(combs),1))
for pairCount, pair in enumerate(combs):
    n1 = np.where(units == pair[0])[0][0]
    n2 = np.where(units == pair[1])[0][0]

    a = np.flip(RFLocMat[n1], axis=0)
    b = np.flip(RFLocMat[n2], axis=0)
    m1, cov1, p = gauss2dParams(a)
    m2, cov2, p2 = gauss2dParams(b)
    BC = bhattCoef2D(m1,m2,cov1,cov2)
    pairLocSimScore[pairCount] = BC


## direction tuning similarity b/w neurons Bhattacharyya Distance
dirTuningMat = np.load('../Direction Tuning/unitsDirTuningMat.npy') #load from directions folder
# unitsBaselineMat = np.load('../Direction Tuning/unitsBaselineMat.npy')
extTunMat = np.concatenate((dirTuningMat[:,3:], dirTuningMat[:,:], 
                            dirTuningMat[:,:3]), axis=1)
angleMat = np.arange(180,900,60)

combs = [i for i in combinations(units, 2)]
pairSimPrefDir = np.zeros((len(combs),1))
pairSimScore = np.zeros((len(combs),1))
for pairCount, pair in enumerate(combs):
    n1 = np.where(units == pair[0])[0][0]
    n2 = np.where(units == pair[1])[0][0]
    n1Max = int(np.where(dirTuningMat[n1] == np.max(dirTuningMat[n1]))[0] + 3)
    n1X = angleMat[n1Max-3:n1Max+4]
    # n1Y = extTunMat[n1][n1Max-3:n1Max+4]
    n1Y = extTunMat[n1][n1Max-3:n1Max+4]/max(extTunMat[n1][n1Max-3:n1Max+4])
    n1XFull = np.linspace(n1X[0],n1X[-1],1000)
    params = gaussFit(n1X, n1Y)
    # n1YFull = gauss(n1XFull, *params)
    m1 = params[2] # mean neuron 1
    v1 = params[3]**2 # var neuron 1
    n1TrueMean = m1 - 360
    
    n2Max = int(np.where(dirTuningMat[n2] == np.max(dirTuningMat[n2]))[0] + 3)
    n2X = angleMat[n2Max-3:n2Max+4]
    # n2Y = extTunMat[n2][n2Max-3:n2Max+4]
    n2Y = extTunMat[n2][n2Max-3:n2Max+4]/max(extTunMat[n2][n2Max-3:n2Max+4])
    n2XFull = np.linspace(n2X[0], n2X[-1],1000)
    params = gaussFit(n2X, n2Y)
    # n2YFull = gauss(n2XFull, *params)
    m2 = params[2]
    v2 = params[3]**2
    n2TrueMean = m2 - 360

    if abs(m1-m2) > 180:
        if m1 > m2:
            m1 = m2-(360-(m1-m2))
        else:
            m2 = m1-(360-(m2-m1))

    # similarity of pref dirs only
    pairSimPrefDir[pairCount] = 1 - (abs(m1-m2)/180)
    # bhattacharyya similarity score 
    BC = bhattCoef(m1, m2, v1, v2)
    pairSimScore[pairCount] = BC


## unit Direction Selectivity 
unitsBaselineMat = np.load('../Direction Tuning/unitsBaselineMat.npy')
unitSelectivity = np.zeros(len(units))
for unit in filterUnits:
    unitIndex = np.where(units == unit)[0][0]
    nMax = np.where(dirTuningMat[unitIndex] == np.max(dirTuningMat[unitIndex]))[0].squeeze() 
    prefDir = dirArray[nMax]
    nullDir = (prefDir + 180) % 360
    nullResp = dirTuningMat[unitIndex][np.where(dirArray==nullDir)[0][0]]
    prefResp = dirTuningMat[unitIndex][nMax]
    baseResp = unitsBaselineMat[unitIndex]
    DS = 1 - ((nullResp-baseResp)/(prefResp-baseResp))
    unitSelectivity[unitIndex] = DS


#pair Selectivity
combs = [i for i in combinations(filterUnits, 2)]
pairSelScore = np.zeros((len(combs),1))
for pairCount, pair in enumerate(combs):
    n1 = np.where(units == pair[0])[0][0]
    n2 = np.where(units == pair[1])[0][0]
    geoMean = np.sqrt((unitSelectivity[n1]*unitSelectivity[n2]))
    pairSelScore[pairCount] = geoMean


## selectivity for direction (similar to Bram)
unitSelectivity = np.zeros((len(units),49)) 
unitSelectivity[:,:] = np.nan
for uCount, unit in enumerate(units):
    for stim in range(49):
        loc0Contrast = stimIndexDict[stim][0]['contrast']
        loc1Contrast = stimIndexDict[stim][1]['contrast']
        if loc0Contrast != 0 and loc0Contrast != 0:
            dir1 = stimIndexDict[stim][0]['direction']
            con1 = stimIndexDict[stim][0]['contrast']
            l1Index = np.where((stimIndexArray[:,0]== dir1) & (stimIndexArray[:,1] == con1)
                              & (stimIndexArray[:,3]==0))
            l1 = meanSpike[uCount][l1Index[0]]
            dir2 = stimIndexDict[stim][1]['direction']
            con2 = stimIndexDict[stim][1]['contrast']
            l2Index = np.where((stimIndexArray[:,2]== dir2) & (stimIndexArray[:,3] == con2)
                              & (stimIndexArray[:,1]==0))
            l2 = meanSpike[uCount][l2Index[0]]
            unitSelectivity[uCount][stim] = (l1-l2)/(l1+l2)


'''
'''
###### mean response of PP, PN, NN to low and high contrast
for unitCount, unit in enumerate(units):
    if unit != 61: ########## REMOVE THIS 
        prefDir, nullDir, bSmooth = unitPrefNullDir(meanSpikeReshaped, unitCount)
        print(unit, prefDir, nullDir)

        for count, i in enumerate([lowContrast, highContrast]):
            colorCount = 0
            if i == lowContrast:
                x = 1
            else:
                x = 2
            for j in [(prefDir, prefDir),(prefDir,nullDir),(nullDir,prefDir),(nullDir,nullDir)]:
                loc0Dir, loc1Dir = j
                bSmoothi = np.where(dirArray==loc0Dir)[0][count]
                bSmoothj = np.where(dirArray==loc1Dir)[0][count]
                # spikeIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) & 
                # (stimIndexDF['loc0 Contrast'] == i) & (stimIndexDF['loc1 Direction'] == loc1Dir)
                # & (stimIndexDF['loc1 Contrast'] == i)][0]
                # avgSpike = meanSpike[unitCount][spikeIndex] * 1000/trueStimDurMS
                # avgSpike = c[spikeIndex]
                avgSpike = bSmooth[bSmoothi][bSmoothj] * 1000/trueStimDurMS
                if colorCount == 0:
                    plt.scatter(x, avgSpike, c='blue', label='pref')
                if colorCount == 1 or colorCount == 2: 
                    plt.scatter(x,avgSpike, c='brown', label='pref+null/null+pref')
                if colorCount == 3:
                    plt.scatter(x,avgSpike, c='red', label='null')
                colorCount +=1 
            if count == 0:
                plt.legend(loc='upper right', prop={'size': 6})
        plt.ylabel('Firing Rate spikes/sec')
        plt.xlabel('Contrast')
        plt.title(f'Firing Rate vs Contrast for unit {unit}')
        plt.xticks([1,2], ['Low Contrast', 'High Contrast'])
        plt.xlim(left=0.5,right=2.5)
        plt.ylim(bottom=0)
        plt.savefig(f'{unit}ContrastResponse.pdf')
        plt.close('all')


# PSTHs for P,N, P+N
yMax = 0
unit2Pref = spikeHists[5,129,:] * 1000/stimIndexCount[129]
unit2Null = spikeHists[5,108,:] * 1000/stimIndexCount[108]
unit2PN = spikeHists[5,111,:] * 1000/stimIndexCount[111]
gaussSmoothPref = gaussian_filter1d(unit2Pref, 10)
gaussSmoothNull = gaussian_filter1d(unit2Null, 10)
gaussSmoothPN = gaussian_filter1d(unit2PN, 10)
if max(gaussSmoothPN) > yMax:
    yMax = max(gaussSmoothPN)
if max(gaussSmoothNull) > yMax:
    yMax = max(gaussSmoothNull)
if max(gaussSmoothPref) > yMax:
    yMax = max(gaussSmoothPref)
plt.plot(gaussSmoothPref, label='pref')   
plt.plot(gaussSmoothNull, label='null') 
plt.plot(gaussSmoothPN, label='p+n')
plt.title('loc0 Pref, loc1 Null: Pref+Pref, Pref+Null, Null+Null')
plt.legend()

plt.xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS],[-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS])
plt.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
plt.ylim([0, yMax*1.1])
plt.xlabel('time (ms)')
plt.ylabel('Firing Rate spikes/sec')
plt.show()
