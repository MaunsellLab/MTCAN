'''
MTN Analysis Script
 - heatmap of normalized responses
 - potentially some correlations
 - normalization equation, parameter fitting
 - PSTHs of normalized responses


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

### testing
allTrials, header = loadMatFile73('Testing', 'Meetz_220621', 'Meetz_220621_MTN.mat')


# load relevant file 
allTrials, header = loadMatFile73('Meetz', '220622_3', 'Meetz_220622_MTN_Spikes.mat')


if not os.path.exists('Normalization'):
    os.makedirs('Normalization')
os.chdir('Normalization/')


# list of indices of correctTrials (non-instruct, valid trialCertify)
corrTrials = correctTrialsMTX(allTrials)

# generate list of unique active units
units = activeUnits('spikeData', allTrials)


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


# generates a dictionary of stim Index and corresponding directions/contrasts
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


# numpy array of stimIndexDict
stimIndexArray = np.zeros((169,4))
for i in range(len(stimIndexDict)):
    stimIndexArray[i][0] = stimIndexDict[i][0]['direction']
    stimIndexArray[i][1] = stimIndexDict[i][0]['contrast']
    stimIndexArray[i][2] = stimIndexDict[i][1]['direction']
    stimIndexArray[i][3] = stimIndexDict[i][1]['contrast']

# pandas dataframe
stimIndexDF = pd.DataFrame(stimIndexArray, columns=['loc0 Direction', 'loc0 Contrast',
                                            'loc1 Direction', 'loc1 Contrast'])


blocksDone = int(allTrials[corrTrials[-1]]['blockStatus']['data']['blocksDone']
                 .tolist()) 
lowContrast,highContrast = stimIndexDF['loc0 Contrast'].unique()[0], \
                           stimIndexDF['loc0 Contrast'].unique()[1]
dirArray = np.array([0,60,120,180,240,300,0,60,120,180,240,300])
spikeCountMat = np.zeros((len(units),blocksDone+2,169))
spikeCountMat[:,0,:] = np.arange(0,169)
spikeCountDF = pd.DataFrame(columns = ['unit', 'stimIndex', 'stimCount', 'spikeCount'])
fixSpikeCount = np.zeros((len(units),len(corrTrials),1))
histPrePostMS = 100
spikeHists = np.zeros((len(units),169, trueStimDurMS+2*histPrePostMS+12))
# spikeCountMat[:,1:,:] = np.nan
stimIndexCount = {}

# insert spike counts into matrix of unique stimulus sets
for trialCount, corrTrial in enumerate(corrTrials):
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
                stimIndexCount[stimIndex] = stimIndexCount.get(stimIndex, 0) + 1
                for unitCount, unit in enumerate(units):
                    if unit in currTrial['spikeData']['unit']:
                        unitIndex = np.where(currTrial['spikeData']['unit'] == unit)[0]
                        unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnTimeS) & 
                                    (unitTimeStamps <= stimOffTimeS))[0]
                        spikeCountMat[unitCount][stimIndexCount[stimIndex]][stimIndex] \
                        = len(stimSpikes)
                        spikeCountDF = spikeCountDF.append({'unit':unit, 
                        'stimIndex':stimIndex, 'stimCount':stimIndexCount[stimIndex],
                        'spikeCount' :len(stimSpikes)}, ignore_index = True)

                        # fixation spikeCount
                        fixSpikes = np.where((unitTimeStamps >= fixateTimeS) & 
                                             (unitTimeStamps <= stim1TimeS))[0]
                        fixSpikeCount[unitCount][trialCount][0] = len(fixSpikes)

                        #PSTHs
                        stimOnPreSNEV = stimOnTimeS - histPrePostMS/1000
                        stimOffPostSNEV = stimOffTimeS + histPrePostMS/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)\
                                    & (unitTimeStamps <= stimOffPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        spikeHists[unitCount, stimIndex, histStimSpikes] += 1

        

meanSpike = np.mean(spikeCountMat[:,1:blocksDone+1,:], axis = 1)
meanSpikeReshaped = np.zeros((len(units),1,169))
for count,i in enumerate(meanSpikeReshaped):
    i[:,0:6] = meanSpike[count][0:6]
    i[:,6:12] = meanSpike[count][36:42]
    i[:,12] = meanSpike[count][156]
    i[:,13:19] = meanSpike[count][6:12]
    i[:,19:25] = meanSpike[count][42:48]
    i[:,25] = meanSpike[count][157]
    i[:,26:32] = meanSpike[count][12:18]
    i[:,32:38] = meanSpike[count][48:54]
    i[:,38] = meanSpike[count][158]
    i[:,39:45] = meanSpike[count][18:24]
    i[:,45:51] = meanSpike[count][54:60]
    i[:,51] = meanSpike[count][159]
    i[:,52:58] = meanSpike[count][24:30]
    i[:,58:64] = meanSpike[count][60:66]
    i[:,64] = meanSpike[count][160]
    i[:,65:71] = meanSpike[count][30:36]
    i[:,71:77] = meanSpike[count][66:72]
    i[:,77] = meanSpike[count][161]
    i[:,78:84] = meanSpike[count][72:78]
    i[:,84:90] = meanSpike[count][108:114]
    i[:,90] = meanSpike[count][162]
    i[:,91:97] = meanSpike[count][78:84]
    i[:,97:103] = meanSpike[count][114:120]
    i[:,103] = meanSpike[count][163]
    i[:,104:110] = meanSpike[count][84:90]
    i[:,110:116] = meanSpike[count][120:126]
    i[:,116] = meanSpike[count][164]
    i[:,117:123] = meanSpike[count][90:96]
    i[:,123:129] = meanSpike[count][126:132]
    i[:,129] = meanSpike[count][165]
    i[:,130:136] = meanSpike[count][96:102]
    i[:,136:142] = meanSpike[count][132:138]
    i[:,142] = meanSpike[count][166]
    i[:,143:149] = meanSpike[count][102:108]
    i[:,149:155] = meanSpike[count][138:144]
    i[:,155] = meanSpike[count][167]
    i[:,156:168] = meanSpike[count][144:156]
    i[:,168] = meanSpike[count][168]

## heatmap of normalization
for unit in range(len(units)):
    a = meanSpikeReshaped[unit]
    b = a.reshape(13,13)
    bSmooth = gaussian_filter(b, sigma=1)
    maxLoc0 = max(bSmooth[12,:])
    maxLoc1 = max(bSmooth[:,12])
    if maxLoc0 > maxLoc1:
        prefDir = dirArray[np.where(bSmooth[12,:]==maxLoc0)[0][0]]
    else:
        prefDir = dirArray[np.where(bSmooth[:,12]==maxLoc1)[0][0]]
    nullDir = (prefDir + 180)%360
    nullIndex = np.where(dirArray==nullDir)[0][0]
    reIndex = (np.array([0,1,2,3,4,5,0,1,2,3,4,5])+nullIndex) % 6

    bSmoothReIndex = bSmooth[:12,:12]
    bSmoothReIndex = bSmoothReIndex[:,reIndex]
    bSmoothReIndex = bSmoothReIndex[reIndex,:]
    b0Blank = bSmooth[12,:12]
    b0Blank = b0Blank[reIndex]
    b1Blank = bSmooth[:12,12]
    b1Blank = b1Blank[reIndex]
    b0L1L = bSmoothReIndex[:6,:6]
    b0L1H = bSmoothReIndex[:6,6:12]
    b0H1L = bSmoothReIndex[6:12,:6]
    b0H1H = bSmoothReIndex[6:12,6:12]
    
    gridspec_kw = {"height_ratios":[6,6,1], "width_ratios" : [6,6,1]}
    heatmapkws = dict(square=False, cbar=False, linewidths=1.0, vmin=0, vmax=np.max(bSmooth))
    asp = bSmooth.shape[0]/float(bSmooth.shape[1])
    figW = 8
    figH = figW*asp
    fig, axes = plt.subplots(ncols=3, nrows=3, figsize=(figW, figH), gridspec_kw=gridspec_kw)
    # left = 0.07; right=0.7
    # bottom = 0.1; top = 0.1
    # plt.subplots_adjust(left=left,right=right,bottom=bottom,top=top,wspace=0.1,hspace=0.1*asp)
    plt.subplots_adjust(wspace=0.1,hspace=0.1*asp)
    sns.heatmap(b0L1L, ax=axes[0,0], xticklabels=True, yticklabels=True, **heatmapkws)
    sns.heatmap(b0L1H, ax=axes[0,1], xticklabels=True, yticklabels=False, **heatmapkws)
    sns.heatmap(b0H1L, ax=axes[1,0], xticklabels=False, yticklabels=True, **heatmapkws)
    sns.heatmap(b0H1H, ax=axes[1,1], xticklabels=False, yticklabels=False, **heatmapkws)
    sns.heatmap(b0Blank[np.newaxis,:6], ax=axes[2,0], xticklabels=False, yticklabels=False, **heatmapkws)
    sns.heatmap(b0Blank[np.newaxis,6:12],ax=axes[2,1], xticklabels=False, yticklabels=False, **heatmapkws)
    sns.heatmap(b1Blank[:6,np.newaxis], ax=axes[0,2], xticklabels=False, yticklabels=False, **heatmapkws)
    sns.heatmap(b1Blank[6:12,np.newaxis],ax=axes[1,2], xticklabels=False, yticklabels=False, **heatmapkws)
    sns.heatmap(bSmooth[12:,12:],ax=axes[2,2], xticklabels=False, yticklabels=False, **heatmapkws)
    
    axes[0,0].xaxis.set_ticks_position("top")
    axes[0,1].xaxis.set_ticks_position("top")
    axes[0,0].set_xticklabels(dirArray[:6][reIndex[:6]], rotation = 45)
    axes[0,1].set_xticklabels(dirArray[:6][reIndex[:6]], rotation = 45)
    axes[0,0].set_yticklabels(dirArray[:6][reIndex[:6]], rotation = 0)
    axes[1,0].set_yticklabels(dirArray[:6][reIndex[:6]], rotation = 0)
    axes[0,0].set_xlabel('loc 0 low contrast')
    axes[0,0].xaxis.set_label_position('top') 
    axes[0,0].set_ylabel('loc 1 low contrast')
    axes[1,0].set_ylabel('loc 1 high contrast')
    axes[0,1].set_xlabel('loc 0 high contrast')
    axes[0,1].xaxis.set_label_position('top') 
    plt.tight_layout()
    plt.show()

    #using seaborn
    ax = sns.heatmap(bSmooth, square=True, linewidths=0.2, vmin=0)
    # ax = plt.imshow(b, cmap='gist_heat', interpolation='bilinear')
    ax.set_xticks(np.arange(13)+0.5)
    ax.set_title(f'heatmap of normalization for {units[unit]}')
    ax.set_xticklabels(['0','60','120','180','240','300','0','60','120',
                        '180','240','300','blank'], rotation = 45)
    ax.xaxis.set_ticks_position("top")

    ax.set_yticks(np.arange(13)+0.5)
    ax.set_yticklabels(['0','60','120','180','240','300','0','60','120',
                        '180','240','300','blank'], rotation = 0)
    plt.tight_layout()
    plt.savefig(f'{unit}.pdf')
    plt.close('all')


## scipy curveFit Normalization parameters
def func(fixed, L_0, L_60, L_120, L_180, L_240, L_300, sL0, sL1, aL, sig):
    c0,c1,l0,l1 = fixed
    L = np.array([L_0, L_60, L_120, L_180, L_240, L_300])
    # return((c0*sL0*L*l0).sum(-1) + (c1*sL1*L*l1).sum(-1))/((aL0*c0[:,0])+(aL1*c1[:,0])+sig)
    # return (((c0*sL0*(l0*L)).sum(-1) + (c1*sL1*(l1*L)).sum(-1))/((aL0*c0[:,0])+(aL1*c1[:,0])+sig))
    return (((c0*sL0*(l0*L)).sum(-1) + (c1*sL1*(l1*L)).sum(-1))/((aL*c0[:,0])+(aL*c1[:,0])+sig))
for unitCount, unit in enumerate(units):
    resp = np.reshape(spikeCountMat[unitCount][1:blocksDone+1,:],(169*blocksDone))
    fixParam = np.tile(np.arange(169), blocksDone)

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
# high high contrast curve_fit
highCIndex = stimIndexDF.index[(stimIndexDF['loc0 Contrast'] == 1) 
             & (stimIndexDF['loc1 Contrast'] == 1)]
for unitCount, unit in enumerate(units):
    respSubMat = spikeCountMat[unitCount,1:blocksDone+1,108:144]
    resp = np.reshape(respSubMat,(36*blocksDone))
    fixParam = np.tile(np.arange(108,144), blocksDone)

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
        (0,0,0,0,0,0,0,0,0,0,0),(np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,1,1,1,1,1)))
    print(unit,pOpt)


## PSTHs for P+P, N+N, P+I, and converse for other location
for unitCount, unit in enumerate(units):
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
    ax2 = []
    for row in range(2):
        for col in range(4):
            ax2.append(plt.subplot2grid((2,4), (row,col)))
    ax2 = np.array(ax2)
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
            smoothPlot = gaussian_filter1d(dirPlot,8)
            if max(smoothPlot) > yMax:
                yMax = max(smoothPlot)
            ax2[plotCount].plot(smoothPlot, label=f'{loc0Dir}+{loc1Dir}')
            ax2[plotCount].set_title(f'loc0 contrast {loc0Con} loc1 contrast {loc1Con}', fontsize= 5)
            ax2[plotCount].set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
            ax2[plotCount].set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
            ax2[plotCount].axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.1)
            if plotCount == 0 or plotCount == 1:
                ax2[plotCount].legend(loc='upper left', prop={'size': 4})
            subCount += 1
            if subCount % 7 == 0:
                plotCount += 1
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    for i in range(8):
        ax2[i].set_ylim([0,yMax*1.1])
    plt.show()


## mean response of PP, PN, NN to low and high contrast
for unitCount, unit in enumerate(units):
    b = meanSpikeReshaped[unitCount].reshape(13,13)
    bSmooth = gaussian_filter(b, sigma=1)
    maxLoc0 = max(bSmooth[12,:])
    maxLoc1 = max(bSmooth[:,12])
    if maxLoc0 > maxLoc1:
        prefDir = dirArray[np.where(bSmooth[12,:]==maxLoc0)[0][0]]
    else:
        prefDir = dirArray[np.where(bSmooth[:,12]==maxLoc1)[0][0]]
    nullDir = (prefDir + 180)%360
    print(unit, prefDir, nullDir)
    otherDir = [i for i in dirArray[:6] if i != prefDir and i != nullDir]

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
    plt.show()


## correlations incomplete 
combs = [i for i in combinations(units, 2)]

#filtered units test
combs = [i for i in combinations(filterUnits, 2)]
corrMat = np.zeros((len(combs),169))

# z-scored spikeCountMat
zSpikeCountMat = stats.zscore(spikeCountMat[:,1:blocksDone+1,:], axis=1, nan_policy='omit')
zSpikeCountMat = np.nan_to_num(zSpikeCountMat)
for count, i in enumerate(combs):
    n1 = np.where(units == i[0])[0][0]
    n2 = np.where(units == i[1])[0][0]
    for j in range(np.shape(spikeCountMat)[2]):
        # stimCorr = stats.pearsonr(zSpikeCountMat[n1,1:blocksDone+1,j],
        #                          zSpikeCountMat[n2,1:blocksDone+1,j])
        stimCorr = stats.pearsonr(spikeCountMat[n1,1:blocksDone+1,j],
                                 spikeCountMat[n2,1:blocksDone+1,j])
        corrMat[count][j] = stimCorr[0]

popCorr = np.mean(np.nanmean(corrMat,axis=1))


#fixation correlations
fixCorrMat = np.zeros(len(combs))
zfixSpikeCount = stats.zscore(fixSpikeCount, axis = 1)
for count, i in enumerate(combs):
    n1 = np.where(units == i[0])[0][0]
    n2 = np.where(units == i[1])[0][0]
    fixCorr = stats.pearsonr(zfixSpikeCount[n1].flatten(),zfixSpikeCount[n2].flatten())
    fixCorrMat[count] = fixCorr[0]

## RF location tuning similarity b/w neurons Bhattacharyya Distance 2D
RFLocMat = np.load('../RFLoc Tuning/unitsRFLocMat.npy')
combs = [i for i in combinations(units, 2)]
pairLocSimScore = np.zeros((len(combs),1))
for pairCount, pair in enumerate(combs):
    n1 = np.where(units == pair[0])[0][0]
    n2 = np.where(units == pair[1])[0][0]

    m1, cov1, p = gauss2dParams(RFLocMat[n1])
    m2, cov2, p2 = gauss2dParams(RFLocMat[n2])
    BC = bhattCoef2D(m1,m2,cov1,cov2)
    pairLocSimScore[pairCount] = BC


## direction tuning similarity b/w neurons Bhattacharyya Distance
dirTuningMat = np.load('../Direction Tuning/unitsDirTuningMat.npy') #load from directions folder
unitsBaselineMat = np.load('../Direction Tuning/unitsBaselineMat.npy')
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
    params = gauss_fit(n1X, n1Y)
    # n1YFull = gauss(n1XFull, *params)
    m1 = params[2] # mean neuron 1
    v1 = params[3]**2 # var neuron 1
    n1TrueMean = m1 - 360
    
    n2Max = int(np.where(dirTuningMat[n2] == np.max(dirTuningMat[n2]))[0] + 3)
    n2X = angleMat[n2Max-3:n2Max+4]
    # n2Y = extTunMat[n2][n2Max-3:n2Max+4]
    n2Y = extTunMat[n2][n2Max-3:n2Max+4]/max(extTunMat[n2][n2Max-3:n2Max+4])
    n2XFull = np.linspace(n2X[0], n2X[-1],1000)
    params = gauss_fit(n2X, n2Y)
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
unitSelectivity = np.zeros((len(units),169)) 
unitSelectivity[:,:] = np.nan
for uCount, unit in enumerate(units):
    for stim in range(169):
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


'''
Normalization across grid PSTHs plots - working
    for x in contrastArray:
    loc0Con, loc1Con = x
    subPlot = plt.subplot2grid((2,4), subPlotList[subPlotCount])
    # loc0 = Pref : loc1 = Null, Intermediate
    # pref pref
    histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == prefDir) & 
        (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == prefDir)
        & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
    dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
    subPlot.plot(gaussian_filter1d(dirPlot,10), label='pref')
    subPlot.set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
    subPlot.set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
    subPlot.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
    # null null
    histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == nullDir) & 
        (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == nullDir)
        & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
    dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
    subPlot.plot(gaussian_filter1d(dirPlot,10), label='null')
    subPlot.set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
    subPlot.set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
    subPlot.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
    # pref null
    histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == prefDir) & 
        (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == nullDir)
        & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
    dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
    subPlot.plot(gaussian_filter1d(dirPlot,10), label='pref+null')
    subPlot.set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
    subPlot.set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
    subPlot.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
    #intermediate Dir
    for i in otherDir:
        histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == prefDir) & 
        (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == i)
        & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
        dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
        subPlot.plot(gaussian_filter1d(dirPlot,10), label=f'pref+I{prefDir,i}')
        subPlot.set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
        subPlot.set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
        subPlot.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
    subPlotCount += 1

    subPlot = plt.subplot2grid((2,4), subPlotList[subPlotCount])
    #loc0 = Null, Intermediate : loc1 = Pref
    histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == prefDir) & 
        (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == prefDir)
        & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
    dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
    subPlot.plot(gaussian_filter1d(dirPlot,10), label='pref')
    subPlot.set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
    subPlot.set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
    subPlot.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
    # null null
    histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == nullDir) & 
        (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == nullDir)
        & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
    dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
    subPlot.plot(gaussian_filter1d(dirPlot,10), label='null')
    subPlot.set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
    subPlot.set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
    subPlot.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
    # null pref
    histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == nullDir) & 
        (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == prefDir)
        & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
    dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
    subPlot.plot(gaussian_filter1d(dirPlot,10), label='null+pref')
    subPlot.set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
    subPlot.set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
    subPlot.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
    #intermediate Dirs
    for i in otherDir:
        histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == i) & 
        (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == prefDir)
        & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
        dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
        subPlot.plot(gaussian_filter1d(dirPlot,10), label=f'I+prefDir{i, prefDir}')
        subPlot.set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
        subPlot.set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
        subPlot.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
    subPlotCount += 1


        ax2 = []
    for row in range(2):
        for col in range(4):
            ax2.append(plt.subplot2grid((2,4), (row,col)))
    ax2 = np.array(ax2)

    fig = plt.figure()
    fig.set_size_inches(8,4)
    subPlotList = [(0,0),(0,1), (0,2), (0,3), (1,0), (1,1), (1,2), (1,3)]
    plotCount = 0
    for x in contrastList:
        subCount = 0
        loc0Con, loc1Con = x
        ax2[plotCount] = plt.subplot2grid((2,4), subPlotList[plotCount])
        for locDir in locDirList:
            loc0Dir, loc1Dir = locDir
            histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) & 
            (stimIndexDF['loc0 Contrast'] == loc0Con) & (stimIndexDF['loc1 Direction'] == loc1Dir)
            & (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
            dirPlot = spikeHists[unitCount,histIndex,:] * 1000/stimIndexCount[histIndex]
            ax2[plotCount].plot(gaussian_filter1d(dirPlot,8), label=f'{loc0Dir}+{loc1Dir}')
            ax2[plotCount].set_xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS])
            ax2[plotCount].set_xticklabels([-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS], fontsize=5)
            ax2[plotCount].axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.1)
            subCount += 1
            if subCount == 7:
                plotCount += 1
                ax2[plotCount] = plt.subplot2grid((2,4), subPlotList[plotCount])
            if subCount == 14:
                plotCount += 1
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    for i in range(8):
        ax2[i].set_ylim([0,300])
    plt.show()
'''
