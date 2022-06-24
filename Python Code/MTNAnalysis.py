'''
MTN Analysis Script
 - heatmap of Norm responses
 - potentially some correlations

to do:
trial certify
rotate norm plots to incorporate the preferrred direction in the middle
convert heatmap to spikes/sec, it's at spikes/stimDurMS
heatmap range starts at 0

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
allTrials, header = loadMatFile73('Meetz', '220622', 'Meetz_220622_MTN_Spikes.mat')


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


# insert spike counts into matrix of unique stimulus sets
stimDurMS = int(header['stimDurationMS']['data'].tolist())
frameRateHz = header['frameRateHz']['data'].tolist()
spikeCountMat = np.zeros((len(units),30,169))
spikeCountMat[:,0,:] = np.arange(0,169)
histPrePostMS = 100
spikeHists = np.zeros((len(units),169, trueStimDurMS+2*histPrePostMS+12))
# spikeCountMat[:,1:,:] = np.nan
stimIndexCount = {}
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    if 'spikeData' in currTrial:
        stimDesc = currTrial['stimDesc']['data']
        stim1TimeS = currTrial['taskEvents']['stimulusOn'][0]['time'].tolist()
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

                        #PSTHs
                        stimOnPreSNEV = stimOnTimeS - histPrePostMS/1000
                        stimOffPostSNEV = stimOffTimeS + histPrePostMS/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)\
                                    & (unitTimeStamps <= stimOffPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        spikeHists[unitCount, stimIndex, histStimSpikes] += 1

# blocks Done 
blocksDone = int(allTrials[corrTrials[-1]]['blockStatus']['data']['blocksDone']
                 .tolist()) 

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

## heatmap of correlations
for unitCount in range(len(units)):
    a = meanSpikeReshaped[unitCount]
    b = a.reshape(13,13)
    bSmooth = gaussian_filter(b, sigma=1)

    #split 13x13 into smaller grids, so different combs(LL, LH,HL, HH) of
    # high and low contrasts are separate grids 
    b0L1L = bSmooth[:6,:6]
    b0H1L = bSmooth[:6,6:12]
    b0L1H = bSmooth[6:12,:6]
    b0H1H = bSmooth[6:12,6:12]
    b1Blank = bSmooth[12,:12]
    b0Blank = bSmooth[:12,12]
    b01Blank = bSmooth[12,12]

    #using seaborn
    ax = sns.heatmap(bSmooth, square=True, linewidths=0.2, vmin=0)
    ax.set_xticks(np.arange(13)+0.5)
    ax.set_title(f'heatmap of normalization for {units[unitCount]}')
    ax.set_xticklabels(['0','60','120','180','240','300','0','60','120',
                        '180','240','300','blank'], rotation = 45)
    ax.xaxis.set_ticks_position("top")

    ax.set_yticks(np.arange(13)+0.5)
    ax.set_yticklabels(['0','60','120','180','240','300','0','60','120',
                        '180','240','300','blank'], rotation = 0)
    plt.tight_layout()
    plt.savefig(f'{units[unitCount]}.pdf')
    plt.close('all')


## correlations incomplete 
combs = [i for i in combinations(units, 2)]
corrMat = np.zeros((len(combs),169))

# z-scored spikeCountMat
zSpikeCountMat = stats.zscore(spikeCountMat[:,1:blocksDone+1,:], axis=2, nan_policy='omit')
zSpikeCountMat = np.nan_to_num(zSpikeCountMat)
for count, i in enumerate(combs):
    n1 = np.where(units == i[0])[0][0]
    n2 = np.where(units == i[1])[0][0]
    for j in range(np.shape(spikeCountMat)[2]):
        stimCorr = stats.pearsonr(zSpikeCountMat[n1,1:blocksDone+1,j],
                                 zSpikeCountMat[n2,1:blocksDone+1,j])
        corrMat[count][j] = stimCorr[0]

popCorr = np.mean(np.nanmean(corrMat,axis=1))


# scipy curveFit Normalization parameters
def func(fixed, L_0, L_60, L_120, L_180, L_240, L_300, aL0, aL1, sig):
    c0,c1,l0,l1 = fixed
    L = np.array([L_0, L_60, L_120, L_180, L_240, L_300])
    return((c0*L*l0).sum(-1) + (c1*L*l1).sum(-1))/((aL0*c0[:,0])+(aL1*c1[:,0])+sig)

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
        #append([c0,l0_oh,c1,l1_oh])
    c0s = np.array(c0s)
    c1s = np.array(c1s)
    l0s = np.array(l0s)
    l1s = np.array(l1s)
    #fix = np.array(fix)
    #print(fix.shape)
    pOpt, pCov = curve_fit(func, (c0s, c1s, l0s, l1s), resp, bounds=(
        (0,0,0,0,0,0,0,0,0),(np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,1,1,1)))
    print(unit,pOpt)


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

#PSTHs for P, P+-60, P+-120, etc
yMax = 0
unitPref = spikeHists[2,115,:] * 1000/stimIndexCount[115]
unitI1 = spikeHists[2,114,:] * 1000/stimIndexCount[114]
unitI2 = spikeHists[2,116,:] * 1000/stimIndexCount[116]
unitI3 = spikeHists[2,117,:] * 1000/stimIndexCount[117]
unitI4 = spikeHists[2,118,:] * 1000/stimIndexCount[118]
unitI5 = spikeHists[2,119,:] * 1000/stimIndexCount[119]
unitNull = spikeHists[2,136,:] * 1000/stimIndexCount[136]
gaussSmoothPref = gaussian_filter1d(unitPref, 10)
gaussSmoothI1 = gaussian_filter1d(unitI1, 10)
gaussSmoothI2 = gaussian_filter1d(unitI2, 10)
gaussSmoothI3 = gaussian_filter1d(unitI3, 10)
gaussSmoothI4 = gaussian_filter1d(unitI4, 10)
gaussSmoothI5 = gaussian_filter1d(unitI5, 10)
gaussSmoothNull = gaussian_filter1d(unitNull, 10)
if max(gaussSmoothPref) > yMax:
    yMax = max(gaussSmoothPref)
plt.plot(gaussSmoothPref, label='pref')   
plt.plot(gaussSmoothNull, label='null') 
plt.plot(gaussSmoothI1, label='P, I1')
plt.plot(gaussSmoothI2, label='P, I2')
plt.plot(gaussSmoothI3, label='P, I3')
plt.plot(gaussSmoothI4, label='P, N')
plt.plot(gaussSmoothI5, label='P, I5')
plt.title('P+P, N+N, P+I')
plt.legend()
plt.xticks([0,histPrePostMS,histPrePostMS+trueStimDurMS,2*histPrePostMS+trueStimDurMS],[-(histPrePostMS), 0, 0+trueStimDurMS, trueStimDurMS+histPrePostMS])
plt.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS, color='grey', alpha=0.2)
plt.ylim([0, yMax*1.1])
plt.xlabel('time (ms)')
plt.ylabel('Firing Rate spikes/sec')
plt.show()

'''
'''

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


## direction tuning similarity b/w neurons 
dirTuningMat = np.load('unitsDirTuningMat.npy') #load from directions folder
extTunMat = np.concatenate((dirTuningMat[:,3:], dirTuningMat[:,:], 
                            dirTuningMat[:,:3]), axis=1)
angleMat = np.arange(180,900,60)

combs = [i for i in combinations(units, 2)]
pairSimScore = np.zeros((len(combs),1))
for pairCount, pair in enumerate(combs):
    n1, n2 = units.index(pair[0]), units.index(pair[1])
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

    # bhattacharyya similarity score 
    BC = bhattCoef(m1, m2, v1, v2)
    pairSimScore[pairCount] = BC