'''
MTN Analysis Script
 - heatmap of Norm responses
 - potentially some correlations

to do:
trial certify
incorp frame render to align spikes
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


# load relevant file 
allTrials, header = loadMatFile73('Meetz', '220509', 'Meetz_220509_MTN_Spikes.mat')


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
                         'contrast': stim['contrast'].tolist()}
                else:
                    if int(stim['stimLoc'].tolist()) not in stimIndexDict[index]:
                        stimIndexDict[index][int(stim['stimLoc'].tolist())] = \
                        {'direction': int(stim['directionDeg'].tolist()),
                         'contrast': stim['contrast'].tolist()}


# numpy array of stimIndexDict
stimIndexArray = np.zeros((169,4))
for i in range(len(stimIndexDict)):
    stimIndexArray[i][0] = stimIndexDict[i][0]['direction']
    stimIndexArray[i][1] = stimIndexDict[i][0]['contrast']
    stimIndexArray[i][2] = stimIndexDict[i][1]['direction']
    stimIndexArray[i][3] = stimIndexDict[i][1]['contrast']

# are there correct trials without spikeData
noSpikeData = []
for trialCount, currTrial in enumerate(allTrials):
    trial = currTrial['trial']['data']
    extendedEOT = currTrial['extendedEOT']['data']
    if extendedEOT == 0 and trial['instructTrial'] != 1:
        if spikeData not in currTrial:
            noSpikeData.append(trialCount)


# generate list of unique active units
units = activeUnits('spikeData', allTrials)

# list of indices of correctTrials (non-instruct, valid trialCertify)
corrTrials = correctTrialsMTX(allTrials)

# insert spike counts into matrix of unique stimulus sets
spikeCountMat = np.zeros((len(units),30,169))
spikeCountMat[:,0,:] = np.arange(0,169)
spikeCountMat[:,1:,:] = np.nan
stimIndexCount = {}
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    if 'spikeData' in currTrial:
        trialStartMS = currTrial['trialStart']['timeMS']
        trialStartSNEV = currTrial['taskEvents']['trialStart']['timeS']
        stimDesc = currTrial['stimDesc']['data']
        for stimCount, stim in enumerate(stimDesc):
            if stim['stimLoc'] == 0 and stim['listType'] == 1:
                stimOnTimeMS = currTrial['stimDesc']['timeMS'][stimCount]
                stimDiffMS = stimOnTimeMS - trialStartMS
                stimOnSNEV = trialStartSNEV + (stimDiffMS / 1000)
                stimIndex = np.int32(stim['stimIndex'])
                stimIndexCount[stimIndex] = stimIndexCount.get(stimIndex, 0) + 1
                for unitCount, unit in enumerate(units):
                    if unit in currTrial['spikeData']['unit']:
                        unitIndex = np.where(currTrial['spikeData']['unit'] == unit)
                        unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnSNEV) & 
                                    (unitTimeStamps <= stimOnSNEV + 493/1000))
                        spikeCountMat[unitCount][stimIndexCount[stimIndex]][stimIndex] \
                        = len(stimSpikes[0])

meanSpike = np.nanmean(spikeCountMat[:,1:,:], axis = 1)

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

# heatmap of correlations
for unitCount in range(len(units)):
    a = meanSpikeReshaped[unitCount]
    b = a.reshape(13,13)

    #using seaborn
    ax = sns.heatmap(b)
    ax.set_xticks(np.arange(13)+0.5)
    ax.set_xticklabels(['0','60','120','180','240','300','0','60','120',
                        '180','240','300','blank'], rotation = 45)
    ax.xaxis.set_ticks_position("top")

    ax.set_yticks(np.arange(13)+0.5)
    ax.set_yticklabels(['0','60','120','180','240','300','0','60','120',
                        '180','240','300','blank'], rotation = 0)
    plt.show()


#correlations
combs = [i for i in combinations(units, 2)]
corrMat = np.zeros((len(combs),1,169))

for count, i in enumerate(combs):
    n1 = units.index(i[0])
    n2 = units.index(i[1])
    for j in range(np.shape(spikeCountMat)[2]):
        stimCorr = ma.corrcoef(ma.masked_invalid(spikeCountMat[n1,1:,j]),\
                               ma.masked_invalid(spikeCountMat[n2,1:,j]))
        corrMat[count][0][j] = stimCorr[0][1]


b = stats.zscore(corrMatMean, axis=0, nan_policy='omit')

corrMatMean = np.mean(corrMat, axis = (1,0))
popCorrMean = ma.mean(ma.masked_invalid(corrMatMean))

print(ma.mean(ma.masked_invalid(corrMat[0][0])))
print(ma.mean(ma.masked_invalid(corrMat[1][0])))
print(ma.mean(ma.masked_invalid(corrMat[2][0])))

# selectivity for direction (similar to Bram)
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

# tuning similarity b/w neurons 
dirTuningMat = np.load('unitsDirTuningMat.npy') #load from directions folder
extTunMat = np.concatenate((dirTuningMat[:,3:], dirTuningMat[:,:], 
                            dirTuningMat[:,:3], axis=1))
angleMat = np.arange(180,900,60)

combs = [i for i in combinations(units, 2)]
pairSimScore = np.zeros((len(combs),1))
for pairCount, pair in enumerate(combs):
    n1, n2 = units.index(pair[0]), units.index(pair[1])
    n1Max = int(np.where(dirTuningMat[n1] == np.max(dirTuningMat[n1]))[0] + 3)
    n1X = angleMat[n1Max-3:n1Max+4]
    n1Y = extTunMat[n1][n1Max-3:n1Max+4]



    n2Max = int(np.where(dirTuningMat[n2] == np.max(dirTuningMat[n2]))[0] + 3)
    n2X = angleMat[n2Max-3:n2Max+4]
    n2Y = extTunMat[n2][n2Max-3:n2Max+4]





tc = np.array([10,23,45,80,37,16,12])
x = np.array([0,60,120,180,240,300,360])
x_full = np.linspace(0, 360, 1000)
params = gauss_fit(x, tcNorm)
y_full_fit = gauss(x_full, *params)

plt.plot(x_full, y_full_fit, '--r', label='fit')
plt.scatter(x, tcNorm, label= 'not fit')
plt.legend()




