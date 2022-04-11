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


# load my file 

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

spikeCountMat = np.zeros((len(units),30,169))
spikeCountMat[:,0,:] = np.arange(0,169)
spikeCountMat[:,1:,:] = np.nan
stimIndexCount = {}

for currTrial in allTrials: 
    trial = currTrial['trial']['data']
    extendedEOT = currTrial['extendedEOT']['data']
    if trial['instructTrial'] != 1 and 'spikeData' in currTrial and extendedEOT == 0:
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