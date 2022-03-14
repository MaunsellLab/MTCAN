from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

allTrials, header = loadMatFile73('testing_220310_Heatmap_GRF_Spikes.mat')

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
                                    (unitTimeStamps <= stimOnSNEV + 200/1000))
                    spikeCountMat[stCount][eleIndex][aziIndex] = len(stimSpikes[0])
                    stimCount[eleIndex][aziIndex] = stimCount[eleIndex][aziIndex] + 1
                    
                    #histograms
                    histSpikes = np.arange(stimOnSNEV - 0.050, stimOnSNEV + \
                                          (stimDurMS+49)/1000, 0.001)
                    for histCount, i in enumerate(range(len(histSpikes))):
                        if np.around(histSpikes[i],3) in unitTimeStamps:
                            spikeHists[histCount][eleIndex][aziIndex] += 1




spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)

heatMap = plt.imshow(spikeCountMean, cmap='hot', interpolation='nearest')
plt.show()


smoothHist = savgol_filter(spikeHist, 100,3)
plt.plot(smoothHist)
plt.show()

#moving point average
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

histSmooth = smooth(spikeHist,5)
plt.plot(histSmooth)
plt.show()




#insert fake spikes for heatmap

#values for fakespikes
fakeAzi = np.random.randint(numAzi)
fakeEle = np.random.randint(numEle)
fakeSpikes = 50 #spikes/sec
baseFR = 10 #spikes/sec
dt = 1/1000

for currTrial in allTrials:
    trial = currTrial['trial']['data']
    trialEnd = currTrial['trialEnd']['data']
    if trial['instructTrial'] != 1 and trialEnd == 0:
        trialStartMS = currTrial['trialStart']['timeMS']
        trialStartSNEV = currTrial['taskEvents']['trialStart']['timeS']
        trialEndSNEV = currTrial['taskEvents']['trialEnd']['timeS']
        trialLen = trialEndSNEV - trialStartSNEV
        stimDesc = currTrial['stimDesc']
        spikeData = currTrial['spikeData']
        spikeData['channel'] = []
        spikeData['unit'] = []
        spikeData['timeStamp'] = []
        for unit in units:
            channelIdentity = int(unit[0:unit.find('_')])
            spikes = np.random.poisson(trialLen*baseFR)
            spikeTimeS = (np.sort(np.random.choice(np.arange(trialLen * 1000),\
                           spikes, replace = False)))/1000
            spikeData['timeStamp'] = np.append(spikeData['timeStamp'], \
                                     trialStartSNEV + spikeTimeS, 0)
            spikeData['unit'] = np.append(spikeData['unit'], \
                                [unit] * len(spikeTimeS), 0)
            spikeData['channel'] = np.append(spikeData['channel'], \
                                   [channelIdentity] * len(spikeTimeS), 0)
        for count, stim in enumerate(stimDesc['data']):
            if stim['stimType'] == 2:
                break
            if stim['gaborIndex'] == 1:
                aziIndex = int(stim['azimuthIndex'])
                eleIndex = int(stim['elevationIndex'])
                stimOnTimeMS = stimDesc['timeMS'][count]
                stimDiffS = (stimOnTimeMS - trialStartMS)/1000
                stimOnSNEV = trialStartSNEV + stimDiffS
                if aziIndex == fakeAzi and eleIndex == fakeEle:
                    for unit in units:
                        channelIdentity = int(unit[0:unit.find('_')])
                        spikeTimeS = (np.sort(np.random.choice(np.arange(stimDurMS),\
                                    int(np.round(fakeSpikes/(1000/stimDurMS))),\
                                    replace = False)))/1000
                        spikeData['timeStamp'] = np.append(spikeData['timeStamp'],\
                                                stimOnSNEV + spikeTimeS, 0)
                        spikeData['unit'] = np.append(spikeData['unit'], \
                                            [unit] * len(spikeTimeS), 0)
                        spikeData['channel'] = np.append(spikeData['channel'], \
                                               [channelIdentity] * len(spikeTimeS), 0)
