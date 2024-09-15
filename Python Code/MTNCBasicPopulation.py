# imports
import matplotlib.pyplot as plt
from usefulFns import *
from normalizationFunctions import *
from itertools import combinations
from scipy import stats
import time
import numpy as np
import pandas as pd

######################## MEETZ ########################
# good sessions: 221110, 221115, 221117, 221124, 221128, 221208, 221229, 230123, 230126
# okay sessions: 221010, 221013, 221108, 221206
# bad sessions: 230124

# for loop to run through all good files
fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
            'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
            'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
            'Meetz_230126']

# for loop to run through all good files
# fileList = ['Meetz_221010']

######################## AKSHAN ########################


t0 = time.time()
totUnits = []
totR2 = []
totR2Bram = []
totR2RFWeight = []
totCombs = []
totFilterUnits = []
totFilterR2 = []
totFilterCombs = []
unitIdentifier = []
unitsNormParams = []
alphaLoc0PrefOnly = []
alphaLoc1PrefOnly = []
alphaLoc0Cond = []
alphaLoc1Cond = []
popTunCurve = []
pairSelFromCond = []
pairNPSuppFromCond = []
stimIndex12Cond = [0, 3, 18, 21, 7, 10,
                   25, 28, 14, 17, 32, 35]
popPrefFR = []
highFRAlpha = []
highFRParams = []
medFRAlpha = []
medFRParams = []
lowFRAlpha = []
lowFRParams = []
gaborSigmaSep = []
# arrays for nn, pn, np, pp, and baseline PSTHs
ppArr = []
pnArr = []
npArr = []
nnArr = []
baseArr = []
# arrays for loc 0, 1 p and n PSTHs
p0Arr = []
p1Arr = []
n0Arr = []
n1Arr = []
# arrays for PN, NP z-scored spike counts
pnZscoSpikeCounts = []
npZscoSpikeCounts = []
# array for population heatmap 7x7 grid of resp to loc0,loc1 6 dirs
popRespHeatmap = []
electrodeArr = np.array(np.arange(0, 32)).reshape(16, 2)
pairDistance = []
# ruff and cohen 2016 analysis arr
allNMIAvg = []
alln1n2CorrAvg = []
n1n2NMISimIndx = []
n1n2ElectrodeDiff = []

for file in fileList:

    # Load relevant file here with pyMat reader
    monkeyName, seshDate = file.split('_')
    fileName = f'{monkeyName}_{seshDate}_MTNC_Spikes.mat'
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

    if not os.path.exists('Normalization'):
        os.makedirs('Normalization')
    os.chdir('Normalization/')

    # list of indices of correctTrials (non-instruct, valid trialCertify)
    corrTrials = correctTrialsMTX(allTrials)

    # generate list of unique active units, and their channel
    units = activeUnits('spikeData', allTrials)
    unitCluster = allTrials[corrTrials[0]]['spikeTempInfo']['cgs']
    unitsChannel = unitsInfo(units, corrTrials, allTrials)

    # adds sessions unit to master list of units across sessions
    for unit in units:
        unitIdentifier.append(f'{file}_{unit}')

    # change stimDesc to be list of dictionaries
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        nStim = len(currTrial['stimDesc']['data']['listType'])
        currTrial['stimDesc']['data'] = [{k: v[i] for k, v in
                                          currTrial['stimDesc']['data'].items()}
                                         for i in range(nStim)]

    # assert: stim sequence list is frozen
    seqList = []
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        stimDesc = currTrial['stimDesc']['data']
        for stim in stimDesc:
            if stim['stimLoc'] == 0 and stim['listType'] == 1:
                if len(seqList) < 49:
                    seqList.append(stim['stimIndex'])
                    seqArr = np.array(seqList)
                    lastIndex = stim['stimIndex']
                else:
                    posLastIndex = np.where(seqArr == lastIndex)[0][0]
                    if posLastIndex == len(seqArr)-1:
                        if stim['stimIndex'] != seqArr[0]:
                            print('out of sequence')
                        else:
                            lastIndex = stim['stimIndex']
                    else:
                        if stim['stimIndex'] != seqArr[posLastIndex+1]:
                            print('out of sequence')
                        else:
                            lastIndex = stim['stimIndex']

    # assert: are there correct trials without spikeData
    noSpikeData = []
    for trialCount, currTrial in enumerate(allTrials):
        trial = currTrial['trial']['data']
        extendedEOT = currTrial['extendedEOT']['data']
        if extendedEOT == 0 and trial['instructTrial'] != 1:
            if 'spikeData' not in currTrial:
                noSpikeData.append(trialCount)

    # assert: frame consistency during stimulus duration
    frameRateHz = header['frameRateHz']['data']
    stimDurFrame = []
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        stimDesc = currTrial['stimDesc']['data']
        for stim in stimDesc:
            if stim['stimLoc'] == 0:
                frameDiff = stim['stimOffFrame'] - stim['stimOnFrame']
                stimDurFrame.append(frameDiff)
    if len(set(stimDurFrame)) != 1:
        print('stimulus frame duration not consistent for mapping stimuli')
    else:
        trueStimDurMS = np.int32(np.around(1000/frameRateHz * stimDurFrame[0]))

    # distance between Gabors in the RF
    stimDesc = allTrials[corrTrials[0]]['stimDesc']['data']
    sigmaDeg = allTrials[corrTrials[0]]['rfGabor']['data']['sigmaDeg']
    loc0X, loc0Y, loc1X, loc1Y = 1, 1, 1, 1
    for stim in stimDesc:
        if stim['stimLoc'] == 0:
            loc0X = stim['azimuthDeg']
            loc0Y = stim['elevationDeg']
        if stim['stimLoc'] == 1:
            loc1X = stim['azimuthDeg']
            loc1Y = stim['elevationDeg']
    sep = np.sqrt(((loc0X - loc1X) ** 2) + ((loc0Y - loc1Y) ** 2))
    sigmaSep = sep/sigmaDeg
    gaborSigmaSep.append(sigmaSep)

    # generates a dictionary, numpy array, and Pandas Dataframe of stim Index
    # and corresponding directions/contrasts
    stimIndexDict = {}
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        stimDesc = currTrial['stimDesc']['data']
        for stim in stimDesc:
            if stim['stimLoc'] != 2:
                index = stim['stimIndex']
                if index not in stimIndexDict:
                    stimIndexDict[index] = {}
                    if stim['stimLoc'] not in stimIndexDict[index]:
                        stimIndexDict[index][stim['stimLoc']] = \
                           {'direction': stim['directionDeg'],
                            'contrast': np.around(stim['contrast'], 2)}
                else:
                    if stim['stimLoc'] not in stimIndexDict[index]:
                        stimIndexDict[index][stim['stimLoc']] = \
                           {'direction': stim['directionDeg'],
                            'contrast': np.around(stim['contrast'], 2)}
    stimIndexArray = np.zeros((49, 4))
    for i in range(len(stimIndexDict)):
        stimIndexArray[i][0] = stimIndexDict[i][0]['direction']
        stimIndexArray[i][1] = stimIndexDict[i][0]['contrast']
        stimIndexArray[i][2] = stimIndexDict[i][1]['direction']
        stimIndexArray[i][3] = stimIndexDict[i][1]['contrast']
    stimIndexDF = pd.DataFrame(stimIndexArray, columns=['loc0 Direction',
                                                        'loc0 Contrast',
                                                        'loc1 Direction',
                                                        'loc1 Contrast'])

    # initialize lists/arrays/dataframes for counting spikeCounts and for analysis
    blocksDone = allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone']
    highContrast, zeroContrast = max(stimIndexDF['loc0 Contrast'].unique()), \
                                 min(stimIndexDF['loc0 Contrast'].unique())
    zeroDir = 0
    dirArray = np.array([0, 60, 120, 180, 240, 300])
    loc0IndexArray = np.array([36, 37, 38, 39, 40, 41])
    loc1IndexArray = np.array([42, 43, 44, 45, 46, 47])
    angleMat = np.arange(180, 900, 60)
    spikeCountMat = np.zeros((len(units), blocksDone+1, 49))
    spikeCountLong = []
    onLatency = 50/1000  # time in MS for counting window latency after stim on
    offLatency = 50/1000  # time in MS for counting window latency after stim off
    histPrePostMS = 100  # 100ms window pre/post stimulus on/off
    spikeHists = np.zeros((len(units), 49, trueStimDurMS + (2*histPrePostMS+1)))
    stimIndexCount = np.zeros(49)

    # insert spike counts into matrix of unique stimulus sets
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        if 'spikeData' in currTrial:
            stimDesc = currTrial['stimDesc']['data']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            fixateTimeS = currTrial['taskEvents']['fixate']['time']
            for stim in stimDesc:
                if stim['stimLoc'] == 0 and stim['listType'] == 1:
                    stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'])
                                   / 1000) + stim1TimeS
                    stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'])
                                    / 1000) + stim1TimeS
                    stimIndex = np.int32(stim['stimIndex'])
                    stCount = int(stimIndexCount[stimIndex])
                    stimIndexCount[stimIndex] += 1
                    for unitCount, unit in enumerate(units):
                        if unit in currTrial['spikeData']['unit']:
                            unitIndex = np.where(currTrial['spikeData']['unit'] == unit)[0]
                            unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                            stimSpikes = np.where((unitTimeStamps >= (stimOnTimeS + onLatency)) &
                                                  (unitTimeStamps <= (stimOffTimeS + offLatency)))[0]
                            spikeCountMat[unitCount][stCount][stimIndex] \
                                = len(stimSpikes)
                            spikeCountLong.append([unit,
                                                   stimIndex,
                                                   stimIndexCount[stimIndex],
                                                   len(stimSpikes)])

                            # # # # # FAKE DATA PUT FAKE SPIKES FOR NULL (0˚), PREF(180˚)
                            # # FAKE DATA PN, NP - (0, 180) and (180,0)
                            # if stimIndex == 21 or stimIndex == 39 or stimIndex == 45:
                            #     unitTimeStamps = poissonArrivals(stimOnTimeS, 100*1000/250, 250)
                            # if stimIndex == 0 or stimIndex == 42 or stimIndex == 36:
                            #     unitTimeStamps = poissonArrivals(stimOnTimeS, 10*1000/250, 250)
                            # if stimIndex == 3 or stimIndex == 18:
                            #     unitTimeStamps = poissonArrivals(stimOnTimeS, 55*1000/250, 250)
                            # # # # # # FAKE DATA ABOVE REMOVE

                            # PSTHs
                            stimOnPreSNEV = stimOnTimeS - (histPrePostMS/1000)
                            stimOffPostSNEV = stimOffTimeS + (histPrePostMS/1000)
                            histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                                            & (unitTimeStamps <= stimOffPostSNEV)
                                                             )] - stimOnPreSNEV
                            histStimSpikes = np.int32(histStimSpikes*1000)
                            spikeHists[unitCount][stimIndex][histStimSpikes] += 1

    # mean, SEM, and reshaping of spikeCount matrices
    # create pandas dataframe of spikeCount with corresponding unit, stimIndex
    spikeCountDF = pd.DataFrame(spikeCountLong, columns=['unit', 'stimIndex',
                                                         'stimCount', 'stimSpikes'])
    meanSpike = np.mean(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSD = np.std(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSEM = spikeCountSD/np.sqrt(blocksDone)
    meanSpikeReshaped = np.zeros((len(units), 1, 49))
    for count, i in enumerate(meanSpikeReshaped):
        # row1 of 7x7 grid when reshaped
        i[:, 0:6] = meanSpike[count][0:6]
        i[:, 6] = meanSpike[count][42]
        # row2
        i[:, 7:13] = meanSpike[count][6:12]
        i[:, 13] = meanSpike[count][43]
        # row3
        i[:, 14:20] = meanSpike[count][12:18]
        i[:, 20] = meanSpike[count][44]
        # row4
        i[:, 21:27] = meanSpike[count][18:24]
        i[:, 27] = meanSpike[count][45]
        # row5
        i[:, 28:34] = meanSpike[count][24:30]
        i[:, 34] = meanSpike[count][46]
        # row6
        i[:, 35:41] = meanSpike[count][30:36]
        i[:, 41] = meanSpike[count][47]
        # row7
        i[:, 42:48] = meanSpike[count][36:42]
        i[:, 48] = meanSpike[count][48]

    # Unit's preferred direction based off of gaussian fit of MTNC stimuli
    # filter units based off of gauss fit > 0.90
    angleMat = np.arange(180, 900, 60)
    unitGaussMean = np.zeros(len(units))
    unitGaussSig = np.zeros(len(units))
    filterUnits = []
    for unitCount, unit in enumerate(units):
        loc0Resp = meanSpike[unitCount][36:42]
        loc1Resp = meanSpike[unitCount][42:48]
        combResp = (loc0Resp + loc1Resp) / 2
        extRespMat = np.concatenate((combResp[3:], combResp, combResp[:3]), axis=0)
        maxIndex = np.where(combResp == np.max(combResp))[0][0] + 3
        x = angleMat[maxIndex-3:maxIndex+4]
        y = extRespMat[maxIndex-3:maxIndex+4]
        params = gaussFit(x, y)
        yPred = gauss(x, *params)
        r2 = r2_score(y, yPred)
        if r2 > 0.90:
            filterUnits.append(unit)
            totFilterR2.append(r2)
        unitGaussMean[unitCount] = params[2] % 360
        unitGaussSig[unitCount] = params[3]
        popPrefFR.append(params[1])

    #population analysis for reduced stim set (how to fit normalziation equations)
    # split the 6x6 matrix into all possible 2x2 matrices and run the normalization fit on
    # all the sub-matrices (technically 7x7 matrix with blank condition)
    rfWeightR2 = []
    bramR2 = []
    subsampleMatrix = [i for i in combinations([0, 1, 2, 3, 4, 5], 2)]
    for subpair in subsampleMatrix:
        unitPairedNormFit = []
        unitPairedNormR2 = np.zeros(len(units))
        unitNormFitEstimate = [[] for i in range(len(units))]
        sli = [subpair[0], subpair[1], 6]
        dir1 = dirArray[sli[0]]
        dir2 = dirArray[sli[1]]
        condArr = np.array([dir1, dir2, dir1, dir2])

        for unitCount, unit in enumerate(filterUnits):

            # unitCount = np.where(units == unit)[0][0]
            prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
            nullDirIndex = np.where(dirArray == nullDir)[0][0]
            prefDirIndex = np.where(dirArray == prefDir)[0][0]
            reIndex = (np.array([0, 1, 2, 3, 4, 5]) + nullDirIndex) % 6

            # responses (dependent variables) - matrix of responses
            b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000 / trueStimDurMS

            # fixed (independent) variables - matrix of corresponding stim Indexes
            stimMat, stimMatReIndex = reIndexedStimMat(reIndex)

            # condensed matrices for 2 directions only
            bCondensed = b[sli, :]
            bCondensed = bCondensed[:, sli]
            stimMatCond = stimMat[sli, :]
            stimMatCond = stimMatCond[:, sli]

            # Bram Normalization (L1+L2)/(1+al2+sig) w.o scalar
            guess0 = np.concatenate((bCondensed[2, :-1],
                                     bCondensed[:-1, 2],
                                     [0.2]), axis=0)
            resp = bCondensed.reshape(9)[:-1]
            fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9)[:-1],
                                                    stimIndexDict, dir1, dir2)
            pOpt, pCov, = curve_fit(genNormCondensed, fixedVals, resp.squeeze(),
                                    p0=guess0, maxfev=10000000)

            y_pred = genNormCondensed(fixedVals, *pOpt)
            r2Bram = r2_score(resp.squeeze(), y_pred)

            # rf weight normalization
            fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9),
                                                    stimIndexDict, dir1, dir2)
            resp = bCondensed.reshape(9)
            pOpt, pCov, = curve_fit(rfWeightCondensed, fixedVals, resp.squeeze(),
                                    p0=guess0, maxfev=10000000)

            y_pred = rfWeightCondensed(fixedVals, *pOpt)
            r2 = r2_score(resp.squeeze(), y_pred)


            # bram trick to keep values > 0
            pOpt = pOpt ** 2

            # Append fit parameters for condensed matrix
            totR2.append(r2)
            rfWeightR2.append(r2)
            bramR2.append(r2Bram)
            totR2Bram.append(r2Bram)
            totR2RFWeight.append(r2)
            # alphaLoc0Cond.append(pOpt[4])
            # alphaLoc1Cond.append(pOpt[5])
            unitPairedNormR2[unitCount] = r2
            unitPairedNormFit.append(pOpt)
            unitNormFitEstimate[unitCount] = pOpt

    os.chdir('../../../Python Code')

# plotting Bram vs RF Weight R2
plotData = [totR2Bram, totR2RFWeight]

fig, axes = plt.subplots(1, 1, figsize=(4, 4))

# Plot the violin plots for genNorm arrays
sns.boxplot(data=plotData, ax=axes)
axes.set_title('Box Plot for bramNorm Arrays')
axes.set_xlabel('Position')
axes.set_ylabel('Values')
axes.set_xticks([0, 1])
axes.set_xticklabels(['Bram Norm', 'RF Weight Norm'])
axes.set_ylim(0, 1.2)
plt.show()
