# imports
import matplotlib.pyplot as plt
from usefulFns import *
from itertools import combinations
from scipy import stats
import time
import numpy as np
import pandas as pd
from normalizationFunctions import *

########################################### MEETZ ###########################################
# good sessions: 221110, 221115, 221117, 221124, 221128, 221208, 221229, 230123, 230126
# okay sessions: 221010, 221013, 221108, 221206
# bad sessions: 230124

# fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
#             'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
#             'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
#             'Meetz_230126']

########################################### AKSHAN ###########################################

# good sessions: 240927, 240930, 241016, 241017
# okay sessions: 240826, 241002 (13 blocks), 241021 (17 blocks)
# bad sessions: 240827, 240828

# fileList = ['Akshan_240826', 'Akshan_240927', 'Akshan_240930, 'Akshan_241002',
#             'Akshan_241016', 'Akshan_241017', 'Akshan_241021']




########################################### Both ###########################################
# # for loop to run through all files
fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
            'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
            'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
            'Meetz_230126', 'Akshan_240826', 'Akshan_240927', 'Akshan_240930',
            'Akshan_241002', 'Akshan_241016', 'Akshan_241017', 'Akshan_241021']

t0 = time.time()
totUnits = []
totR2 = []
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

for fileIterator in fileList:

    # Load relevant file here with pyMat reader
    monkeyName, seshDate = fileIterator.split('_')
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
        unitIdentifier.append(f'{fileIterator}_{unit}')

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
    offLatency = 100/1000  # time in MS for counting window latency after stim off
    histPrePostMS = 100  # 100ms window pre/post stimulus on/off
    spikeHists = np.zeros((len(units), 49, trueStimDurMS + (2*histPrePostMS+1)))
    stimIndexCount = np.zeros(49)
    filterUnits = []

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
    highFRUnits = []
    lowFRUnits = []
    tightUnits = []
    wideUnits = []
    daysPrefFR = []
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
        if r2 > 0.90:  # working version
        # if r2 > 0.75:
            filterUnits.append(unit)
            totFilterR2.append(r2)
        unitGaussMean[unitCount] = params[2] % 360
        unitGaussSig[unitCount] = params[3]
        popPrefFR.append(params[1])
        daysPrefFR.append(params[1])
        if params[3] < 40:
            tightUnits.append(unit)
        else:
            wideUnits.append(unit)
        # if params[1] > 5.775:
        #     highFRUnits.append(unit)
        # else:
        #     lowFRUnits.append(unit)

    for unitCount, unit in enumerate(units):
        if daysPrefFR[unitCount] > np.median(daysPrefFR):
            highFRUnits.append(unit)
        else:
            lowFRUnits.append(unit)

    totUnits.append(len(units))
    totFilterUnits.append(len(filterUnits))

    # the different combinations of neuron pairs from total units
    combs = [i for i in combinations(units, 2)]
    totCombs.append(len(combs))
    # filtered pairs for pairs of units having R2 > 0.75
    filterCombs = [i for i in combinations(filterUnits, 2)]
    totFilterCombs.append(len(filterCombs))
    # combinations of high FR units
    highFRCombs = [i for i in combinations(highFRUnits, 2)]
    # combinations of low FR units
    lowFRCombs = [i for i in combinations(lowFRUnits, 2)]
    # combinations of units with tight tuning
    tightCombs = [i for i in combinations(tightUnits, 2)]
    # combinations of units with wide tuning
    wideCombs = [i for i in combinations(wideUnits, 2)]
    # combinations of units separated by more than 250ums
    distCombs = []
    # for i in combs:  #### working version
    for i in filterCombs:
        n1 = unitsChannel[np.where(units == i[0])[0][0]]
        n2 = unitsChannel[np.where(units == i[1])[0][0]]
        n1ElectrodePos = np.where(electrodeArr == n1)[0][0]
        n2ElectrodePos = np.where(electrodeArr == n2)[0][0]
        pairDistOnElectrode = abs(n1ElectrodePos - n2ElectrodePos)
        pairDistance.append(pairDistOnElectrode * 50)
        # if pairDistOnElectrode >= 5: ### working version
        if pairDistOnElectrode >= 5:
            distCombs.append(i)

    # get p, n, p+n psth for every unit - normalized and aligned
    # to each unit's pref direction. This will also extract spike counts
    # for PN, NP condition to create distribution of variance from the mean (z-score)
    sli = [0, 3, 6]
    for unitCount, unit in enumerate(units):

        # unitCount = np.where(units == unit)[0][0]
        # find direction tested that is closest to the pref dir
        # and reindex around this so that it is in the middle of the grid
        prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])

        nullDirIndex = np.where(dirArray == nullDir)[0][0]
        reIndex = (np.array([0, 1, 2, 3, 4, 5]) + nullDirIndex) % 6

        # raw data reindex to have null in the top left corner
        b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS
        bReIndex = reIndexedRespMat(b, reIndex)
        bReIndexNorm = bReIndex / np.max(bReIndex)
        popRespHeatmap.append(bReIndexNorm.reshape(49))

        # matrix of indices for b and bReIndexed
        stimMat, stimMatReIndex = reIndexedStimMat(reIndex)

        # condensed matrices for 2 directions only
        stimMatCond = stimMatReIndex[sli, :]
        stimMatCond = stimMatCond[:, sli]

        # matrix of each unit's normalized resp to p,n, p+n (9 conditions)
        unitNormHist = []
        yMax = 0
        for count, i in enumerate(stimMatCond.reshape(9)):
            dirPlot = spikeHists[unitCount, int(i), :] * 1000 / stimIndexCount[int(i)]
            smoothPlot = gaussian_filter1d(dirPlot, 5)
            if max(smoothPlot) > yMax:
                yMax = max(smoothPlot)
            unitNormHist.append(dirPlot)
            if count == 1:
                tempMat = spikeCountMat[unitCount, :blocksDone, int(i)]
                tempMatZsco = stats.zscore(tempMat)
                pnZscoSpikeCounts.append(tempMatZsco)
            if count == 3:
                tempMat = spikeCountMat[unitCount, :blocksDone, int(i)]
                tempMatZsco = stats.zscore(tempMat)
                npZscoSpikeCounts.append(tempMatZsco)

        # unitNormHist = spikeHists[unitCount, stimMatCond.reshape(9).astype(int), :] *1000
        # unitNormHist /= stimIndexCount[stimMatCond.reshape(9).astype(int)].reshape(9, 1)
        # yMax = unitNormHist.max()
        # unitNormHist /= yMax
        # assert(unitNormHist.max() == 1)

        # normalize response
        unitNormHist = np.array(unitNormHist) / yMax

        # trial shuffle PN, NP condition
        # if np.random.rand(1) >= 0.5:
        #     pnArr.append(unitNormHist[1])  # normal condition
        #     npArr.append(unitNormHist[3])  # normal condition
        # else:
        #     pnArr.append(unitNormHist[3])  # shuffling
        #     npArr.append(unitNormHist[1])  # shuffling
        # if sigmaSep > 4.300 and sigmaSep < 4.53:  # median split for gabor separation
        pnArr.append(unitNormHist[1])
        npArr.append(unitNormHist[3])
        nnArr.append(unitNormHist[0])
        ppArr.append(unitNormHist[4])
        p0Arr.append(unitNormHist[7])
        p1Arr.append(unitNormHist[5])
        n0Arr.append(unitNormHist[6])
        n1Arr.append(unitNormHist[2])
        baseArr.append(unitNormHist[8])

    # initialize lists for paired Gabors
    pairPairedCorr = []
    pairDoubleCov = []
    pairDoubleSD = []
    pairSelectivityIndex = []
    pairNonPrefSuppIndex = []
    pairNMIIndex = []
    pairSimpleNMIIndex = []
    pairRawSelectivityIndex = []
    pairNewSelectivityIndex = []
    # initialize lists for single Gabor presentations
    pairSingleCorr = []
    pairSingleCov = []
    pairSingleSD =[]
    pairSingleSelectivityIndex = []
    pairSingleNonPrefSuppIndex = []
    pairSingleNMI = []
    pairSingleSimpleNMI = []
    pairSingleRawSelectivity = []
    singleNewSelectivity = []
    # initialize lists for blank Gabor presentations
    pairBlankCorr = []
    # initialize lists for ruff and cohen 2016 (NMI vs corr - zscored)
    dayNMI = []
    dayCorr = []

    # # subsampling the 6x6 matrix of directions into unique 2x2 matrices
    subsampleMatrix = [i for i in combinations([0, 1, 2, 3, 4, 5], 2)]
    for subpair in subsampleMatrix:
        unitPairedNormFit = []
        unitPairedNormR2 = np.zeros(len(units))
        unitNormFitEstimate = [[] for i in range(len(units))]
        sli = [subpair[0], subpair[1], 6]
        dir1 = dirArray[sli[0]]
        dir2 = dirArray[sli[1]]
        condArr = np.array([dir1, dir2, dir1, dir2])

        for unitCount, unit in enumerate(units):

            # unitCount = np.where(units == unit)[0][0]
            prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
            nullDirIndex = np.where(dirArray == nullDir)[0][0]
            prefDirIndex = np.where(dirArray == prefDir)[0][0]
            reIndex = (np.array([0, 1, 2, 3, 4, 5]) + nullDirIndex) % 6

            # responses (dependent variables) - matrix of responses
            b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS

            # fixed (independent) variables - matrix of corresponding stim Indexes
            stimMat, stimMatReIndex = reIndexedStimMat(reIndex)

            # condensed matrices for 2 directions only
            bCondensed = b[sli, :]
            bCondensed = bCondensed[:, sli]
            stimMatCond = stimMat[sli, :]
            stimMatCond = stimMatCond[:, sli]

            # RF Weighted Normalization (L1+L2)/(1+w+sig)
            guess0 = np.concatenate((bCondensed[2, :-1],
                                     bCondensed[:-1, 2],
                                     [0.2], [0.5]), axis=0)
            resp = bCondensed.reshape(9)[:-1]
            fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9)[:-1],
                                                    stimIndexDict, dir1, dir2)
            pOpt, pCov, = curve_fit(rfWeightedNormCondensed, fixedVals, resp.squeeze(),
                                    maxfev=10000000)

            y_pred = rfWeightedNormCondensed(fixedVals, *pOpt)
            r2 = r2_score(resp.squeeze(), y_pred)
            print(unit, r2)

            # bram trick to keep values > 0
            pOpt = pOpt ** 2

            # Append fit parameters for condensed matrix
            totR2.append(r2)

            unitPairedNormR2[unitCount] = r2
            unitPairedNormFit.append(pOpt)
            unitNormFitEstimate[unitCount] = pOpt

            # reindex the single presentation conditions to align to direction
            # that is closest to unit's preferred direction

            prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
            nullDirIndex = np.where(dirArray == nullDir)[0][0]
            reIndex = (np.array([0, 1, 2, 3, 4, 5])+nullDirIndex) % 6
            loc1ReIndex = b[:6, 6][reIndex]
            loc0ReIndex = b[6, :6][reIndex]
            respAvg = ((loc1ReIndex + loc0ReIndex) / 2)
            baselineNormalized = b[6, 6] / np.max(respAvg)
            respNormalized = respAvg / np.max(respAvg)
            respNormalized = np.concatenate((respNormalized,
                                             respNormalized[:1],
                                             [np.array(unitGaussMean[unitCount])],
                                             [np.array(baselineNormalized)]),
                                            axis=0).tolist()
            popTunCurve.append(respNormalized)

        # generate paired stimulus correlation, selectivity index,
        # suppression, and NMI index for each gabor pair for each pair
        # of neurons (similar to Bram fig 2c)
        tempMat = stats.zscore(spikeCountMat[:, :blocksDone, :], axis=1, nan_policy='omit')
        stimIndexMat = np.arange(36).reshape(6, 6)
        upperTriangle = upperTriMasking(stimIndexMat)

        for pairCount, pair in enumerate(distCombs):   ### working version
        # for pairCount, pair in enumerate(filterCombs):
            n1 = np.where(units == pair[0])[0][0]
            n2 = np.where(units == pair[1])[0][0]

            # unit's preferred directions
            n1PrefDir = unitGaussMean[n1]
            n2PrefDir = unitGaussMean[n2]

            # ruff and cohen 2016 nmi vs corr reproduction
            n1NMIArr = []
            n2NMIArr = []
            n1n2Corr = []

            # indices and correlations for paired Gabor stimuli
            for count, i in enumerate(np.int_(stimMatCond[:2, :2].reshape(4))):

                # correlation for that Gabor pair b/w 2 units excluding trials where
                # spike counts exceeded 3 SD from mean
                n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
                pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                # extract directions of the gabor pair
                loc0Dir = stimIndexDict[i][0]['direction']
                loc1Dir = stimIndexDict[i][1]['direction']

                # n1 NMI and selectivity (raw spike counts)
                b = meanSpikeReshaped[n1].reshape(7, 7) * 1000/trueStimDurMS
                bCondensed = b[sli, :]
                bCondensed = bCondensed[:, sli]
                pairedMat = bCondensed[:2, :2].reshape(4)
                loc0Single = bCondensed[2, :2]
                loc1Single = bCondensed[:2, 2]
                loc0SingleResp = loc0Single[np.where(condArr == loc0Dir)[0][0]]
                loc1SingleResp = loc1Single[np.where(condArr == loc1Dir)[0][0]]
                rp, rn = max([loc0SingleResp, loc1SingleResp]), min(
                             [loc0SingleResp, loc1SingleResp])
                rx = pairedMat[count]
                # n1NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[count]) / (
                #       (loc0SingleResp + loc1SingleResp) + pairedMat[count])
                n1NMI = (rp - rx) / (rp + rx)
                n1SimpleNMI = (loc0SingleResp + loc1SingleResp) / pairedMat[count]
                n1RawSelectivity = (loc0SingleResp - loc1SingleResp) / (
                                    loc0SingleResp + loc1SingleResp)
                # if n1RawSelectivity >= 0:
                #     n1NMI = 0.5 * (-n1NMI + 1)
                # else:
                #     n2NMI = 0.5 * (n1NMI + 1)

                # n2 NMI and selectivity (raw spike counts)
                b = meanSpikeReshaped[n2].reshape(7, 7) * 1000/trueStimDurMS
                bCondensed = b[sli, :]
                bCondensed = bCondensed[:, sli]
                pairedMat = bCondensed[:2, :2].reshape(4)
                loc0Single = bCondensed[2, :2]
                loc1Single = bCondensed[:2, 2]
                loc0SingleResp = loc0Single[np.where(condArr == loc0Dir)[0][0]]
                loc1SingleResp = loc1Single[np.where(condArr == loc1Dir)[0][0]]
                rp, rn = max([loc0SingleResp, loc1SingleResp]), min(
                             [loc0SingleResp, loc1SingleResp])
                rx = pairedMat[count]
                # n2NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[count]) / (
                #       (loc0SingleResp + loc1SingleResp) + pairedMat[count])
                n2NMI = (rp - rx) / (rp + rx)
                n2SimpleNMI = (loc0SingleResp + loc1SingleResp) / pairedMat[count]
                n2RawSelectivity = (loc0SingleResp - loc1SingleResp) / (
                                    loc0SingleResp + loc1SingleResp)
                # if n2RawSelectivity >= 0:
                #     n2NMI = 0.5 * (-n2NMI + 1)
                # else:
                #     n2NMI = 0.5 * (n2NMI + 1)

                # bram selectivity recreation
                # n1 selectivity and suppression index
                loc0Resp = unitNormFitEstimate[n1][np.where(condArr == loc0Dir)[0][0]]
                loc1Resp = unitNormFitEstimate[n1][np.where(condArr == loc1Dir)[0][0] + 2]
                al0 = 1
                al1 = unitNormFitEstimate[n1][4]
                n1Selectivity, n1NonPrefSupp = getSelAndSuppIndx(loc0Resp, loc1Resp, al0, al1)

                # n2 selectivity and suppression index
                loc0Resp = unitNormFitEstimate[n2][np.where(condArr == loc0Dir)[0][0]]
                loc1Resp = unitNormFitEstimate[n2][np.where(condArr == loc1Dir)[0][0] + 2]
                al0 = 1
                al1 = unitNormFitEstimate[n2][4]
                n2Selectivity, n2NonPrefSupp = getSelAndSuppIndx(loc0Resp, loc1Resp, al0, al1)

                # # pair selectivity, suppression, and NMI index
                pairSelectivity = (np.sign(n1Selectivity) * np.sign(n2Selectivity) *
                                   np.sqrt(abs(n1Selectivity) * abs(n2Selectivity)))
                pairNewSelectivity = 1
                # pairSuppression = (abs(n1NonPrefSupp - n2NonPrefSupp)) / (
                #                       n1NonPrefSupp + n2NonPrefSupp)
                pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)

                pairNMI = (n1NMI + n2NMI) / 2
                pairSimpleNMI = abs(n1SimpleNMI - n2SimpleNMI) / (
                                    n1SimpleNMI + n2SimpleNMI)
                pairRawSelectivity = (np.sign(n1RawSelectivity) *
                                      np.sign(n2RawSelectivity) *
                                      np.sqrt(abs(n1RawSelectivity) *
                                              abs(n2RawSelectivity)))

                pairPairedCorr.append(pairStimCorr)
                pairDoubleCov.append(pairDCov[0][1])
                pairDoubleSD.append(pairDSD)
                pairSelectivityIndex.append(pairSelectivity)
                pairNewSelectivityIndex.append(pairNewSelectivity)
                pairNonPrefSuppIndex.append(pairSuppression)
                pairNMIIndex.append(pairNMI)
                pairSimpleNMIIndex.append(pairSimpleNMI)
                pairRawSelectivityIndex.append(pairRawSelectivity)

                # ruff and cohen 2016
                n1n2Corr.append(pairStimCorr)
                n1NMIArr.append(n1SimpleNMI)
                n2NMIArr.append(n2SimpleNMI)

                if i in stimIndex12Cond:
                    pairSelFromCond.append(pairSelectivity)
                    pairNPSuppFromCond.append(pairSuppression)

                # loc 0 single Gabor corr, excluding trials where spike counts > 3 SD
                stimIndex = np.int_(stimMatCond[2, :][np.where(condArr == loc0Dir)[0][0]])

                # # INSERT FAKE DATA FOR LOC 0 SINGLE GABOR
                # spikeCountMat[n1, :blocksDone, stimIndex], spikeCountMat[n2, :blocksDone, stimIndex] = generateTwoCorrArrays(blocksDone, 0.03)
                # # FAKE DATA ABOVE REMOVE DURING REAL ANALYSIS

                n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
                singleStimLoc0Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                pairSingleCov.append(pairSCov[0][1])
                pairSingleSD.append(pairSSD)
                pairSingleCorr.append(singleStimLoc0Corr)
                pairSingleSelectivityIndex.append(pairSelectivity)
                pairSingleNonPrefSuppIndex.append(pairSuppression)
                pairSingleSimpleNMI.append(pairSimpleNMI)
                pairSingleNMI.append(pairNMI)
                pairSingleRawSelectivity.append(pairRawSelectivity)
                singleNewSelectivity.append(pairNewSelectivity)

                # loc 1 single Gabor corr, excluding trials where spike counts > 3 SD
                stimIndex = np.int_(stimMatCond[:, 2][np.where(condArr == loc1Dir)[0][0]])

                # # INSERT FAKE DATA FOR LOC 1 SINGLE GABOR
                # spikeCountMat[n1, :blocksDone, stimIndex], spikeCountMat[n2, :blocksDone, stimIndex] = generateTwoCorrArrays(blocksDone, 0.03)
                # # FAKE DATA ABOVE REMOVE DURING REAL ANALYSIS

                n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
                singleStimLoc1Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                pairSingleCov.append(pairSCov[0][1])
                pairSingleSD.append(pairSSD)
                pairSingleCorr.append(singleStimLoc1Corr)
                pairSingleSelectivityIndex.append(pairSelectivity)
                pairSingleNonPrefSuppIndex.append(pairSuppression)
                pairSingleSimpleNMI.append(pairSimpleNMI)
                pairSingleNMI.append(pairNMI)
                pairSingleRawSelectivity.append(pairRawSelectivity)
                singleNewSelectivity.append(pairNewSelectivity)

                # pair blank corr
                blankIndex = 48
                skipTrials = []
                n1SpikeMat = spikeCountMat[n1, :blocksDone, blankIndex]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, blankIndex]
                blankCorr, blankCov, blankSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                pairBlankCorr.append(blankCorr)

            n1n2NMIAvg = np.mean(np.array([n1NMIArr, n2NMIArr]).flatten())
            n1n2CorrAvg = np.nanmean(n1n2Corr)

            # n1NMI n2NMI similarity index (n1-n2/n1+n2)
            n1n2NMISim = abs(np.mean(n1NMIArr) - np.mean(n2NMIArr)) / (
                             np.mean(n1NMIArr) + np.mean(n2NMIArr))

            n1Chan = unitsChannel[n1]
            n2Chan = unitsChannel[n2]
            n1ElectrodePos = np.where(electrodeArr == n1Chan)[0][0]
            n2ElectrodePos = np.where(electrodeArr == n2Chan)[0][0]
            pairDistOnElectrode = abs(n1ElectrodePos - n2ElectrodePos)

            n1n2NMISimIndx.append(n1n2NMISim)
            n1n2ElectrodeDiff.append(pairDistOnElectrode * 50)
            dayNMI.append(n1n2NMIAvg)
            dayCorr.append(n1n2CorrAvg)

    # ruff and cohen 2016
    dayNMIZscore = stats.zscore(dayNMI)
    dayCorrZscore = stats.zscore(dayCorr)
    for i in range(len(dayNMIZscore)):
        allNMIAvg.append(dayNMIZscore[i])
        alln1n2CorrAvg.append(dayCorrZscore[i])

    pairDoubleCov = np.array(pairDoubleCov)
    pairDoubleSd = np.array(pairDoubleSD)
    pairPairedCorr = np.array(pairPairedCorr)
    pairSelectivityIndex = np.array(pairSelectivityIndex)
    pairNonPrefSuppIndex = np.array(pairNonPrefSuppIndex)
    pairNMIIndex = np.array(pairNMIIndex)
    pairSimpleNMIIndex = np.array(pairSimpleNMIIndex)
    pairRawSelectivityIndex = np.array(pairRawSelectivityIndex)
    pairNewSelectivityIndex = np.array(pairNewSelectivityIndex)

    pairSingleCov = np.array(pairSingleCov)
    pairSingleSD = np.array(pairSingleSD)
    pairSingleCorr = np.array(pairSingleCorr)
    pairSingleSelectivityIndex = np.array(pairSingleSelectivityIndex)
    pairSingleNonPrefSuppIndex = np.array(pairSingleNonPrefSuppIndex)
    pairSingleSimpleNMI = np.array(pairSingleSimpleNMI)
    pairSingleNMI = np.array(pairSingleNMI)
    pairSingleRawSelectivity = np.array(pairSingleRawSelectivity)
    singleNewSelectivity = np.array(singleNewSelectivity)

    pairBlankCorr = np.array(pairBlankCorr)

    # compile similarity scores and correlation matrix into one matrix
    gaborPairCompiledArr = np.array([pairPairedCorr,
                                     pairSelectivityIndex,
                                     pairNonPrefSuppIndex,
                                     pairNMIIndex,
                                     pairSimpleNMIIndex,
                                     pairRawSelectivityIndex,
                                     pairDoubleCov,
                                     pairDoubleSD,
                                     pairBlankCorr,
                                     pairNewSelectivityIndex])
    gaborSingleCompiledArr = np.array([pairSingleCorr,
                                       pairSingleSelectivityIndex,
                                       pairSingleNonPrefSuppIndex,
                                       pairSingleSimpleNMI,
                                       pairSingleNMI,
                                       pairSingleRawSelectivity,
                                       pairSingleCov,
                                       pairSingleSD,
                                       singleNewSelectivity])
    np.save(f'../../gaborSingleCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborSingleCompiledArr)
    np.save(f'../../gaborPairCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborPairCompiledArr)
    os.chdir('../../../Python Code')

print(time.time()-t0)

