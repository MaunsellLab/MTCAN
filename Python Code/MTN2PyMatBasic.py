# imports
import matplotlib.pyplot as plt
from usefulFns import *
from itertools import combinations
from scipy import stats
import time
import numpy as np
import pandas as pd

########################################### MEETZ ###########################################
# good sessions: 221110, 221115, 221117, 221124, 221128, 221208, 221229, 230123, 230126
# okay sessions: 221010, 221013, 221108, 221206
# bad sessions: 230124

# # for loop to run through all files
# fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
#             'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
#             'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
#             'Meetz_230126']

########################################### AKSHAN ###########################################

# for loop to run through all files
# fileList = ['Akshan_240826']


########################################### Both ###########################################
# for loop to run through all files
fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
            'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
            'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
            'Meetz_230126', 'Akshan_240826']


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
        if r2 > 0.90:
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
    for i in combs:  #### working version
    # for i in filterCombs:
        n1 = unitsChannel[np.where(units == i[0])[0][0]]
        n2 = unitsChannel[np.where(units == i[1])[0][0]]
        n1ElectrodePos = np.where(electrodeArr == n1)[0][0]
        n2ElectrodePos = np.where(electrodeArr == n2)[0][0]
        pairDistOnElectrode = abs(n1ElectrodePos - n2ElectrodePos)
        pairDistance.append(pairDistOnElectrode * 50)
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

        # # FAKE DATA TEST
        # prefDir = 180
        # nullDir = 0
        # # FAKE DATA TEST

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

    # temp mats for ranking each unit's alpha based off FR for each axis
    # alphaRankMat = np.zeros((len(units), 6, 3))
    # alphaRankMat = np.zeros((len(units), 7, 3))

    # # for loop to run through the 3 orthogonal axes
    # for orthoAxis in range(3):
    #     # SCIPY CURVEFIT/Minimize
    #     # scipy curvefit Normalization parameters will also generate population
    #     # tuning average aligned to preferred direction
    #     unitPairedNormFit = []
    #     unitPairedNormR2 = np.zeros(len(units))
    #     unitNormFitEstimate = [[] for i in range(len(units))]
    #
    #     # # put null as first element
    #     # prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
    #     # nullDirIndex = np.where(dirArray == nullDir)[0][0]
    #     # temp = [orthoAxis, orthoAxis+3]
    #     # tempDiff = abs(temp - nullDirIndex)
    #     # nullIndex = np.where(tempDiff == min(tempDiff))[0][0]
    #     # prefIndex = np.where(tempDiff == max(tempDiff))[0][0]
    #     # sli = [temp[nullIndex], temp[prefIndex], 6]
    #
    #     sli = [orthoAxis, orthoAxis+3, 6]
    #     dir1 = dirArray[sli[0]]
    #     dir2 = dirArray[sli[1]]
    #     condArr = np.array([dir1, dir2, dir1, dir2])

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

            # Generic Normalization (L1+L2)/(1+al2+sig) w.o scalar
            guess0 = np.concatenate((bCondensed[2, :-1],
                                     bCondensed[:-1, 2],
                                     [0.2]), axis=0)
            resp = bCondensed.reshape(9)[:-1]
            fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9)[:-1],
                                                    stimIndexDict, dir1, dir2)

            # # bram implementation
            # lb = np.sqrt([0.01, 0.01, 0.01, 0.01, 0.1])
            # ub = np.sqrt([20, 20, 20, 20, 1])
            # params = []
            # sse = []
            # attempts = 0
            # while attempts < 10:
            #     startPoint = (ub - lb) * np.random.uniform() + lb
            #     guess = guess0 + startPoint
            #     pOpt, pCov, = curve_fit(genNormCondensed, fixedVals, resp.squeeze(),
            #                             p0=guess, maxfev=10000000)
            #     y_pred = genNormCondensed(fixedVals, *pOpt)
            #     print(r2_score(resp.squeeze(), y_pred))
            #     attemptSSE = np.sum((y_pred - resp.squeeze()) ** 2)
            #     params.append(pOpt ** 2)
            #     sse.append(attemptSSE)
            #     attempts += 1

            # pOpt, pCov, = curve_fit(genNormCondensed, fixedVals, resp.squeeze(),
            #                         p0=guess0,
            #                         bounds=((0, 0, 0, 0, 0),
            #                                 (np.inf, np.inf, np.inf, np.inf,
            #                                  6)))
            pOpt, pCov, = curve_fit(genNormCondensed, fixedVals, resp.squeeze(),
                                    p0=guess0, maxfev=10000000)

            y_pred = genNormCondensed(fixedVals, *pOpt)
            r2 = r2_score(resp.squeeze(), y_pred)
            print(unit, r2)

            # bram trick to keep values > 0
            pOpt = pOpt ** 2

            # Append fit parameters for condensed matrix
            totR2.append(r2)
            # alphaLoc0Cond.append(pOpt[4])
            # alphaLoc1Cond.append(pOpt[5])
            unitPairedNormR2[unitCount] = r2
            unitPairedNormFit.append(pOpt)
            unitNormFitEstimate[unitCount] = pOpt

            # if subpair[0] == nullDirIndex and subpair[1] == prefDirIndex:
            #     alphaLoc0PrefOnly.append(pOpt[4])
            #     alphaLoc1PrefOnly.append(pOpt[5])
            # if subpair[0] == prefDirIndex and subpair[1] == nullDirIndex:
            #     alphaLoc0PrefOnly.append(pOpt[4])
            #     alphaLoc1PrefOnly.append(pOpt[5])

            # alphaRankMat[unitCount, :, orthoAxis] = pOpt

            # # fit using scipy.optimize.minimize for condensed mat
            # p0 = np.concatenate((bCondensed[2, :-1], bCondensed[:-1, 2],
            #                      [0.0, 0.05]), axis=0)
            # bnds = ((0, None), (0, None), (0, None),
            #         (0, None), (0, 5), (0, None))
            # resp = bCondensed.reshape(9)
            # fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9),
            #                                         stimIndexDict, dir1, dir2)
            #
            # # # res with bounds
            # res = scipy.optimize.minimize(driverFuncCondensend, p0, args=(fixedVals, resp),
            #                               method='Nelder-Mead', bounds=bnds,
            #                               options={'maxfev': 1*(10**20),
            #                                        'fatol': 0.11,
            #                                        'maxiter': 10000,
            #                                        'xatol': 0.0001})
            # y_pred = genNormCondensed(fixedVals, *res.x)
            # r2_1 = r2_score(resp.squeeze(), y_pred)
            #
            # res2 = scipy.optimize.minimize(driverFuncCondensend1, p0, args=(fixedVals, resp),
            #                               method='Nelder-Mead', bounds=bnds,
            #                               options={'maxfev': 1*(10**20),
            #                                        'fatol': 0.11,
            #                                        'maxiter': 10000,
            #                                        'xatol': 0.0001})
            # y_pred = genNormCondensed1(fixedVals, *res2.x)
            # r2_2 = r2_score(resp.squeeze(), y_pred)
            #
            # if r2_1 > r2_2:
            #     unitPrefSite.append(0)
            #     param = res.x
            #     r2 = r2_1
            # else:
            #     unitPrefSite.append(1)
            #     param = res2.x
            #     r2 = r2_2
            #
            # print(r2)
            # # res without bounds
            # res = scipy.optimize.minimize(driverFuncCondensend, p0, args=(fixedVals, resp),
            #                               method='Nelder-Mead',
            #                               options={'maxfev': 1*(10**20),
            #                                        'fatol': 0.11,
            #                                        'maxiter': 10000,
            #                                        'xatol': 0.0001})


            # # Append fit parameters for condensed matrix
            # totR2.append(r2)
            # # alphaLoc0Cond.append(res.x[4])
            # # alphaLoc1Cond.append(res.x[5])
            # unitPairedNormR2[unitCount] = r2
            # unitPairedNormFit.append(param)
            # unitNormFitEstimate[unitCount] = param
            # unitsNormParams.append(param)

            # # append fit params to alpha ranking mat
            # alphaRankMat[unitCount, :, orthoAxis] = res.x

            # if r2 > 0.75:
            #     filterUnits.append(unit)
            #     totFilterR2.append(r2)

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

        for pairCount, pair in enumerate(distCombs):
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

    # # alphaRankMat reorder to rank alpha by response rate
    # maxResp = alphaRankMat[:, :4, :].max(1)
    # maxRespSorted = np.argsort(maxResp, 1)
    # sortedAlphas = np.array([alphaRankMat[i, 4, si] for i, si in enumerate(maxRespSorted)])
    # sortedParams = np.array([alphaRankMat[i, :, si] for i, si in enumerate(maxRespSorted)])
    # highFRAlpha.append(sortedAlphas[:, 2].tolist())
    # medFRAlpha.append(sortedAlphas[:, 1].tolist())
    # lowFRAlpha.append(sortedAlphas[:, 0].tolist())
    # for i in range(len(sortedParams)):
    #     highFRParams.append(sortedParams[i, 2, :].tolist())
    #     medFRParams.append(sortedParams[i, 1, :].tolist())
    #     lowFRParams.append(sortedParams[i, 0, :].tolist())

print(time.time()-t0)

"""
                 STOP HERE
"""

"""
extra code

# following code is for selectivity based off of neuron tuning 
            # n1 selectivity and suppression index
            loc0Dir = stimIndexDict[i][0]['direction']
            loc1Dir = stimIndexDict[i][1]['direction']
            prefDir = unitGaussMean[n1]
            loc0Resp = unitNormFitEstimate[n1][np.where(dirArray == loc0Dir)[0][0]]
            loc1Resp = unitNormFitEstimate[n1][np.where(dirArray == loc1Dir)[0][0] + 6]
            if loc0Resp > loc1Resp:
                if abs(loc0Dir - prefDir) > 180:
                    if loc0Dir > prefDir:
                        loc0Dir = prefDir - (360 - (loc0Dir - prefDir))
                    else:
                        prefDir = loc0Dir - (360 - (prefDir - loc0Dir))
                n1Selectivity = 1 - (abs(loc0Dir - prefDir) / 180)
                n1NonPrefSupp = unitNormFitEstimate[n1][13] / (
                    unitNormFitEstimate[n1][13] + unitNormFitEstimate[n1][12])
                n1PrefResp = 0
            else:
                if abs(loc1Dir - prefDir) > 180:
                    if loc1Dir > prefDir:
                        loc1Dir = prefDir - (360 - (loc1Dir - prefDir))
                    else:
                        prefDir = loc1Dir - (360 - (prefDir - loc1Dir))
                n1Selectivity = 1 - (abs(loc1Dir - prefDir) / 180)
                n1NonPrefSupp = unitNormFitEstimate[n1][12] / (
                    unitNormFitEstimate[n1][13] + unitNormFitEstimate[n1][12])
                n1PrefResp = 1
            if n1Selectivity < 0.5:
                n1Selectivity = -n1Selectivity
            else:
                n1Selectivity = n1Selectivity - 0.5

            # n2 selectivity and suppression index
            loc0Dir = stimIndexDict[i][0]['direction']
            loc1Dir = stimIndexDict[i][1]['direction']
            prefDir = unitGaussMean[n2]
            loc0Resp = unitNormFitEstimate[n2][np.where(dirArray == loc0Dir)[0][0]]
            loc1Resp = unitNormFitEstimate[n2][np.where(dirArray == loc1Dir)[0][0] + 6]
            if loc0Resp > loc1Resp:
                if abs(loc0Dir - prefDir) > 180:
                    if loc0Dir > prefDir:
                        loc0Dir = prefDir - (360 - (loc0Dir - prefDir))
                    else:
                        prefDir = loc0Dir - (360 - (prefDir - loc0Dir))
                n2Selectivity = 1 - (abs(loc0Dir - prefDir) / 180)
                n2NonPrefSupp = unitNormFitEstimate[n2][13] / (
                    unitNormFitEstimate[n2][13] + unitNormFitEstimate[n2][12])
                n2PrefResp = 0
            else:
                if abs(loc1Dir - prefDir) > 180:
                    if loc1Dir > prefDir:
                        loc1Dir = prefDir - (360 - (loc1Dir - prefDir))
                    else:
                        prefDir = loc1Dir - (360 - (prefDir - loc1Dir))
                n2Selectivity = 1 - (abs(loc1Dir - prefDir) / 180)
                n2NonPrefSupp = unitNormFitEstimate[n2][12] / (
                    unitNormFitEstimate[n2][13] + unitNormFitEstimate[n2][12])
                n2PrefResp = 1
            if n2Selectivity < 0.5:
                n2Selectivity = -n2Selectivity
            else:
                n2Selectivity = n2Selectivity - 0.5

            loc0Dir = stimIndexDict[i][0]['direction']
            loc1Dir = stimIndexDict[i][1]['direction']
"""

# Plotting Code

# distribution of spike counts around PN, NP mean
pnZscoSpikeCounts = [ele for i in pnZscoSpikeCounts for ele in i]
npZscoSpikeCounts = [ele for i in npZscoSpikeCounts for ele in i]
plt.hist(pnZscoSpikeCounts, bins=50); plt.show()
plt.hist(npZscoSpikeCounts, bins=50); plt.show()


# SUPERPLOT OF PSTH, Normalization (pref+null+blank) heatmap, and bar plot
# Create subset of Pandas DF for P, N, P+N conditions for all units
pnContrastPlotDF = pd.DataFrame(columns=['unit', 'stimIndex', 'stimCount',
                                'stimSpikes', 'contrast', 'prefNullStr'])
for unitCount, unit in enumerate(units):
    # find direction tested that is closest to unit's preferred and null direction
    prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
    orientCount = 0
    orientList = ['pref+pref', 'pref+null', 'null+pref', 'null+null',
                 'pref+blank', 'null+blank','blank+pref', 'blank+null']
    for j in [(prefDir, prefDir), (prefDir, nullDir),
              (nullDir, prefDir), (nullDir, nullDir),]:
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

for unitCount, unit in enumerate(units):
    b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS
    bSmooth = gaussian_filter(b, sigma=1)
    prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])

    ## figure
    # fig = plt.figure(constrained_layout=True)
    fig = plt.figure()
    fig.set_size_inches(17, 7)

    gs0 = gridspec.GridSpec(1, 5)

    heights = [1, 1, 1]
    widths = [4, 4, 4]
    gs00 = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs0[0],
                                            width_ratios=widths,
                                            height_ratios=heights)
    gs01 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[1:3])
    gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[3])
    gs03 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[4])

    # PSTHs Plot
    locDirList = [(nullDir, nullDir), (prefDir, nullDir), ('blank', nullDir),
                  (nullDir, prefDir), (prefDir, prefDir), ('blank', prefDir),
                  (nullDir, 'blank'), (prefDir, 'blank'), ('blank', 'blank')]
    yMax = 0
    plotCount = 0
    for row in range(3):
        for col in range(3):
            ax = fig.add_subplot(gs00[row,col])
            locDir = locDirList[plotCount]
            loc0Con = highContrast
            loc1Con = highContrast
            loc0Dir, loc1Dir = locDir
            loc0Title = loc0Dir
            loc1Title = loc1Dir
            if loc0Dir == 'blank':
                loc0Dir = zeroDir
                loc0Con = zeroContrast
            if loc1Dir == 'blank':
                loc1Dir = zeroDir
                loc1Con = zeroContrast
            histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                          (stimIndexDF['loc0 Contrast'] == loc0Con) &
                                          (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                          (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
            dirPlot = spikeHists[unitCount, histIndex, :] * 1000/stimIndexCount[histIndex]
            smoothPlot = gaussian_filter1d(dirPlot, 5)
            if max(smoothPlot) > yMax:
                yMax = max(smoothPlot)
            ax.plot(smoothPlot)
            ax.set_title(f'loc0: {loc0Title}, loc1: {loc1Title}', fontsize=5)
            ax.set_xticks([0,
                           histPrePostMS,
                           histPrePostMS+trueStimDurMS,
                           2*histPrePostMS+trueStimDurMS])
            ax.set_xticklabels([])
            ax.set_xlim([0, trueStimDurMS+(2*histPrePostMS+1)])
            ax.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS,
                       color='grey', alpha=0.1)
            ax.axhline(y=meanSpike[unitCount][48]*1000/trueStimDurMS,
                       linestyle='--', color='grey')
            if plotCount == 6:
                ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax.set_xlabel('Stimulus Duration (ms)', fontsize=7)
                ax.set_xticklabels([-histPrePostMS, 0, 0+trueStimDurMS,
                                    trueStimDurMS+histPrePostMS],
                                   fontsize=7)
            plotCount += 1

        axes = fig.get_axes()
        for ax in axes:
            ax.set_ylim([0, yMax*1.1])


    # Normalization Plot
    # ReIndexing to have pref direction in the middle Gauss Smooth
    nullIndex = np.where(dirArray == nullDir)[0][0]
    reIndex = (np.array([0, 1, 2, 3, 4, 5])+nullIndex) % 6
    tickLabels = np.array(['0', '60', '120', '180', '240', '300'])[reIndex]
    tickLabels = np.append(tickLabels, ['blank'])

    # smoothed data
    bSmoothReIndex = np.zeros((7, 7))
    tempMain = bSmooth[:6, :6][:, reIndex]
    tempMain = tempMain[:6, :6][reIndex, :]
    temp0Blank = bSmooth[:6, 6][reIndex]
    temp1Blank = bSmooth[6, :6][reIndex]
    bSmoothReIndex[:6, :6] = tempMain
    bSmoothReIndex[:6, 6] = temp0Blank
    bSmoothReIndex[6, :6] = temp1Blank
    bSmoothReIndex[6, 6] = bSmooth[6, 6]

    vMax = np.max(bSmoothReIndex)

    # raw data
    bReIndex = np.zeros((7, 7))
    tempMain = b[:6, :6][:, reIndex]
    tempMain = tempMain[:6, :6][reIndex, :]
    temp0Blank = b[:6, 6][reIndex]
    temp1Blank = b[6, :6][reIndex]
    bReIndex[:6, :6] = tempMain
    bReIndex[:6, 6] = temp0Blank
    bReIndex[6, :6] = temp1Blank
    bReIndex[6, 6] = b[6, 6]

    if np.max(bReIndex) > vMax:
        vMax = np.max(bReIndex)

    # fitted data
    stimMat = np.zeros((7, 7))
    stimMat[:6, :6] = np.arange(36).reshape(6, 6)
    stimMat[6, :6] = np.arange(36, 42)
    stimMat[:, 6] = np.arange(42, 49)

    stimMatReIndex = np.zeros((7, 7))
    tempMain = stimMat[:6, :6][:, reIndex]
    tempMain = np.squeeze(tempMain[:6, :6][reIndex, :])
    temp0Blank = stimMat[:6, 6][reIndex]
    temp1Blank = stimMat[6, :6][reIndex]
    stimMatReIndex[:6, :6] = tempMain
    stimMatReIndex[:6, 6] = temp0Blank
    stimMatReIndex[6, :6] = temp1Blank
    stimMatReIndex[6, 6] = stimMat[6, 6]

    y_pred = np.zeros(49)
    # fitting Gaussian EMS to paired stimuli
    fixedVals = fixedValsForGenericNorm(stimMatReIndex.reshape(49),
                                        stimIndexDict)

    y_pred = genericNormNoScalar(fixedVals, *unitNormFitEstimate[unitCount])
    y_pred = y_pred.reshape((7, 7))

    if np.max(y_pred) > vMax:
        vMax = np.max(y_pred)

    # residual (diff between real vs predicted resp)
    residual = abs(y_pred - bReIndex)
    if np.max(residual) > vMax:
        vMax = np.max(residual)

    # scatter plot of fit responses vs real responses
    respReal = bReIndex.reshape(49)
    respFit = y_pred.reshape(49)
    ax2 = fig.add_subplot(gs01[0, 1])
    ax2.scatter(respFit, respReal)
    ax2.set_ylabel('Real Responses (spikes/sec)', fontsize=7)
    ax2.set_xlabel('Fit Responses (spikes/sec)', fontsize=7)
    ax2.xaxis.set_label_position('top')
    ax2.set_ylim([-2, np.max(respReal)*1.10])
    ax2.set_xlim([-2, np.max(respReal)*1.10])
    line = lines.Line2D([0, 1], [0, 1], color='red')
    transform = ax2.transAxes
    line.set_transform(transform)
    ax2.add_line(line)
    # ax2 = plt.axis('equal')

    ax3 = fig.add_subplot(gs01[0, 0])
    ax3 = sns.heatmap(bReIndex, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
    ax3.set_xticks(np.arange(7)+0.5)
    ax3.set_title(f'Raw Data', y=-0.1, fontsize=7)
    ax3.set_xlabel('Location 0 Stimulus Direction', fontsize=7)
    ax3.xaxis.set_label_position('top')
    ax3.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax3.set_ylabel('Location 1 Stimulus Direction', fontsize=7)
    ax3.xaxis.set_ticks_position("top")
    ax3.set_yticks(np.arange(7)+0.5)
    ax3.set_yticklabels(tickLabels, rotation=0, fontsize=7)

    ax6 = fig.add_subplot(gs01[1, 0])
    ax6 = sns.heatmap(y_pred, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
    ax6.set_xticks(np.arange(7) + 0.5)
    ax6.set_title(f'Predicted', y=-0.1, fontsize=7)
    ax6.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax6.set_ylabel('Location 1 Stimulus Direction', fontsize=7)
    ax6.xaxis.set_ticks_position("top")
    ax6.set_yticks(np.arange(7) + 0.5)
    ax6.set_yticklabels(tickLabels, rotation=0, fontsize=7)

    ax7 = fig.add_subplot(gs01[1, 1])
    ax7 = sns.heatmap(residual, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
                      # cbar=True, cbar_kws={'orientation': 'horizontal'})
    ax7.set_xticks(np.arange(7) + 0.5)
    ax7.set_title(f'Residual (Predicted-Raw)', y=-0.1, fontsize=7)
    ax7.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax7.xaxis.set_ticks_position("top")
    ax7.set_yticks(np.arange(7) + 0.5)
    ax7.set_yticklabels(tickLabels, rotation=0, fontsize=7)


    # bar plot
    ax4 = fig.add_subplot(gs02[0, 0])
    unitDF = pnContrastPlotDF.loc[pnContrastPlotDF['unit'] == unit]

    #spikeCounts in spikes/sec
    unitDF['stimSpikes'] = unitDF['stimSpikes'] * 1000/trueStimDurMS

    offset = lambda p: transforms.ScaledTranslation(p/72., 0, plt.gcf().dpi_scale_trans)
    trans = plt.gca().transData
    offsetCount = 0

    colorList = ['b', 'g', 'r', 'c']
    # plot category by category
    yMax = 0
    for pos in [0.35, 0.65]:
        if pos == 0.35:
            for count, indx in enumerate(['pref+blank', 'null+blank',
                                          'blank+pref', 'blank+null']):
                mean = unitDF.loc[(unitDF['prefNullStr'] == indx)].mean()[3]
                sem = unitDF.loc[(unitDF['prefNullStr'] == indx)].sem()[3]
                ax4 = plt.scatter(pos, mean, transform=trans+offset(offsetCount),
                                  label=indx, color=colorList[count])
                ax4 = plt.errorbar(pos, mean, yerr=sem, fmt='o',
                                   transform=trans+offset(offsetCount),
                                   color=colorList[count])
                offsetCount += 5
                stimIndex = unitDF.loc[(unitDF['prefNullStr'] == indx), 'stimIndex'].iloc[0]

                # # using EMS generic for L0-L6
                # fixedVals = fixedValsForEMSGen(stimIndex, stimIndexDict)
                # if scaledLoc[unitCount] == 0:
                #     y_pred = emsNormGen1(fixedVals, *unitNormFitEstimate[unitCount])
                # else:
                #     y_pred = emsNormGen0(fixedVals, *unitNormFitEstimate[unitCount])

                # using gaussian EMS for single stim
                fixedVals = fixedValsForGenericNorm(stimIndex, stimIndexDict)
                y_pred = genericNormNoScalar(fixedVals, *unitNormFitEstimate[unitCount])
                ax4 = plt.scatter(pos, y_pred, transform=trans+offset(offsetCount),
                                  label=indx, color=lightenColor(colorList[count], 0.5))
                offsetCount += 5
                if mean > yMax:
                    yMax = mean
        else:
            offsetCount = 0
            for count, indx in enumerate(['pref+pref', 'pref+null',
                                          'null+pref', 'null+null']):
                mean = unitDF.loc[(unitDF['prefNullStr'] == indx)].mean()[3]
                sem = unitDF.loc[(unitDF['prefNullStr'] == indx)].sem()[3]
                ax4 = plt.scatter(pos,mean, transform=trans+offset(offsetCount),
                                  label=indx, color=colorList[count])
                ax4 = plt.errorbar(pos, mean, yerr=sem, fmt="o",
                                   transform=trans+offset(offsetCount),
                                   color=colorList[count])
                offsetCount += 5
                stimIndex = unitDF.loc[(unitDF['prefNullStr'] == indx), 'stimIndex'].iloc[0]
                fixedVals = fixedValsForGenericNorm(stimIndex, stimIndexDict)
                y_pred = genericNormNoScalar(fixedVals, *unitNormFitEstimate[unitCount])
                ax4 = plt.scatter(pos, y_pred, transform=trans+offset(offsetCount),
                                  label=indx, color=lightenColor(colorList[count], 0.5))
                offsetCount += 5
                if mean > yMax:
                    yMax = mean
    ax4 = plt.axhline(y=meanSpike[unitCount][48]*1000/trueStimDurMS, linestyle='--', color='grey')
    ax4 = plt.ylim([0, yMax*1.5])
    ax4 = plt.xticks([0.35, 0.65], ['Singular Stimulus', 'Dual Stimulus'])
    ax4 = plt.xlim(left=0.2, right=0.8)
    ax4 = plt.ylabel('Firing Rate spikes/sec')
    ax4 = plt.legend(loc='upper right', prop={'size': 6}, bbox_to_anchor=(1.25, 1.0))

    # norm fit parameters with EMS generic (l1-l6)
    ax5 = fig.add_subplot(gs03[0, 0])
    ax5.text(0.5,0.5, f'L0_0: {unitNormFitEstimate[unitCount][0]:.2f}\n\
    L0_60: {unitNormFitEstimate[unitCount][1]:.2f}\n\
    L0_120: {unitNormFitEstimate[unitCount][2]:.2f}\n\
    L0_180: {unitNormFitEstimate[unitCount][3]:.2f}\n\
    L0_240: {unitNormFitEstimate[unitCount][4]:.2f}\n\
    L0_300: {unitNormFitEstimate[unitCount][5]:.2f}\n\
    L1_0: {unitNormFitEstimate[unitCount][6]:.2f}\n\
    L1_60: {unitNormFitEstimate[unitCount][7]:.2f}\n\
    L1_120: {unitNormFitEstimate[unitCount][8]:.2f}\n\
    L1_180: {unitNormFitEstimate[unitCount][9]:.2f}\n\
    L1_240: {unitNormFitEstimate[unitCount][10]:.2f}\n\
    L1_300: {unitNormFitEstimate[unitCount][11]:.2f}\n\
    alpha loc0: {unitNormFitEstimate[unitCount][12]:.2f}\n\
    alpha loc1: {unitNormFitEstimate[unitCount][13]:.2f}\n\
    fit R2: {unitPairedNormR2[unitCount]:.2f}', size=10, ha='center', transform=ax5.transAxes)
    ax5.axis('off')

    ax8 = fig.add_subplot(gs03[1, 0], polar='True')
    theta = np.radians(np.arange(0, 420, 60))
    sponTheta = np.radians(np.arange(0, 360, 360/100))
    sponTheta = np.append(sponTheta, sponTheta[0])

    # loc 0
    r0 = (np.append(meanSpike[unitCount][36:42], meanSpike[unitCount][36])) \
          * 1000/trueStimDurMS
    er0 = (np.append(spikeCountSEM[unitCount][36:42], spikeCountSEM[unitCount][36])) \
          * 1000/trueStimDurMS
    ax8.plot(theta, r0, markersize=2, color='green', label='location 0')
    ax8.errorbar(theta, r0, yerr=er0, fmt='o', ecolor='green',
                 color='green', markersize=2)

    # loc 1
    r1 = (np.append(meanSpike[unitCount][42:48], meanSpike[unitCount][42])) \
        * 1000/trueStimDurMS
    er1 = (np.append(spikeCountSEM[unitCount][42:48], spikeCountSEM[unitCount][42])) \
        * 1000/trueStimDurMS
    ax8.plot(theta, r1, markersize=2, color='red', label='location 1')
    ax8.errorbar(theta, r1, yerr=er1, fmt='x', ecolor='red',
                 color='red', markersize=2)
    ax8.set_theta_zero_location('W')
    ax8.set_title('Direction tuning at both locations')
    spon = np.array([meanSpike[unitCount][48] * 1000/trueStimDurMS] * len(sponTheta))
    ax8.plot(sponTheta, spon, linestyle='--', color='blue',
             label='spontaneous rate')
    ax8.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0,
               prop={'size': 6})

    plt.tight_layout()

    plt.savefig(f'{unit}EMS2PartsuperPlot.pdf')
    plt.close('all')


# population tuning curve based off of MTNC stimuli aligned to preferred
# stimulus
popTunCurve = np.array(popTunCurve)
meanPopTunResp = np.mean(popTunCurve[:, :-2], axis=0)
semPopTunResp = stats.sem(popTunCurve[:, :-2], axis=0)
meanBaseResp = np.mean(popTunCurve[:, -1], axis=0)
semBaseResp = stats.sem(popTunCurve[:, -1], axis=0)
xDeg = np.arange(-180, 240, 60)
params = gaussFit(xDeg+180, meanPopTunResp)
params[2] = params[2] - 180
xFull = np.linspace(-180, 180, 1000)
yFull = gauss(xFull, *params)

# figure
fig = plt.figure()
fig.set_size_inches(15, 6)
gs0 = gridspec.GridSpec(1, 3)
gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0])
gs01 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[2])

ax1 = fig.add_subplot(gs00[0, 0])
ax1.plot(xDeg, meanPopTunResp, label='raw data', color='black', alpha=0.1)
ax1.errorbar(xDeg, meanPopTunResp, yerr=semPopTunResp, color='black', alpha=0.1)
ax1.plot(xFull, yFull, linestyle='--', color='green', label='fitted Gaussian')
ax1.axhline(y=meanBaseResp, linestyle='--', color='black')
y = np.repeat(meanBaseResp, 7)
ax1.fill_between(xDeg, (y+semBaseResp), (y-semBaseResp), color='b', alpha=.1)
ax1.set_xticks(np.arange(-180, 240, 60))
ax1.set_xticklabels(np.arange(-180, 240, 60))
ax1.set_ylim(bottom=0)
ax1.set_ylabel('Normalized Firing Rate')
ax1.set_xlabel('Stimulus Direction Relative to Preferred')
ax1.set_title('Population Average Tuning Curve Aligned to Preferred Direction',
              fontsize=7)
plt.legend()

# polar aligned to preferred direction
gaussFitR2 = []
filteredGaussFitR2 = []
ax2 = fig.add_subplot(gs01[0, 0], polar='True')
for i in popTunCurve:
    params = gaussFit(xDeg+180, i[:-2])
    yPred = gauss(xDeg+180, *params)
    r2 = r2_score(i[:-2], yPred)
    gaussFitR2.append(r2)
    if r2 > 0.90:
        filteredGaussFitR2.append(r2)
        params[2] = params[2] - 180
        xFull = np.linspace(-180, 180, 7)
        yFull = gauss(xFull, *params)
        theta = np.radians(np.linspace(-180, 180, 7))
        ax2.plot(theta, yFull, alpha=0.3)
        ax2.set_theta_zero_location('W')
        ax2.set_title('Individual Gaussian Fits Aligned to Preferred Direction',
                      fontsize=7)
ax2.plot(theta, meanPopTunResp, alpha=1, color='black')
ax2.errorbar(theta, meanPopTunResp, yerr=semPopTunResp, fmt='o',
             ecolor='black', color='black', markersize=0.5, alpha=1)
sponTheta = np.radians(np.linspace(0, 360, 360))
spon = np.array([meanBaseResp] * len(sponTheta))
ax2.plot(sponTheta, spon, alpha=1, color='black', linestyle='--')

# polar plot of direction frequency distribution
ax3 = fig.add_subplot(gs02[0, 0], polar='True')
circularHist(ax3, popTunCurve[:, -2], bins=24, density=False)
ax3.set_theta_zero_location('W')
ax3.set_title('Distribution of Preferred Directions Sampled',
              fontsize=7)
ax3.grid(alpha=0.4)

plt.tight_layout()
plt.show()

# plotting alpha as a function of FR
highFRAlpha = np.array([item for sublist in highFRAlpha for item in sublist])
medFRAlpha = np.array([item for sublist in medFRAlpha for item in sublist])
lowFRAlpha = np.array([item for sublist in lowFRAlpha for item in sublist])
highFRParams = np.array(highFRParams)
medFRParams = np.array(medFRParams)
lowFRParams = np.array(lowFRParams)
unitIdentifier = np.array(unitIdentifier)

highFRParamsUnitID = np.hstack((highFRParams, unitIdentifier[:, np.newaxis]))
medFRParamsUnitID = np.hstack((medFRParams, unitIdentifier[:, np.newaxis]))
lowFRParamsUnitID = np.hstack((lowFRParams, unitIdentifier[:, np.newaxis]))

plt.hist(highFRAlpha)
plt.xlabel('alpha value')
plt.ylabel('frequency')
plt.title(f'High FR median alpha: {np.median(highFRAlpha)}')
plt.show()

plt.hist(medFRAlpha)
plt.xlabel('alpha value')
plt.ylabel('frequency')
plt.title(f'Med FR median alpha: {np.median(medFRAlpha)}')
plt.show()

plt.hist(lowFRAlpha)
plt.xlabel('alpha value')
plt.ylabel('frequency')
plt.title(f'Low FR median alpha: {np.median(lowFRAlpha)}')
plt.show()

## backup
# n1 selectivity and suppression index
loc0Resp = unitNormFitEstimate[n1][np.where(condArr == loc0Dir)[0][0]]
loc1Resp = unitNormFitEstimate[n1][np.where(condArr == loc1Dir)[0][0] + 2]
n1Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
n1NonPrefSupp = (1 - unitNormFitEstimate[n1][4]) / (
        1 + unitNormFitEstimate[n1][4])
if n1Selectivity >= 0:
    n1NonPrefSupp = 0.5 * (-n1NonPrefSupp + 1)
else:
    n1NonPrefSupp = 0.5 * (n1NonPrefSupp + 1)

# n2 selectivity and suppression index
loc0Resp = unitNormFitEstimate[n2][np.where(condArr == loc0Dir)[0][0]]
loc1Resp = unitNormFitEstimate[n2][np.where(condArr == loc1Dir)[0][0] + 2]
n2Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
n2NonPrefSupp = (1 - unitNormFitEstimate[n2][4]) / (
        1 + unitNormFitEstimate[n2][4])
if n2Selectivity >= 0:
    n2NonPrefSupp = 0.5 * (-n2NonPrefSupp + 1)
else:
    n2NonPrefSupp = 0.5 * (n2NonPrefSupp + 1)

# # indices and correlations for single Gabor stimuli
# singleStimIndex = np.int_(np.concatenate((stimMatCond[2, :2],
#                                           stimMatCond[:2, 2]),
#                                          axis=0))
# for count, i in enumerate(singleStimIndex):
#     # correlation for the single Gabor b/w 2 units excluding trials where
#     # spike counts exceeded 3 SD from mean
#     skipTrials = []
#     n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
#     n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
#     n1Zscore = stats.zscore(n1SpikeMat)
#     n2Zscore = stats.zscore(n2SpikeMat)
#     n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
#     n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
#     for x in n1SkipTrials:
#         skipTrials.append(x)
#     for x in n2SkipTrials:
#         skipTrials.append(x)
#     # print(skipTrials)
#     goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
#     singleStimCorr = stats.pearsonr(spikeCountMat[n1, goodTrials, i],
#                                   spikeCountMat[n2, goodTrials, i])
#
#     # extract direction and contrast of the gabor
#     loc0Dir = stimIndexDict[i][0]['direction']
#     loc1Dir = stimIndexDict[i][1]['direction']
#     loc0Con = stimIndexDict[i][0]['contrast']
#     loc1Con = stimIndexDict[i][1]['contrast']
#
#     # n1 selectivity and suppression index based of Norm params
#     loc0Dir0Resp = unitNormFitEstimate[n1][0]
#     loc0Dir1Resp = unitPairedNormFit[n1][1]
#     loc1Dir0Resp = unitNormFitEstimate[n1][2]
#     loc1Dir1Resp = unitNormFitEstimate[n1][3]
#     minRespN1 = np.where(unitNormFitEstimate[n1][:4] ==
#                          min(unitNormFitEstimate[n1][:4]))[0]
#     maxRespN1 = np.where(unitNormFitEstimate[n1][:4] ==
#                          max(unitNormFitEstimate[n1][:4]))[0]
#     if maxRespN1 < 2:
#         n1PrefResp = 0
#     else:
#         n1PrefResp = 1
#     if minRespN1 < 2:
#         n1NonPrefSupp = unitNormFitEstimate[n1][4] / (
#             unitNormFitEstimate[n1][4] + unitNormFitEstimate[n1][5])
#     else:
#         n1NonPrefSupp = unitNormFitEstimate[n1][5] / (
#             unitNormFitEstimate[n1][4] + unitNormFitEstimate[n1][5])
#     if loc0Con == 1:
#         if loc0Dir0Resp > loc0Dir1Resp:
#             n1Selectivity = (loc0Dir0Resp - loc0Dir1Resp) / (
#                              loc0Dir0Resp + loc0Dir1Resp)
#         else:
#             n1Selectivity = (loc0Dir1Resp - loc0Dir0Resp) / (
#                              loc0Dir0Resp + loc0Dir1Resp)
#     else:
#         if loc1Dir0Resp > loc1Dir1Resp:
#             n1Selectivity = (loc1Dir0Resp - loc1Dir1Resp) / (
#                              loc1Dir0Resp + loc1Dir1Resp)
#         else:
#             n1Selectivity = (loc1Dir1Resp - loc1Dir0Resp) / (
#                              loc1Dir0Resp + loc1Dir1Resp)
#
#     # n2 selectivity and suppression index based of Norm params
#     loc0Dir0Resp = unitNormFitEstimate[n2][0]
#     loc0Dir1Resp = unitPairedNormFit[n2][1]
#     loc1Dir0Resp = unitNormFitEstimate[n2][2]
#     loc1Dir1Resp = unitNormFitEstimate[n2][3]
#     minRespN2 = np.where(unitNormFitEstimate[n2][:4] ==
#                          min(unitNormFitEstimate[n2][:4]))[0]
#     maxRespN2 = np.where(unitNormFitEstimate[n2][:4] ==
#                          max(unitNormFitEstimate[n2][:4]))[0]
#     if maxRespN2 < 2:
#         n2PrefResp = 0
#     else:
#         n2PrefResp = 1
#     if minRespN2 < 2:
#         n2NonPrefSupp = unitNormFitEstimate[n2][4] / (
#                 unitNormFitEstimate[n2][4] + unitNormFitEstimate[n2][5])
#     else:
#         n2NonPrefSupp = unitNormFitEstimate[n2][5] / (
#                 unitNormFitEstimate[n2][4] + unitNormFitEstimate[n2][5])
#     if loc0Con == 1:
#         if loc0Dir0Resp > loc0Dir1Resp:
#             n2Selectivity = (loc0Dir0Resp - loc0Dir1Resp) / (
#                     loc0Dir0Resp + loc0Dir1Resp)
#         else:
#             n2Selectivity = (loc0Dir1Resp - loc0Dir0Resp) / (
#                     loc0Dir0Resp + loc0Dir1Resp)
#     else:
#         if loc1Dir0Resp > loc1Dir1Resp:
#             n2Selectivity = (loc1Dir0Resp - loc1Dir1Resp) / (
#                     loc1Dir0Resp + loc1Dir1Resp)
#         else:
#             n2Selectivity = (loc1Dir1Resp - loc1Dir0Resp) / (
#                     loc1Dir0Resp + loc1Dir1Resp)
#
#     # pair selectivity, suppression index for single gabors
#     pairSelectivity = np.sqrt(n1Selectivity * n2Selectivity)
#     pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)
#     if n1PrefResp != n2PrefResp:
#         pairSelectivity = -pairSelectivity
#
#     pairSingleCorr.append(singleStimCorr[0])
#     pairSingleSelectivityIndex.append(pairSelectivity)
#     pairSingleNonPrefSuppIndex.append(pairSuppression)

# # indices and correlations for paired Gabor stimuli
# for i in range(36):
#
#     # correlation for that Gabor pair b/w 2 units excluding trials where
#     # spike counts exceeded 3 SD from mean
#     skipTrials = []
#     n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
#     n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
#     n1Zscore = stats.zscore(n1SpikeMat)
#     n2Zscore = stats.zscore(n2SpikeMat)
#     n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
#     n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
#     for x in n1SkipTrials:
#         skipTrials.append(x)
#     for x in n2SkipTrials:
#         skipTrials.append(x)
#     # print(skipTrials)
#     goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
#     pairStimCorr = stats.pearsonr(spikeCountMat[n1, goodTrials, i],
#                                   spikeCountMat[n2, goodTrials, i])
#     # print(pairStimCorr)
#     #
#     # pairStimCorr = stats.pearsonr(spikeCountMat[n1, :blocksDone, i],
#     #                               spikeCountMat[n2, :blocksDone, i])
#     # print(pairStimCorr)
#
#     # extract directions of the gabor pair
#     loc0Dir = stimIndexDict[i][0]['direction']
#     loc1Dir = stimIndexDict[i][1]['direction']
#
#     # # n1 selectivity and suppression based off of average spike counts
#     # b = meanSpikeReshaped[n1][0] * 1000/trueStimDurMS
#     # pairedMat = b[:36]
#     # loc0Resp = b[loc0IndexArray[np.where(dirArray == loc0Dir)[0][0]]]
#     # loc1Resp = b[loc1IndexArray[np.where(dirArray == loc1Dir)[0][0]]]
#     # if loc0Resp > loc1Resp:
#     #     # n1Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
#     #     n1Selectivity = (loc0Resp - loc1Resp) / loc0Resp
#     #     n1NonPrefSupp = abs(loc0Resp - pairedMat[i]) / (loc0Resp + pairedMat[i])
#     #     n1PrefResp = 0
#     # else:
#     #     # n1Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
#     #     n1Selectivity = (loc1Resp - loc0Resp) / loc1Resp
#     #     n1NonPrefSupp = abs(loc1Resp - pairedMat[i]) / (loc1Resp + pairedMat[i])
#     #     n1PrefResp = 1
#     #
#     # # n2 selectivity and suppression bvased off of average spike counts
#     # b = meanSpikeReshaped[n2][0] * 1000/trueStimDurMS
#     # pairedMat = b[:36]
#     # loc0Resp = b[loc0IndexArray[np.where(dirArray == loc0Dir)[0][0]]]
#     # loc1Resp = b[loc1IndexArray[np.where(dirArray == loc1Dir)[0][0]]]
#     # if loc0Resp > loc1Resp:
#     #     # n2Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
#     #     n2Selectivity = (loc0Resp - loc1Resp) / loc0Resp
#     #     n2NonPrefSupp = abs(loc0Resp - pairedMat[i]) / (loc0Resp + pairedMat[i])
#     #     n2PrefResp = 0
#     # else:
#     #     # n2Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
#     #     n2Selectivity = (loc1Resp - loc0Resp) / loc1Resp
#     #     n2NonPrefSupp = abs(loc1Resp - pairedMat[i]) / (loc1Resp + pairedMat[i])
#     #     n2PrefResp = 1
#
#     # n1 selectivity and suppression index based of Normalization params
#     loc0Resp = unitNormFitEstimate[n1][np.where(dirArray == loc0Dir)[0][0]]
#     loc1Resp = unitNormFitEstimate[n1][np.where(dirArray == loc1Dir)[0][0] + 6]
#     if loc0Resp > loc1Resp:
#         n1Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
#         n1NonPrefSupp = unitNormFitEstimate[n1][13] / (
#                         unitNormFitEstimate[n1][13] + unitNormFitEstimate[n1][12])
#         n1PrefResp = 0
#         if abs(loc0Resp - n1PrefDir) > 180:
#             n1Diff = 360 - abs(loc0Resp - n1PrefDir)
#         else:
#             n1Diff = abs(loc0Resp - n1PrefDir)
#     else:
#         n1Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
#         n1NonPrefSupp = unitNormFitEstimate[n1][12] / (
#                         unitNormFitEstimate[n1][13] + unitNormFitEstimate[n1][12])
#         n1PrefResp = 1
#         if abs(loc1Resp - n1PrefDir) > 180:
#             n1Diff = 360 - abs(loc1Resp - n1PrefDir)
#         else:
#             n1Diff = abs(loc1Resp - n1PrefDir)
#
#     # n2 selectivity and suppression index
#     loc0Resp = unitNormFitEstimate[n2][np.where(dirArray == loc0Dir)[0][0]]
#     loc1Resp = unitNormFitEstimate[n2][np.where(dirArray == loc1Dir)[0][0] + 6]
#     if loc0Resp > loc1Resp:
#         n2Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
#         n2NonPrefSupp = unitNormFitEstimate[n2][13] / (
#                         unitNormFitEstimate[n2][13] + unitNormFitEstimate[n2][12])
#         n2PrefResp = 0
#         if abs(loc0Resp - n2PrefDir) > 180:
#             n2Diff = 360 - abs(loc0Resp - n2PrefDir)
#         else:
#             n2Diff = abs(loc0Resp - n2PrefDir)
#     else:
#         n2Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
#         n2NonPrefSupp = unitNormFitEstimate[n2][12] / (
#                         unitNormFitEstimate[n2][13] + unitNormFitEstimate[n2][12])
#         n2PrefResp = 1
#         if abs(loc1Resp - n2PrefDir) > 180:
#             n2Diff = 360 - abs(loc1Resp - n2PrefDir)
#         else:
#             n2Diff = abs(loc1Resp - n2PrefDir)
#
#     # n1 NMI
#     b = meanSpikeReshaped[n1].reshape(7, 7) * 1000/trueStimDurMS
#     pairedMat = b[:6, :6].reshape(36)
#     loc0Single = b[6, :6]
#     loc1Single = b[:6, 6]
#     loc0SingleResp = loc0Single[np.where(dirArray == loc0Dir)[0]][0]
#     loc1SingleResp = loc1Single[np.where(dirArray == loc1Dir)[0]][0]
#     n1NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[i]) / (
#           (loc0SingleResp + loc1SingleResp) + pairedMat[i])
#     n1SimpleNMI = pairedMat[i] / (loc0SingleResp + loc1SingleResp)
#
#     # n2 NMI
#     b = meanSpikeReshaped[n2].reshape(7, 7) * 1000/trueStimDurMS
#     pairedMat = b[:6, :6].reshape(36)
#     loc0Single = b[6, :6]
#     loc1Single = b[:6, 6]
#     loc0SingleResp = loc0Single[np.where(dirArray == loc0Dir)[0]][0]
#     loc1SingleResp = loc1Single[np.where(dirArray == loc1Dir)[0]][0]
#     n2NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[i]) / (
#           (loc0SingleResp + loc1SingleResp) + pairedMat[i])
#     n2simpleNMI = pairedMat[i] / (loc0SingleResp + loc1SingleResp)
#
#     # pair selectivity, suppression, and NMI index
#     pairSelectivity = np.sqrt(n1Selectivity * n2Selectivity)
#     pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)
#     if n1PrefResp != n2PrefResp:
#         pairSelectivity = -pairSelectivity
#     pairNMI = (n1NMI + n2NMI) / 2
#     pairSimpleNMI = np.sqrt(n1SimpleNMI * n2simpleNMI)
#
#     pairPairedCorr.append(pairStimCorr[0])
#     pairSelectivityIndex.append(pairSelectivity)
#     pairNonPrefSuppIndex.append(pairSuppression)
#     pairNMIIndex.append(pairNMI)
#     pairSimpleNMIIndex.append(pairSimpleNMI)
#
#     # appending pair Selectivity and pair NMI to get average
#     # for the single Gabor figure
#     neuronPairAvgSelectivity.append(pairSelectivity)
#     if i in upperTriangle:
#         neuronPairAvgNMIloc0.append(pairNMI)
#     else:
#         neuronPairAvgNMIloc1.append(pairNMI)

# # indices and correlations for single Gabor b/w pair of neuron
# for j in range(36, 48):
#     pairSingCorr = stats.pearsonr(spikeCountMat[n1, :blocksDone, j],
#                                   spikeCountMat[n2, :blocksDone, j])
#
#     # extract location and direction of tested gabor
#     loc0Dir = stimIndexDict[j][0]['direction']
#     loc0Con = stimIndexDict[j][0]['contrast']
#     loc1Dir = stimIndexDict[j][1]['direction']
#     loc1Con = stimIndexDict[j][1]['contrast']
#     if loc0Con == 1:
#         avgNMI = np.mean(neuronPairAvgNMIloc1)
#     else:
#         avgNMI = np.mean(neuronPairAvgNMIloc0)
#
#     pairSingleGaborNMI.append(avgNMI)
#     pairSingleGaborSelectivity.append(np.mean(neuronPairAvgSelectivity))
#     pairSingleCorr.append(pairSingCorr[0])

# # n1 selectivity and suppression index based of Normalization params
# loc0Resp = unitNormFitEstimate[n1][np.where(condArr == loc0Dir)[0][0]]
# loc1Resp = unitNormFitEstimate[n1][np.where(condArr == loc1Dir)[0][0] + 2]
# if loc0Resp > loc1Resp:
#     n1Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
#     n1NonPrefSupp = unitNormFitEstimate[n1][5] / (
#                     unitNormFitEstimate[n1][5] + unitNormFitEstimate[n1][4])
#     n1PrefResp = 0
# else:
#     n1Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
#     n1NonPrefSupp = unitNormFitEstimate[n1][4] / (
#                     unitNormFitEstimate[n1][4] + unitNormFitEstimate[n1][5])
#     n1PrefResp = 1
#
# # n2 selectivity and suppression index
# loc0Resp = unitNormFitEstimate[n2][np.where(condArr == loc0Dir)[0][0]]
# loc1Resp = unitNormFitEstimate[n2][np.where(condArr == loc1Dir)[0][0] + 2]
# if loc0Resp > loc1Resp:
#     n2Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
#     n2NonPrefSupp = unitNormFitEstimate[n2][5] / (
#                     unitNormFitEstimate[n2][5] + unitNormFitEstimate[n2][4])
#     n2PrefResp = 0
# else:
#     n2Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
#     n2NonPrefSupp = unitNormFitEstimate[n2][4] / (
#                     unitNormFitEstimate[n2][4] + unitNormFitEstimate[n2][5])
#     n2PrefResp = 1
#
# # n1 NMI
# b = meanSpikeReshaped[n1].reshape(7, 7) * 1000/trueStimDurMS
# bCondensed = b[sli, :]
# bCondensed = bCondensed[:, sli]
# pairedMat = bCondensed[:2, :2].reshape(4)
# loc0Single = bCondensed[2, :2]
# loc1Single = bCondensed[:2, 2]
# loc0SingleResp = loc0Single[np.where(condArr == loc0Dir)[0][0]]
# loc1SingleResp = loc1Single[np.where(condArr == loc1Dir)[0][0]]
# n1NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[count]) / (
#       (loc0SingleResp + loc1SingleResp) + pairedMat[count])
# n1SimpleNMI = pairedMat[count] / (loc0SingleResp + loc1SingleResp)
#
# # n2 NMI
# b = meanSpikeReshaped[n2].reshape(7, 7) * 1000/trueStimDurMS
# bCondensed = b[sli, :]
# bCondensed = bCondensed[:, sli]
# pairedMat = bCondensed[:2, :2].reshape(4)
# loc0Single = bCondensed[2, :2]
# loc1Single = bCondensed[:2, 2]
# loc0SingleResp = loc0Single[np.where(condArr == loc0Dir)[0][0]]
# loc1SingleResp = loc1Single[np.where(condArr == loc1Dir)[0][0]]
# n2NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[count]) / (
#       (loc0SingleResp + loc1SingleResp) + pairedMat[count])
# n2simpleNMI = pairedMat[count] / (loc0SingleResp + loc1SingleResp)
#
# # pair selectivity, suppression, and NMI index
# pairSelectivity = np.sqrt(n1Selectivity * n2Selectivity)
# pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)
# if n1PrefResp != n2PrefResp:
#     pairSelectivity = -pairSelectivity
# pairNMI = (n1NMI + n2NMI) / 2
# pairSimpleNMI = np.sqrt(n1SimpleNMI * n2simpleNMI)
#
# pairPairedCorr.append(pairStimCorr[0])
# pairSelectivityIndex.append(pairSelectivity)
# pairNonPrefSuppIndex.append(pairSuppression)
# pairNMIIndex.append(pairNMI)
# pairSimpleNMIIndex.append(pairSimpleNMI)
#
# if i in stimIndex12Cond:
#     pairSelFromCond.append(pairSelectivity)
#     pairNPSuppFromCond.append(pairSuppression)
#
# # indices and correlations for single Gabor stimuli
# singleStimIndex = np.int_(np.concatenate((stimMatCond[2, :2],
#                                           stimMatCond[:2, 2]),
#                                          axis=0))
#
# # loc 0 single Gabor corr
# stimIndex = np.int_(stimMatCond[2, :][np.where(condArr == loc0Dir)[0][0]])
# # correlation for that Gabor pair b/w 2 units excluding trials where
# # spike counts exceeded 3 SD from mean
# skipTrials = []
# n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
# n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
# n1Zscore = stats.zscore(n1SpikeMat)
# n2Zscore = stats.zscore(n2SpikeMat)
# n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
# n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
# for x in n1SkipTrials:
#     skipTrials.append(x)
# for x in n2SkipTrials:
#     skipTrials.append(x)
# goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
# singleStimLoc0Corr = stats.pearsonr(spikeCountMat[n1, goodTrials, stimIndex],
#                                     spikeCountMat[n2, goodTrials, stimIndex])
# pairSingleCorr.append(singleStimLoc0Corr[0])
# pairSingleSelectivityIndex.append(pairSelectivity)
# pairSingleNonPrefSuppIndex.append(pairSuppression)
#
# # loc 1 single Gabor corr
# stimIndex = np.int_(stimMatCond[:, 2][np.where(condArr == loc1Dir)[0][0]])
# # correlation for that Gabor pair b/w 2 units excluding trials where
# # spike counts exceeded 3 SD from mean
# skipTrials = []
# n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
# n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
# n1Zscore = stats.zscore(n1SpikeMat)
# n2Zscore = stats.zscore(n2SpikeMat)
# n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
# n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
# for x in n1SkipTrials:
#     skipTrials.append(x)
# for x in n2SkipTrials:
#     skipTrials.append(x)
# goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
# singleStimLoc1Corr = stats.pearsonr(spikeCountMat[n1, goodTrials, stimIndex],
#                                     spikeCountMat[n2, goodTrials, stimIndex])
# pairSingleCorr.append(singleStimLoc1Corr[0])
# pairSingleSelectivityIndex.append(pairSelectivity)
# pairSingleNonPrefSuppIndex.append(pairSuppression)