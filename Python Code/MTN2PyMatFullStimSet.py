# imports
import matplotlib.pyplot as plt
from usefulFns import *
from itertools import combinations
from scipy import stats
import time
import numpy as np
import pandas as pd

# good sessions: 221110, 221115, 221117, 221124, 221128, 221208, 221229, 230123
# okay sessions: 221010, 221013, 221108, 221206
# bad sessions: 230124

# for loop to run through all files
fileList = ['221010', '221013', '221108', '221110', '221115', '221117',
            '221124', '221128', '221206', '221208', '221229', '230123',
            '230126']
t0 = time.time()
totUnits = []
totR2 = []
totCombs = []
totFilterUnits = []
totFilterR2 = []
totFilterCombs = []
unitIdentifier = []
popTunCurve = []
pairSel12FromFull = []
pairNPSupp12FromFull = []
stimIndex12Cond1 = [0, 3, 18, 21]
stimIndex12Cond2 = [7, 10, 25, 28]
stimIndex12Cond3 = [14, 17, 32, 35]
loc0to1RespRatio = []
alphaLoc0 = []
alphaLoc1 = []
loc0PrefNormalized = []
loc1PrefNormalized = []
loc0NullNMI = []
loc1NullNMI = []
electrodeArr = np.array(np.arange(0, 32)).reshape(16, 2)

for fileIterator in fileList:

    # Load relevant file here with pyMat reader
    monkeyName = 'Meetz'
    seshDate = fileIterator
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
        currTrial['stimDesc']['data'] = [{k:v[i] for k,v in currTrial['stimDesc']['data'].items()}
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
    sponSpikeCountLong = []
    onLatency = 50/1000  # time in MS for counting window latency after stim on
    offLatency = 150/1000  # time in MS for counting window latency after stim off
    histPrePostMS = 100  # 100ms window pre/post stimulus on/off
    sponWindowMS = 50  # 50ms window before stimulus onset
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
                            # added 50ms onset latency for spike counts (100 for offset)
                            unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                            stimSpikes = np.where((unitTimeStamps >= (stimOnTimeS + onLatency)) &
                                         (unitTimeStamps <= (stimOffTimeS + offLatency)))[0]
                            spikeCountMat[unitCount][stCount][stimIndex] \
                                = len(stimSpikes)
                            spikeCountLong.append([unit,
                                                   stimIndex,
                                                   stimIndexCount[stimIndex],
                                                   len(stimSpikes)])

                            # Spontaneous Spikes
                            sponSpikes = np.where((unitTimeStamps >= (stimOnTimeS-(sponWindowMS/1000)))
                                                  & (unitTimeStamps <= stimOnTimeS))[0]
                            sponSpikeCountLong.append([unit, len(sponSpikes)])

                            # PSTHs
                            stimOnPreSNEV = stimOnTimeS - (histPrePostMS/1000)
                            stimOffPostSNEV = stimOffTimeS + (histPrePostMS/1000)
                            histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                                            & (unitTimeStamps <= stimOffPostSNEV)
                                                             )] - stimOnPreSNEV
                            histStimSpikes = np.int32(histStimSpikes*1000)
                            # spikeHists[unitCount, stimIndex, histStimSpikes] += 1
                            spikeHists[unitCount][stimIndex][histStimSpikes] += 1

    # mean, SEM, and reshaping of spikeCount matrices
    # create pandas dataframe of spikeCount with corresponding unit, stimIndex
    spikeCountDF = pd.DataFrame(spikeCountLong, columns=['unit', 'stimIndex',
                                                         'stimCount', 'stimSpikes'])
    sponSpikeCountDF = pd.DataFrame(sponSpikeCountLong, columns=['unit', 'sponSpikes'])
    sponSpikesMean = np.zeros(len(units))  # spikes in 50ms window
    sponSpikesSEM = np.zeros(len(units))
    for unitCount, unit in enumerate(units):
        sponSpikesMean[unitCount] = sponSpikeCountDF.loc[sponSpikeCountDF['unit'] == unit].mean()[1]
        sponSpikesSEM[unitCount] = sponSpikeCountDF.loc[sponSpikeCountDF['unit'] == unit].sem()[1]
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

    # Unit's preferred direction based off of
    # gaussian fit of MTNC stimuli:
    angleMat = np.arange(180, 900, 60)
    unitGaussMean = np.zeros(len(units))
    unitGaussSig = np.zeros(len(units))
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

    # SCIPY CURVEFIT/Minimize
    # scipy curvefit Normalization parameters will also generate population
    # tuning average aligned to preferred direction
    unitPairedNormFit = []
    unitPairedNormR2 = np.zeros(len(units))
    unitNormFitEstimate = [[] for i in range(len(units))]
    for unitCount, unit in enumerate(units):
        b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS

        # fixed (independent) variables - matrix of corresponding stim Indexes
        stimMat = np.zeros((7, 7))
        stimMat[:6, :6] = np.arange(36).reshape(6, 6)
        stimMat[6, :6] = np.arange(36, 42)
        stimMat[:, 6] = np.arange(42, 49)

        # Generic Normalization (L1+L2)/(al1+al2+sig) w.o scalar
        # fits L0-L6 for both locations separately loc0 and loc1 will
        # have different L0-L6 values
        guess0 = np.concatenate((b[6, :6], b[:6, 6], [0.4]), axis=0)
        resp = b.reshape(49)[:-1]
        fixedVals = fixedValsForGenericNorm(stimMat.reshape(49)[:-1], stimIndexDict)
        # pOpt, pCov = curve_fit(genericNormNoScalar, fixedVals, resp.squeeze(), p0=guess0,
        #                        bounds=((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        #                                (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 6)))
        pOpt, pCov = curve_fit(genericNormNoScalar, fixedVals, resp.squeeze(), p0=guess0,
                               maxfev=10000000)
        y_pred = genericNormNoScalar(fixedVals, *pOpt)
        r2 = r2_score(resp.squeeze(), y_pred)
        print(unit, r2)

        # # fit using scipy.optimize.minimize
        # p0 = np.concatenate((b[6, :6], b[:6, 6], [0.2, 0.2]), axis=0)
        # bnds = ((0, None), (0, None), (0, None), (0, None), (0, None),
        #         (0, None), (0, None), (0, None), (0, None), (0, None),
        #         (0, None), (0, None), (0, 5), (0, 5))
        # resp = b.reshape(49)[:-1]
        # fixedVals = fixedValsForGenericNorm(stimMat.reshape(49)[:-1], stimIndexDict)
        # res = scipy.optimize.minimize(driverFunc, p0, args=(fixedVals, resp),
        #                               method='Nelder-Mead', bounds=bnds,
        #                               options={'maxfev': 1 * (10 ** 20),
        #                                        'fatol': 0.11,
        #                                        'maxiter': 10000,
        #                                        'xatol': 0.0001})
        # y_pred = genericNormNoScalar(fixedVals, *res.x)
        # r2 = r2_score(resp.squeeze(), y_pred)
        # print(unit, r2)

        # Append fit parameters for full matrix
        totR2.append(r2)

        # bram trick to keep things non negative
        pOpt = pOpt ** 2

        # alphaLoc0.append(pOpt[12])
        # alphaLoc1.append(pOpt[13])
        # alphaLoc0.append(res.x[12])
        # alphaLoc1.append(res.x[13])
        unitPairedNormR2[unitCount] = r2
        unitPairedNormFit.append(pOpt)
        unitNormFitEstimate[unitCount] = pOpt
        # unitPairedNormFit.append(res.x)
        # unitNormFitEstimate[unitCount] = res.x

        # if r2 > 0.75:
        #     filterUnits.append(unit)
        #     totFilterR2.append(r2)

        # reindex the single presentation conditions to align to direction
        # that is closest to unit's preferred direction
        prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
        nullDirIndex = np.where(dirArray == nullDir)[0][0]
        reIndex = (np.array([0, 1, 2, 3, 4, 5])+nullDirIndex) % 6
        loc0ReIndex = b[6, :6][reIndex]
        loc1ReIndex = b[:6, 6][reIndex]
        loc0PrefNormalized.append(loc0ReIndex[3] / np.max(b))
        loc1PrefNormalized.append(loc1ReIndex[3] / np.max(b))
        respAvg = ((loc1ReIndex + loc0ReIndex) / 2)
        baselineNormalized = b[6, 6] / np.max(respAvg)
        respNormalized = respAvg / np.max(respAvg)
        respNormalized = np.concatenate((respNormalized,
                                         respNormalized[:1],
                                         [np.array(unitGaussMean[unitCount])],
                                         [np.array(baselineNormalized)]),
                                        axis=0).tolist()
        popTunCurve.append(respNormalized)
        loc0to1RespRatio.append(loc0ReIndex[3]/loc1ReIndex[3])

        # raw data reindex to have null in the top left corner
        bReIndex = np.zeros((7, 7))
        tempMain = b[:6, :6][:, reIndex]
        tempMain = tempMain[:6, :6][reIndex, :]
        temp0Blank = b[:6, 6][reIndex]
        temp1Blank = b[6, :6][reIndex]
        bReIndex[:6, :6] = tempMain
        bReIndex[:6, 6] = temp0Blank
        bReIndex[6, :6] = temp1Blank
        bReIndex[6, 6] = b[6, 6]
        loc0NMI = (bReIndex[3, 3] + bReIndex[0, 0]) / bReIndex[3, 0]
        loc1NMI = (bReIndex[3, 3] + bReIndex[0, 0]) / bReIndex[0, 3]
        loc0NullNMI.append(loc0NMI)
        loc1NullNMI.append(loc1NMI)

    # the different pairs of neurons from total units
    combs = [i for i in combinations(units, 2)]
    totCombs.append(len(combs))

    # filtered pairs for pairs of units having R2 > 0.75
    filterCombs = [i for i in combinations(filterUnits, 2)]
    totFilterCombs.append(len(filterCombs))

    # combinations of units separated by more than 250ums
    distCombs = []
    for i in combs:
        n1 = unitsChannel[np.where(units == i[0])[0][0]]
        n2 = unitsChannel[np.where(units == i[1])[0][0]]
        n1ElectrodePos = np.where(electrodeArr == n1)[0][0]
        n2ElectrodePos = np.where(electrodeArr == n2)[0][0]
        if abs(n1ElectrodePos - n2ElectrodePos) >= 3:
            distCombs.append(i)

    # generate paired stimulus correlation, selectivity index,
    # suppression, and NMI index for each gabor pair for each pair
    # of neurons (similar to Bram fig 2c)
    # initialize lists for paired Gabors
    pairPairedCorr = []
    pairSelectivityIndex = []
    pairNonPrefSuppIndex = []
    pairNMIIndex = []
    pairSimpleNMIIndex = []

    # initialize lists for single Gabor presentations
    pairSingleCorr = []
    pairSingleSelectivityIndex = []
    pairSingleNonPrefSuppIndex = []

    # arrays for generating 12 orthogonal stim pairings from full matrix
    tempSel12FromFull1 = []
    tempSel12FromFull2 = []
    tempSel12FromFull3 = []
    tempSupp12FromFull1 = []
    tempSupp12FromFull2 = []
    tempSupp12FromFull3 = []

    tempMat = stats.zscore(spikeCountMat[:, :blocksDone, :], axis=1, nan_policy='omit')
    stimIndexMat = np.arange(36).reshape(6, 6)
    upperTriangle = upperTriMasking(stimIndexMat)

    for pairCount, pair in enumerate(distCombs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]

        # indices and correlations for paired Gabor stimuli
        for i in range(36):
            # extract directions of the gabor pair
            loc0Dir = stimIndexDict[i][0]['direction']
            loc1Dir = stimIndexDict[i][1]['direction']

            # correlation for that Gabor pair b/w 2 units excluding trials where
            # spike counts exceeded 3 SD from mean
            n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
            n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
            pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

            # n1 selectivity and suppression index
            loc0Resp = unitNormFitEstimate[n1][np.where(dirArray == loc0Dir)[0][0]]
            loc1Resp = unitNormFitEstimate[n1][np.where(dirArray == loc1Dir)[0][0] + 6]
            al0 = 1
            al1 = unitNormFitEstimate[n1][12]
            n1Selectivity, n1NonPrefSupp = getSelAndSuppIndx(loc0Resp, loc1Resp, al0, al1)
            # n1Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
            # if n1Selectivity >= 0:
            #     n1NonPrefSupp = (unitNormFitEstimate[n1][12]) / (
            #             1 + unitNormFitEstimate[n1][12])
            #     # n1NonPrefSupp = (unitNormFitEstimate[n1][13]) / (
            #     #         unitNormFitEstimate[n1][13] + unitNormFitEstimate[n1][12])
            # else:
            #     n1NonPrefSupp = 1 / (
            #             1 + unitNormFitEstimate[n1][12])
            #     # n1NonPrefSupp = (unitNormFitEstimate[n1][12]) / (
            #     #         unitNormFitEstimate[n1][13] + unitNormFitEstimate[n1][12])

            # n2 selectivity and suppression index
            loc0Resp = unitNormFitEstimate[n2][np.where(dirArray == loc0Dir)[0][0]]
            loc1Resp = unitNormFitEstimate[n2][np.where(dirArray == loc1Dir)[0][0] + 6]
            al0 = 1
            al1 = unitNormFitEstimate[n2][12]
            n2Selectivity, n2NonPrefSupp = getSelAndSuppIndx(loc0Resp, loc1Resp, al0, al1)

            # pair selectivity, suppression, and NMI index
            pairSelectivity = (np.sign(n1Selectivity) * np.sign(n2Selectivity) *
                               np.sqrt(abs(n1Selectivity) * abs(n2Selectivity)))
            pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)
            # pairSuppression = abs(n1NonPrefSupp - n2NonPrefSupp) / (
            #                      (n1NonPrefSupp + n2NonPrefSupp))
            pairNMI = 1
            pairSimpleNMI = 1

            pairPairedCorr.append(pairStimCorr)
            pairSelectivityIndex.append(pairSelectivity)
            pairNonPrefSuppIndex.append(pairSuppression)
            pairNMIIndex.append(pairNMI)
            pairSimpleNMIIndex.append(pairSimpleNMI)

            if i in stimIndex12Cond1:
                tempSel12FromFull1.append(pairSelectivity)
                tempSupp12FromFull1.append(pairSuppression)

            if i in stimIndex12Cond2:
                tempSel12FromFull2.append(pairSelectivity)
                tempSupp12FromFull2.append(pairSuppression)

            if i in stimIndex12Cond3:
                tempSel12FromFull3.append(pairSelectivity)
                tempSupp12FromFull3.append(pairSuppression)

            # loc 0 single Gabor corr
            stimIndex = 36 + np.where(dirArray == loc0Dir)[0][0]
            # correlation for that Gabor pair b/w 2 units excluding trials where
            # spike counts exceeded 3 SD from mean
            n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
            n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
            singleStimLoc0Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

            pairSingleCorr.append(singleStimLoc0Corr)
            pairSingleSelectivityIndex.append(pairSelectivity)
            pairSingleNonPrefSuppIndex.append(pairSuppression)

            # loc 1 single Gabor corr
            stimIndex = 42 + np.where(dirArray == loc1Dir)[0][0]
            # correlation for that Gabor pair b/w 2 units excluding trials where
            # spike counts exceeded 3 SD from mean
            n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
            n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
            singleStimLoc1Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

            pairSingleCorr.append(singleStimLoc1Corr)
            pairSingleSelectivityIndex.append(pairSelectivity)
            pairSingleNonPrefSuppIndex.append(pairSuppression)

    pairSel12FromFull.append(tempSel12FromFull1)
    pairSel12FromFull.append(tempSel12FromFull2)
    pairSel12FromFull.append(tempSel12FromFull3)
    pairNPSupp12FromFull.append(tempSupp12FromFull1)
    pairNPSupp12FromFull.append(tempSupp12FromFull2)
    pairNPSupp12FromFull.append(tempSupp12FromFull3)

        # # indices and correlations for single Gabor b/w pair of neuron
        # for j in range(36, 48):
        #
        #     # correlation for that Gabor pair b/w 2 units excluding trials where
        #     # spike counts exceeded 3 SD from mean
        #     skipTrials = []
        #     n1SpikeMat = spikeCountMat[n1, :blocksDone, j]
        #     n2SpikeMat = spikeCountMat[n2, :blocksDone, j]
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
        #     singleStimCorr = stats.pearsonr(spikeCountMat[n1, goodTrials, j],
        #                                   spikeCountMat[n2, goodTrials, j])
        #
        #     # extract contrast of the gabor at each site
        #     loc0Con = stimIndexDict[j][0]['contrast']
        #     loc1Con = stimIndexDict[j][1]['contrast']
        #
        #     # n1 selectivity and suppression index based of Norm params
        #     loc0MaxResp = max(unitNormFitEstimate[n1][:6])
        #     loc0MinResp = min(unitNormFitEstimate[n1][:6])
        #     loc1MaxResp = max((unitNormFitEstimate[n1][6:12]))
        #     loc1MinResp = min(unitNormFitEstimate[n1][6:12])
        #     minRespN1 = np.where(unitNormFitEstimate[n1][:12] ==
        #                          min(unitNormFitEstimate[n1][:12]))[0]
        #     maxRespN1 = np.where(unitNormFitEstimate[n1][:12] ==
        #                          max(unitNormFitEstimate[n1][:12]))[0]
        #     if maxRespN1 < 6:
        #         n1PrefResp = 0
        #     else:
        #         n1PrefResp = 1
        #     if minRespN1 < 6:
        #         n1NonPrefSupp = unitNormFitEstimate[n1][12] / (
        #                 unitNormFitEstimate[n1][12] + unitNormFitEstimate[n1][13])
        #     else:
        #         n1NonPrefSupp = unitNormFitEstimate[n1][13] / (
        #                 unitNormFitEstimate[n1][12] + unitNormFitEstimate[n1][13])
        #     if loc0Con == 1:
        #         n1Selectivity = (loc0MaxResp - loc0MinResp) / (
        #                          loc0MaxResp + loc0MinResp)
        #     else:
        #         n1Selectivity = (loc1MaxResp - loc1MinResp) / (
        #                 loc1MaxResp + loc1MinResp)
        #
        #     # n2 selectivity and suppression index based of Norm params
        #     loc0MaxResp = max(unitNormFitEstimate[n2][:6])
        #     loc0MinResp = min(unitNormFitEstimate[n2][:6])
        #     loc1MaxResp = max((unitNormFitEstimate[n2][6:12]))
        #     loc1MinResp = min(unitNormFitEstimate[n2][6:12])
        #     minRespN2 = np.where(unitNormFitEstimate[n2][:12] ==
        #                          min(unitNormFitEstimate[n2][:12]))[0]
        #     maxRespN2 = np.where(unitNormFitEstimate[n2][:12] ==
        #                          max(unitNormFitEstimate[n2][:12]))[0]
        #     if maxRespN2 < 6:
        #         n2PrefResp = 0
        #     else:
        #         n2PrefResp = 1
        #     if minRespN2 < 6:
        #         n2NonPrefSupp = unitNormFitEstimate[n2][12] / (
        #                 unitNormFitEstimate[n2][12] + unitNormFitEstimate[n2][13])
        #     else:
        #         n2NonPrefSupp = unitNormFitEstimate[n2][13] / (
        #                 unitNormFitEstimate[n2][12] + unitNormFitEstimate[n2][13])
        #     if loc0Con == 1:
        #         n2Selectivity = (loc0MaxResp - loc0MinResp) / (
        #                          loc0MaxResp + loc0MinResp)
        #     else:
        #         n2Selectivity = (loc1MaxResp - loc1MinResp) / (
        #                 loc1MaxResp + loc1MinResp)
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

    pairPairedCorr = np.array(pairPairedCorr)
    pairSelectivityIndex = np.array(pairSelectivityIndex)
    pairNonPrefSuppIndex = np.array(pairNonPrefSuppIndex)
    pairNMIIndex = np.array(pairNMIIndex)
    pairSimpleNMIIndex = np.array(pairSimpleNMIIndex)
    pairSingleCorr = np.array(pairSingleCorr)
    pairSingleSelectivityIndex = np.array(pairSingleSelectivityIndex)
    pairSingleNonPrefSuppIndex = np.array(pairSingleNonPrefSuppIndex)

    # compile similarity scores and correlation matrix into one matrix
    gaborPairCompiledArr = np.array([pairPairedCorr,
                                     pairSelectivityIndex,
                                     pairNonPrefSuppIndex,
                                     pairNMIIndex,
                                     pairSimpleNMIIndex])
    gaborSingleCompiledArr = np.array([pairSingleCorr,
                                       pairSingleSelectivityIndex,
                                       pairSingleNonPrefSuppIndex])
    np.save(f'../../gaborSingleCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborSingleCompiledArr)
    np.save(f'../../gaborPairCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborPairCompiledArr)
    os.chdir('../../../Python Code')

    totUnits.append(len(units))
    totFilterUnits.append(len(filterUnits))

print(time.time()-t0)


"""
STOP HERE
"""

pairSel12FromFull = np.array([item for sublist in pairSel12FromFull for item in sublist])
pairNPSupp12FromFull = np.array([item for sublist in pairNPSupp12FromFull for item in sublist])

# plotting sel/sel and supp/supp from both matrices
from scipy.stats import gaussian_kde
plt.scatter(pairSelFromCond, pairSel12FromFull)
plt.xlabel('selectivity for each pair, condensed matrix')
plt.ylabel('selectivity for each pair, from full matrix')
plt.title(f'{stats.pearsonr(pairSelFromCond, pairSel12FromFull)[0]}')
plt.show()

plt.scatter(pairNPSuppFromCond, pairNPSupp12FromFull)
plt.xlabel('NP Supp for each pair from condensed matrix')
plt.ylabel('NP Supp for each pair from full matrix')
plt.title(f'{stats.pearsonr(pairNPSuppFromCond, pairNPSupp12FromFull)[0]}')
plt.show()

# density plots
x = np.array(pairSelFromCond)
y = np.array(pairSel12FromFull)
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=100)
plt.show()

x = pairNPSuppFromCond
y = pairNPSupp12FromFull
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=100)
plt.show()

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
#     goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
#     pairStimCorr = stats.pearsonr(spikeCountMat[n1, goodTrials, i],
#                                   spikeCountMat[n2, goodTrials, i])
#
#     # extract directions of the gabor pair
#     loc0Dir = stimIndexDict[i][0]['direction']
#     loc1Dir = stimIndexDict[i][1]['direction']
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
#     if i in stimIndex12Cond1:
#         tempSel12FromFull1.append(pairSelectivity)
#         tempSupp12FromFull1.append(pairSuppression)
#
#     if i in stimIndex12Cond2:
#         tempSel12FromFull2.append(pairSelectivity)
#         tempSupp12FromFull2.append(pairSuppression)
#
#     if i in stimIndex12Cond3:
#         tempSel12FromFull3.append(pairSelectivity)
#         tempSupp12FromFull3.append(pairSuppression)
#
#     # loc 0 single Gabor corr
#     stimIndex = 36 + np.where(dirArray == loc0Dir)[0][0]
#     # correlation for that Gabor pair b/w 2 units excluding trials where
#     # spike counts exceeded 3 SD from mean
#     skipTrials = []
#     n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
#     n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
#     n1Zscore = stats.zscore(n1SpikeMat)
#     n2Zscore = stats.zscore(n2SpikeMat)
#     n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
#     n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
#     for x in n1SkipTrials:
#         skipTrials.append(x)
#     for x in n2SkipTrials:
#         skipTrials.append(x)
#     goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
#     singleStimLoc0Corr = stats.pearsonr(spikeCountMat[n1, goodTrials, stimIndex],
#                                   spikeCountMat[n2, goodTrials, stimIndex])
#     pairSingleCorr.append(singleStimLoc0Corr[0])
#     pairSingleSelectivityIndex.append(pairSelectivity)
#     pairSingleNonPrefSuppIndex.append(pairSuppression)
#
#     # loc 1 single Gabor corr
#     stimIndex = 42 + np.where(dirArray == loc1Dir)[0][0]
#     # correlation for that Gabor pair b/w 2 units excluding trials where
#     # spike counts exceeded 3 SD from mean
#     skipTrials = []
#     n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
#     n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
#     n1Zscore = stats.zscore(n1SpikeMat)
#     n2Zscore = stats.zscore(n2SpikeMat)
#     n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
#     n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
#     for x in n1SkipTrials:
#         skipTrials.append(x)
#     for x in n2SkipTrials:
#         skipTrials.append(x)
#     goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
#     singleStimLoc1Corr = stats.pearsonr(spikeCountMat[n1, goodTrials, stimIndex],
#                                   spikeCountMat[n2, goodTrials, stimIndex])
#     pairSingleCorr.append(singleStimLoc1Corr[0])
#     pairSingleSelectivityIndex.append(pairSelectivity)
#     pairSingleNonPrefSuppIndex.append(pairSuppression)
