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
alphaLoc0 = []
alphaLoc1 = []
popTunCurve = []

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
    onLatency = 25/1000  # time in MS for counting window latency after stim on
    offLatency = 100/1000  # time in MS for counting window latency after stim off
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
        # if r2 > 0.90:
        #     filterUnits.append(unit)
        #     totFilterR2.append(r2)
        unitGaussMean[unitCount] = params[2] % 360
        unitGaussSig[unitCount] = params[3]

    # for loop to run through the 3 orthogonal axes
    # initialize lists for paired Gabors
    pairPairedCorr = []
    pairSelectivityIndex = []
    pairNonPrefSuppIndex = []
    pairNMIIndex = []
    pairSimpleNMIIndex = []
    # initialize lists for single Gabor presentations
    pairSingleCorr = []
    pairSingleGaborNMI = []
    pairSingleGaborSelectivity = []

    for orthoAxis in range(3):

        # SCIPY CURVEFIT/Minimize
        # scipy curvefit Normalization parameters will also generate population
        # tuning average aligned to preferred direction
        unitPairedNormFit = []
        unitPairedNormR2 = np.zeros(len(units))
        unitNormFitEstimate = [[] for i in range(len(units))]
        sli = [orthoAxis, orthoAxis+3, 6]
        dir1 = dirArray[sli[0]]
        dir2 = dir1 + 180
        condArr = np.array([dir1, dir2, dir1, dir2])
        for unitCount, unit in enumerate(units):
            b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS
            bReIndex = np.zeros((7, 7))

            # fixed (independent) variables - matrix of corresponding stim Indexes
            stimMat = np.zeros((7, 7))
            stimMat[:6, :6] = np.arange(36).reshape(6, 6)
            stimMat[6, :6] = np.arange(36, 42)
            stimMat[:, 6] = np.arange(42, 49)

            # condensed matrices for 2 orthogonal directions only
            bCondensed = b[sli, :]
            bCondensed = bCondensed[:, sli]
            stimMatCond = stimMat[sli, :]
            stimMatCond = stimMatCond[:, sli]

            # Generic Normalization (L1+L2)/(al1+al2+sig) w.o scalar
            # fits L0-L6 for both locations separately loc0 and loc1 will
            # have different L0-L6 values
            # resp = b.reshape(49)[:-1]
            # fixedVals = fixedValsForGenericNorm(stimMat.reshape(49)[:-1], stimIndexDict)
            # pOpt, pCov = curve_fit(genericNormNoScalar, fixedVals, resp.squeeze(),
            #                        bounds=((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            #                                (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
            #                                 np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
            #                                 np.inf, np.inf, np.inf)))
            #
            # y_pred = genericNormNoScalar(fixedVals, *pOpt)
            # r2 = r2_score(resp.squeeze(), y_pred)
            # print(unit, r2)

            # fit using scipy.optimize.minimize
            # p0 = np.concatenate((b[6, :6], b[:6, 6], [1, 1], [0.1]), axis=0)
            # bnds = ((.01, None), (.01, None), (.01, None), (.01, None), (.01, None),
            #         (.01, None), (.01, None), (.01, None), (.01, None), (.01, None),
            #         (.01, None), (.01, None), (0, None), (0, None), (0.01, None))
            # resp = b.reshape(49)[:-1]
            # fixedVals = fixedValsForGenericNorm(stimMat.reshape(49)[:-1], stimIndexDict)
            # res = scipy.optimize.minimize(driverFunc, p0, args=(fixedVals, resp),
            #                               method='Nelder-Mead', bounds=bnds,
            #                               options={'maxfev': 100000})
            # y_pred = genericNormNoScalar(fixedVals, *res.x)
            # r2 = r2_score(resp.squeeze(), y_pred)
            # print(unit, r2)

            # # Append fit parameters for full matrix
            # totR2.append(r2)
            # # alphaLoc0.append(pOpt[12])
            # # alphaLoc1.append(pOpt[13])
            # alphaLoc0.append(res.x[12])
            # alphaLoc1.append(res.x[13])
            # unitPairedNormR2[unitCount] = r2
            # # unitPairedNormFit.append(pOpt)
            # # unitNormFitEstimate[unitCount] = pOpt
            # unitPairedNormFit.append(res.x)
            # unitNormFitEstimate[unitCount] = res.x

            # fit using scipy.optimize.minimize for condensed mat
            # (2 orthogonal dirs only)
            p0 = np.concatenate((bCondensed[2, :-1], bCondensed[:-1, 2],
                                 [1, 1], [0.1]), axis=0)
            bnds = ((0, None), (0, None), (0, None), (0, None),
                    (0, None), (0, None), (0, .15))
            resp = bCondensed.reshape(9)
            fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9),
                                                    stimIndexDict, dir1, dir2)
            res = scipy.optimize.minimize(driverFunc, p0, args=(fixedVals, resp),
                                          method='Nelder-Mead', bounds=bnds,
                                          options={'maxfev': 100000})

            y_pred = genNormCondensed(fixedVals, *res.x)
            r2 = r2_score(resp.squeeze(), y_pred)
            print(unit, r2)

            # Append fit parameters for condensed matrix
            totR2.append(r2)
            alphaLoc0.append(res.x[4])
            alphaLoc1.append(res.x[4])
            unitPairedNormR2[unitCount] = r2
            unitPairedNormFit.append(res.x)
            unitNormFitEstimate[unitCount] = res.x

            if r2 > 0.75:
                filterUnits.append(unit)
                totFilterR2.append(r2)


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

        # the different pairs of neurons from total units
        combs = [i for i in combinations(units, 2)]
        totCombs.append(len(combs))

        # filtered pairs for pairs of units having R2 > 0.75
        filterCombs = [i for i in combinations(filterUnits, 2)]
        totFilterCombs.append(len(filterCombs))

        # generate paired stimulus correlation, selectivity index,
        # suppression, and NMI index for each gabor pair for each pair
        # of neurons (similar to Bram fig 2c)
        tempMat = stats.zscore(spikeCountMat[:, :blocksDone, :], axis=1, nan_policy='omit')
        stimIndexMat = np.arange(36).reshape(6, 6)
        upperTriangle = upperTriMasking(stimIndexMat)

        for pairCount, pair in enumerate(combs):
            n1 = np.where(units == pair[0])[0][0]
            n2 = np.where(units == pair[1])[0][0]

            # unit's preferred directions
            n1PrefDir = unitGaussMean[n1]
            n2PrefDir = unitGaussMean[n2]

            # indices and correlations for paired Gabor stimuli
            for count, i in enumerate(np.int_(stimMatCond[:2, :2].reshape(4))):

                # correlation for that Gabor pair b/w 2 units excluding trials where
                # spike counts exceeded 3 SD from mean
                skipTrials = []
                n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
                n1Zscore = stats.zscore(n1SpikeMat)
                n2Zscore = stats.zscore(n2SpikeMat)
                n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
                n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
                for x in n1SkipTrials:
                    skipTrials.append(x)
                for x in n2SkipTrials:
                    skipTrials.append(x)
                # print(skipTrials)
                goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
                pairStimCorr = stats.pearsonr(spikeCountMat[n1, goodTrials, i],
                                              spikeCountMat[n2, goodTrials, i])

                # extract directions of the gabor pair
                loc0Dir = stimIndexDict[i][0]['direction']
                loc1Dir = stimIndexDict[i][1]['direction']

                # n1 selectivity and suppression index based of Normalization params
                loc0Resp = unitNormFitEstimate[n1][np.where(condArr == loc0Dir)[0][0]]
                loc1Resp = unitNormFitEstimate[n1][np.where(condArr == loc1Dir)[0][0] + 2]
                if loc0Resp > loc1Resp:
                    n1Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
                    n1NonPrefSupp = unitNormFitEstimate[n1][5] / (
                                    unitNormFitEstimate[n1][5] + unitNormFitEstimate[n1][4])
                    n1PrefResp = 0
                else:
                    n1Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
                    n1NonPrefSupp = unitNormFitEstimate[n1][4] / (
                                    unitNormFitEstimate[n1][4] + unitNormFitEstimate[n1][5])
                    n1PrefResp = 1

                # n2 selectivity and suppression index
                loc0Resp = unitNormFitEstimate[n2][np.where(condArr == loc0Dir)[0][0]]
                loc1Resp = unitNormFitEstimate[n2][np.where(condArr == loc1Dir)[0][0] + 2]
                if loc0Resp > loc1Resp:
                    n2Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
                    n2NonPrefSupp = unitNormFitEstimate[n2][5] / (
                                    unitNormFitEstimate[n2][5] + unitNormFitEstimate[n2][4])
                    n2PrefResp = 0
                else:
                    n2Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
                    n2NonPrefSupp = unitNormFitEstimate[n2][4] / (
                                    unitNormFitEstimate[n2][4] + unitNormFitEstimate[n2][5])
                    n2PrefResp = 1

                # n1 NMI
                b = meanSpikeReshaped[n1].reshape(7, 7) * 1000/trueStimDurMS
                bCondensed = b[sli, :]
                bCondensed = bCondensed[:, sli]
                pairedMat = bCondensed[:2, :2].reshape(4)
                loc0Single = bCondensed[2, :2]
                loc1Single = bCondensed[:2, 2]
                loc0SingleResp = loc0Single[np.where(condArr == loc0Dir)[0][0]]
                loc1SingleResp = loc1Single[np.where(condArr == loc1Dir)[0][0]]
                n1NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[count]) / (
                      (loc0SingleResp + loc1SingleResp) + pairedMat[count])
                n1SimpleNMI = pairedMat[count] / (loc0SingleResp + loc1SingleResp)

                # n2 NMI
                b = meanSpikeReshaped[n2].reshape(7, 7) * 1000/trueStimDurMS
                bCondensed = b[sli, :]
                bCondensed = bCondensed[:, sli]
                pairedMat = bCondensed[:2, :2].reshape(4)
                loc0Single = bCondensed[2, :2]
                loc1Single = bCondensed[:2, 2]
                loc0SingleResp = loc0Single[np.where(condArr == loc0Dir)[0][0]]
                loc1SingleResp = loc1Single[np.where(condArr == loc1Dir)[0][0]]
                n2NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[count]) / (
                      (loc0SingleResp + loc1SingleResp) + pairedMat[count])
                n2simpleNMI = pairedMat[count] / (loc0SingleResp + loc1SingleResp)

                # pair selectivity, suppression, and NMI index
                pairSelectivity = np.sqrt(n1Selectivity * n2Selectivity)
                pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)
                if n1PrefResp != n2PrefResp:
                    pairSelectivity = -pairSelectivity
                pairNMI = (n1NMI + n2NMI) / 2
                pairSimpleNMI = np.sqrt(n1SimpleNMI * n2simpleNMI)

                pairPairedCorr.append(pairStimCorr[0])
                pairSelectivityIndex.append(pairSelectivity)
                pairNonPrefSuppIndex.append(pairSuppression)
                pairNMIIndex.append(pairNMI)
                pairSimpleNMIIndex.append(pairSimpleNMI)

                # # appending pair Selectivity and pair NMI to get average
                # # for the single Gabor figure
                # neuronPairAvgSelectivity.append(pairSelectivity)
                # if i in upperTriangle:
                #     neuronPairAvgNMIloc0.append(pairNMI)
                # else:
                #     neuronPairAvgNMIloc1.append(pairNMI)

            # indices and correlations for single Gabor stimuli
            singleStimIndex = np.int_(np.concatenate((stimMatCond[2, :2],
                                                      stimMatCond[:2, 2]),
                                                     axis=0))
            for count, i in enumerate(singleStimIndex):
                # correlation for the single Gabor b/w 2 units excluding trials where
                # spike counts exceeded 3 SD from mean
                skipTrials = []
                n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
                n1Zscore = stats.zscore(n1SpikeMat)
                n2Zscore = stats.zscore(n2SpikeMat)
                n1SkipTrials = np.where(abs(n1Zscore) > 3)[0].tolist()
                n2SkipTrials = np.where(abs(n2Zscore) > 3)[0].tolist()
                for x in n1SkipTrials:
                    skipTrials.append(x)
                for x in n2SkipTrials:
                    skipTrials.append(x)
                # print(skipTrials)
                goodTrials = [x for x in range(blocksDone) if x not in skipTrials]
                pairStimCorr = stats.pearsonr(spikeCountMat[n1, goodTrials, i],
                                              spikeCountMat[n2, goodTrials, i])

                # extract directions of the gabor pair
                loc0Dir = stimIndexDict[i][0]['direction']
                loc1Dir = stimIndexDict[i][1]['direction']



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

    pairPairedCorr = np.array(pairPairedCorr)
    pairSelectivityIndex = np.array(pairSelectivityIndex)
    pairNonPrefSuppIndex = np.array(pairNonPrefSuppIndex)
    pairNMIIndex = np.array(pairNMIIndex)
    pairSimpleNMIIndex = np.array(pairSimpleNMIIndex)
    pairSingleCorr = np.array(pairSingleCorr)
    pairSingleGaborNMI = np.array(pairSingleGaborNMI)
    pairSingleGaborSelectivity = np.array(pairSingleGaborSelectivity)

    # compile similarity scores and correlation matrix into one matrix
    gaborPairCompiledArr = np.array([pairPairedCorr,
                                     pairSelectivityIndex,
                                     pairNonPrefSuppIndex,
                                     pairNMIIndex,
                                     pairSimpleNMIIndex])
    gaborSingleCompiledArr = np.array([pairSingleCorr,
                                       pairSingleGaborNMI,
                                       pairSingleGaborSelectivity])
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

# Plotting

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
    ax5 = fig.add_subplot(gs03[0,0])
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
