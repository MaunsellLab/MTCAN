# imports
import matplotlib.pyplot as plt
from usefulFns import *
from normalizationFunctions import *
from itertools import combinations
from scipy import stats
import time
import numpy as np
import pandas as pd
from scipy.stats import f

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
totR2EMS = []
totR2WeightAlpha = []
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
# population PSTH
popPSTH = []
popPSTHBase = []
PSTH7by7arrays = [[] for i in range(49)]
# population SpikeCounts
popSpikeCountMat7by7 = [[] for i in range(49)]

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
    fullDirArray = np.array([0, 60, 120, 180, 240, 300, 0, 60, 120, 180, 240, 300])
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

    # the different combinations of neuron pairs from total units
    combs = [i for i in combinations(units, 2)]
    # filtered pairs for pairs of units having R2 > 0.75 for gaussian fit to direction tuning
    filterCombs = [i for i in combinations(filterUnits, 2)]
    # combinations of units separated by more than 250ums
    distCombs = []
    for i in combs:
        n1 = unitsChannel[np.where(units == i[0])[0][0]]
        n2 = unitsChannel[np.where(units == i[1])[0][0]]
        n1ElectrodePos = np.where(electrodeArr == n1)[0][0]
        n2ElectrodePos = np.where(electrodeArr == n2)[0][0]
        pairDistOnElectrode = abs(n1ElectrodePos - n2ElectrodePos)
        pairDistance.append(pairDistOnElectrode * 50)
        if pairDistOnElectrode >= 5:
            distCombs.append(i)

    # SCIPY CURVEFIT/Minimize
    # scipy curvefit Normalization parameters will also generate population
    # tuning average aligned to preferred direction
    unitPairedNormFit = []
    unitPairedNormR2 = np.zeros(len(units))
    unitNormFitEstimate = [[] for i in range(len(units))]
    unitWeights = [[] for i in range(len(units))]
    # for unit in filterUnits:
    for unit in units:

        unitCount = np.where(units == unit)[0][0]
        b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS

        # fixed (independent) variables - matrix of corresponding stim Indexes
        stimMat = np.zeros((7, 7))
        stimMat[:6, :6] = np.arange(36).reshape(6, 6)
        stimMat[6, :6] = np.arange(36, 42)
        stimMat[:, 6] = np.arange(42, 49)

        # What I used for Meetz's data
        # Generic Normalization (L1+L2)/(al1+al2+sig) w.o scalar
        # fits L0-L6 for both locations separately loc0 and loc1 will
        # have different L0-L6 values
        # guess0 = np.concatenate((b[6, :6], b[:6, 6], [0.4]), axis=0)
        # resp = b.reshape(49)[:-1]
        # fixedVals = fixedValsForGenericNorm(stimMat.reshape(49)[:-1], stimIndexDict)
        # pOpt, pCov = curve_fit(genericNormNoScalar, fixedVals, resp.squeeze(), p0=guess0,
        #                        bounds=((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        #                                (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 6)))
        # # pOpt, pCov = curve_fit(genericNormNoScalar, fixedVals, resp.squeeze(), p0=guess0,
        # #                        maxfev=10000000)
        # y_pred = genericNormNoScalar(fixedVals, *pOpt)
        # r2 = r2_score(resp.squeeze(), y_pred)

        # bram trick to keep things non negative
        # pOpt = pOpt ** 2

        # Bram Normalization (L1+L2)/(al1+al2+sig)
        # fits L0-L6 for both locations separately loc0 and loc1 will
        # have different L0-L6 values, al1 = 1 in this case
        guess0 = np.concatenate((b[6, :6], b[:6, 6], [0.8], [0.8], [0.1], [b[6, 6]]), axis=0)
        resp = b.reshape(49)
        fixedVals = fixedValsForGenericNorm(stimMat.reshape(49), stimIndexDict)
        # pOpt, pCov = curve_fit(bramNormMTNC, fixedVals, resp.squeeze(), p0=guess0,
        #                        bounds=((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        #                                (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 np.inf, 0.1, np.inf)))
        pOpt, pCov = curve_fit(bramNormMTNC, fixedVals, resp.squeeze(), p0=guess0,
                               maxfev=1000000)
        y_pred = bramNormMTNC(fixedVals, *pOpt)
        r2 = r2_score(resp.squeeze(), y_pred)

        # EMS Normalization (L1/c1+c2al2+sig) + (L2/c1al1+c2+sig)
        # fits L0-L6 for both locations with the same value,
        # one L0-L6 value for both locations
        guess0 = np.concatenate((np.mean((b[6, :6], b[:6, 6]), axis=0),
                                 [1], [1], [0.1], [b[6, 6]]), axis=0)
        resp = b.reshape(49)
        fixedVals = fixedValsForGenericNorm(stimMat.reshape(49), stimIndexDict)
        # pOpt, pCov = curve_fit(EMSNormMTNC, fixedVals, resp.squeeze(), p0=guess0,
        #                        bounds=((-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,
        #                                 -10, -10, -np.inf, -np.inf),
        #                                (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 10, 10, np.inf, np.inf)))
        pOpt, pCov = curve_fit(EMSNormMTNC, fixedVals, resp.squeeze(), p0=guess0,
                               maxfev=1000000)
        y_pred = EMSNormMTNC(fixedVals, *pOpt)
        r2EMS = r2_score(resp.squeeze(), y_pred)

        # # RF Weighted Norm
        prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
        prefIndex = np.where(dirArray == prefDir)[0][0]
        rfWeightLoc1 = b[prefIndex, 6] / b[6, prefIndex]
        resp = b.reshape(49)
        fixedVals = fixedValsForRFWeightMTNC(stimMat.reshape(49), stimIndexDict, rfWeightLoc1)
        # fixedVals = fixedValsForGenericNorm(stimMat.reshape(49), stimIndexDict)
        guess0 = np.concatenate((b[6, :6], b[:6, 6], [0.1], [b[6, 6]]), axis=0)
        # pOptRFWeight, pCov = curve_fit(rfWeightMTNC, fixedVals, resp.squeeze(), p0=guess0,
        #                        bounds=((0, 0, 0, 0, 0, 0,
        #                                 0, 0, 0, 0, 0, 0,
        #                                 0, 0),
        #                                (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 0.1, np.inf)))
        pOptRFWeight, pCov = curve_fit(rfWeightMTNC, fixedVals, resp.squeeze(), p0=guess0,
                               maxfev=10000000)
        y_pred = rfWeightMTNC(fixedVals, *pOptRFWeight)
        r2RFWeight = r2_score(resp.squeeze(), y_pred)

        # # RF Weighted Norm with Alpha (extra parameter)
        prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
        prefIndex = np.where(dirArray == prefDir)[0][0]
        rfWeightLoc1 = b[prefIndex, 6] / b[6, prefIndex]
        resp = b.reshape(49)
        fixedVals = fixedValsForRFWeightMTNC(stimMat.reshape(49), stimIndexDict, rfWeightLoc1)
        guess0 = np.concatenate((b[6, :6], b[:6, 6], [0.1], [1], [b[6, 6]]), axis=0)
        # pOpt, pCov = curve_fit(rfWeightMTNC, fixedVals, resp.squeeze(), p0=guess0,
        #                        bounds=((-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,
        #                                 -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,
        #                                 0, -np.inf),
        #                                (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
        #                                 0.1, np.inf)))
        pOpt, pCov = curve_fit(rfWeightWithAlphaMTNC, fixedVals, resp.squeeze(), p0=guess0,
                               maxfev=10000000)
        y_pred = rfWeightWithAlphaMTNC(fixedVals, *pOpt)
        r2RFWeightAlpha = r2_score(resp.squeeze(), y_pred)

        totR2Bram.append(r2)
        totR2RFWeight.append(r2RFWeight)
        totR2WeightAlpha.append(r2RFWeightAlpha)
        totR2EMS.append(r2EMS)

        # bram trick to keep values > 0
        pOptRFWeight = pOptRFWeight ** 2

        # Append fit parameters for full matrix
        totR2.append(r2)
        unitPairedNormR2[unitCount] = r2
        unitPairedNormFit.append(pOpt)
        unitNormFitEstimate[unitCount] = pOptRFWeight
        unitWeights[unitCount] = [1, rfWeightLoc1]

    # generate paired stimulus correlation, selectivity index,
    # and suppression index for each gabor pair for each pair
    # of neurons (similar to Bram fig 2c)

    # initialize lists for paired Gabors
    pairPairedCorr = []
    pairDoubleCov = []
    pairDoubleSD = []
    pairSelectivityIndex = []
    pairNonPrefSuppIndex = []
    # initialize lists for single Gabor presentations
    pairSingleCorr = []
    pairSingleCov = []
    pairSingleSD = []
    pairSingleSelectivityIndex = []
    pairSingleNonPrefSuppIndex = []
    # for pairCount, pair in enumerate(filterCombs):
    for pairCount, pair in enumerate(distCombs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]

        # unit's preferred directions
        n1PrefDir = unitGaussMean[n1]
        n2PrefDir = unitGaussMean[n2]

        # indices and correlations for paired Gabor stimuli
        for count, i in enumerate(np.int_(stimMat[:6, :6].reshape(36))):

            # correlation for that Gabor pair b/w 2 units excluding trials where
            # spike counts exceeded 3 SD from mean
            n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
            n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
            pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

            # extract directions of the gabor pair
            loc0Dir = stimIndexDict[i][0]['direction']
            loc1Dir = stimIndexDict[i][1]['direction']

            # bram selectivity recreation
            # n1 selectivity and suppression index
            loc0Resp = unitNormFitEstimate[n1][np.where(fullDirArray == loc0Dir)[0][0]]
            loc1Resp = unitNormFitEstimate[n1][np.where(fullDirArray == loc1Dir)[0][1]]
            weightLoc0 = 1
            weightLoc1 = unitWeights[n1][1]
            n1Selectivity, n1NonPrefSupp = getSelAndSuppIndx(loc0Resp, loc1Resp, weightLoc0, weightLoc1)

            # n2 selectivity and suppression index
            loc0Resp = unitNormFitEstimate[n2][np.where(fullDirArray == loc0Dir)[0][0]]
            loc1Resp = unitNormFitEstimate[n2][np.where(fullDirArray == loc1Dir)[0][1]]
            weightLoc0 = 1
            weightLoc1 = unitWeights[n2][1]
            n2Selectivity, n2NonPrefSupp = getSelAndSuppIndx(loc0Resp, loc1Resp, weightLoc0, weightLoc1)

            # # pair selectivity, suppression, and NMI index
            pairSelectivity = (np.sign(n1Selectivity) * np.sign(n2Selectivity) *
                               np.sqrt(abs(n1Selectivity) * abs(n2Selectivity)))
            # pairSuppression = (abs(n1NonPrefSupp - n2NonPrefSupp)) / (
            #                       n1NonPrefSupp + n2NonPrefSupp)
            pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)

            pairPairedCorr.append(pairStimCorr)
            pairDoubleCov.append(pairDCov[0][1])
            pairDoubleSD.append(pairDSD)
            pairSelectivityIndex.append(pairSelectivity)
            pairNonPrefSuppIndex.append(pairSuppression)

            # loc 0 single Gabor corr, excluding trials where spike counts > 3 SD
            stimIndex = np.int_(stimMat[6, :][np.where(fullDirArray == loc0Dir)[0][0]])

            n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
            n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
            singleStimLoc0Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

            pairSingleCov.append(pairSCov[0][1])
            pairSingleSD.append(pairSSD)
            pairSingleCorr.append(singleStimLoc0Corr)
            pairSingleSelectivityIndex.append(pairSelectivity)
            pairSingleNonPrefSuppIndex.append(pairSuppression)

            # loc 1 single Gabor corr, excluding trials where spike counts > 3 SD
            stimIndex = np.int_(stimMat[:, 6][np.where(fullDirArray == loc1Dir)[0][0]])

            n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
            n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
            singleStimLoc1Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

            pairSingleCov.append(pairSCov[0][1])
            pairSingleSD.append(pairSSD)
            pairSingleCorr.append(singleStimLoc1Corr)
            pairSingleSelectivityIndex.append(pairSelectivity)
            pairSingleNonPrefSuppIndex.append(pairSuppression)

    pairDoubleCov = np.array(pairDoubleCov)
    pairDoubleSd = np.array(pairDoubleSD)
    pairPairedCorr = np.array(pairPairedCorr)
    pairSelectivityIndex = np.array(pairSelectivityIndex)
    pairNonPrefSuppIndex = np.array(pairNonPrefSuppIndex)

    pairSingleCov = np.array(pairSingleCov)
    pairSingleSD = np.array(pairSingleSD)
    pairSingleCorr = np.array(pairSingleCorr)
    pairSingleSelectivityIndex = np.array(pairSingleSelectivityIndex)
    pairSingleNonPrefSuppIndex = np.array(pairSingleNonPrefSuppIndex)

    # compile similarity scores and correlation matrix into one matrix
    gaborPairCompiledArr = np.array([pairPairedCorr,
                                     pairSelectivityIndex,
                                     pairNonPrefSuppIndex,
                                     pairDoubleCov,
                                     pairDoubleSD])
    gaborSingleCompiledArr = np.array([pairSingleCorr,
                                       pairSingleSelectivityIndex,
                                       pairSingleNonPrefSuppIndex,
                                       pairSingleCov,
                                       pairSingleSD])

    np.save(f'../../gaborSingleCorrMasterRFWeight/pairCorrelationsAndIndices{seshDate}',
            gaborSingleCompiledArr)
    np.save(f'../../gaborPairCorrMasterRFWeight/pairCorrelationsAndIndices{seshDate}',
            gaborPairCompiledArr)

    # add normalized PSTH for each condition from each filtered unit to get population PSTH
    # this will also generate a 7x7 matrix of PSTH's where each axis is the direction at loc 0 and 1.
    # this will enable us to view the PSTH for single vs paired conditions (the 7x7 matrix is aligned
    # to each unit's preferred direction)
    for unit in filterUnits:
        unitCount = np.where(units == unit)[0][0]
        maxResp = 0
        for i in range(49):
            dirPlot = spikeHists[unitCount, i, :] * 1000 / stimIndexCount[i]
            if max(dirPlot) > maxResp:
                maxResp = max(dirPlot)
        for i in range(49):
            # normalized by max response condition for that neuron
            dirPlot = (spikeHists[unitCount, i, :] * 1000 / stimIndexCount[i]) / maxResp
            popPSTH.append(dirPlot)
            if i == 48:
                # smoothPlot = gaussian_filter1d(dirPlot, 5)
                # smoothPlot = smoothPlot / np.max(smoothPlot)
                popPSTHBase.append(dirPlot)

        # reindex the single presentation conditions to align to direction
        # that is closest to unit's preferred direction (PSTH)
        b = spikeHists[unitCount].reshape(7, 7, trueStimDurMS + (2*histPrePostMS+1))
        c = stimIndexCount.reshape(7, 7)
        c = c[:, :, np.newaxis]
        b = b * 1000 / c

        prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
        nullDirIndex = np.where(dirArray == nullDir)[0][0]
        reIndex = (np.array([0, 1, 2, 3, 4, 5])+nullDirIndex) % 6

        # Create a new 7x7x648 array filled with zeros
        bReIndex = np.zeros_like(b)
        # Step 1: Reindex the main 6x6x648 block
        tempMain = b[:6, :6, :][:, reIndex, :]
        tempMain = tempMain[reIndex, :, :]
        # Step 2: Reindex the rightmost column (6th column of the first 6 rows) across all depths
        temp0Blank = b[:6, 6, :][reIndex, :]
        # Step 3: Reindex the bottom row (6th row of the first 6 columns) across all depths
        temp1Blank = b[6, :6, :][reIndex, :]
        # Step 4: Assign the reindexed blocks back into the new array
        bReIndex[:6, :6, :] = tempMain
        bReIndex[:6, 6, :] = temp0Blank
        bReIndex[6, :6, :] = temp1Blank
        # Assign the bottom-right corner across all depths
        bReIndex[6, 6, :] = b[6, 6, :]
        # normalized bReIndex
        bReIndex = bReIndex / maxResp #np.max(meanSpikeReshaped[unitCount] * 1000/trueStimDurMS)

        # reIndex the spikeCount
        b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000 / trueStimDurMS
        bReIndexSpikeCount = np.zeros_like(b)

        tempMain = b[:6, :6][:, reIndex]
        tempMain = tempMain[:6, :6][reIndex, :]
        temp0Blank = b[:6, 6][reIndex]
        temp1Blank = b[6, :6][reIndex]
        bReIndexSpikeCount[:6, :6] = tempMain
        bReIndexSpikeCount[:6, 6] = temp0Blank
        bReIndexSpikeCount[6, :6] = temp1Blank
        bReIndexSpikeCount[6, 6] = b[6, 6]
        bReIndexSpikeCount = bReIndexSpikeCount / np.max(bReIndexSpikeCount)

        count = 0
        for x in range(7):
            for y in range(7):
                PSTH7by7arrays[count].append(bReIndex[x, y])
                popSpikeCountMat7by7[count].append(bReIndexSpikeCount[x, y])
                count += 1

    os.chdir('../../../Python Code')

print(time.time()-t0)

############################################################

#################### Plotting ##############################

############################################################

# plotting Bram, RF Weight, and EMS R2 values
totR2EMS = np.array(totR2EMS)
totR2Bram = np.array(totR2Bram)
totR2RFWeight = np.array(totR2RFWeight)
totR2WeightAlpha = np.array(totR2WeightAlpha)
for j in range(1):
    plotData = [totR2EMS, totR2Bram, totR2RFWeight, totR2WeightAlpha]
    # plotData = [totR2EMS, [], totR2RFWeight]
    fig, axes = plt.subplots(1, 3, figsize=(12, 6))

    # Plot the box plots for Norm R2 values
    sns.boxplot(data=plotData, ax=axes[0])
    axes[0].set_title('Normalization Model Comparisons for MTNC')
    axes[0].set_xlabel('Model')
    axes[0].set_ylabel('R2 Value Boxplots')
    axes[0].set_xticks([0, 1, 2, 3])
    axes[0].set_xticklabels(['EMS Norm', 'Bram Norm', 'RF Weight Norm', 'RF Weight w. Alpha'], rotation=90)
    axes[0].set_ylim(0, 1.2)

    for i, dataset in enumerate(plotData):
        median = np.median(dataset)
        axes[0].text(i, median, f'{median:.2f}', horizontalalignment='center', color='black', weight='bold')

    # plot scatter of RF weight R2 vs EMS R2
    indx = np.where(totR2RFWeight >= 0)[0]
    a = totR2RFWeight[indx]
    b = totR2EMS[indx]
    axes[1].scatter(a, b)
    axes[1].set_xlabel('RF Weight R2 values')
    axes[1].set_ylabel('EMS R2 values')
    min_val = min(a.min(), b.min())
    max_val = max(a.max(), b.max())
    axes[1].plot([min_val, max_val], [min_val, max_val], 'r--')

    # plot scatter of RF weight R2 vs RF Weight w. Alpha R2
    indx = np.where(totR2RFWeight >= 0)[0]
    a = totR2RFWeight[indx]
    b = totR2WeightAlpha[indx]
    axes[2].scatter(a, b)
    axes[2].set_xlabel('RF Weight R2 values')
    axes[2].set_ylabel('RF Weight w. Alpha R2 values')
    min_val = min(a.min(), b.min())
    max_val = max(a.max(), b.max())
    axes[2].plot([min_val, max_val], [min_val, max_val], 'r--')

    plt.tight_layout()
    plt.show()

    # Perform Wilcoxon rank-sum test (Mann-Whitney U test) between totR2EMS and totR2RFWeight
    stat, p = mannwhitneyu(totR2WeightAlpha, totR2RFWeight)
    # Display test results
    print(f'Mann-Whitney U test between totR2EMS and totR2RFWeight: U={stat}, p-value={p}')

    # perform an f-test to see whether additional parameter in RF weight statistically
    # improves the fit
    numNeurons = len(totR2RFWeight)
    RSS1 = (1 - totR2RFWeight) * numNeurons
    RSS2 = (1 - totR2WeightAlpha) * numNeurons

    # Degrees of freedom for each model
    # Simpler model has p1 parameters, more complex model has p2 parameters
    p1 = 14
    p2 = 15

    numerator = (RSS1 - RSS2) / (p2 - p1)
    denominator = RSS2 / (numNeurons - p2)
    F_statistic = numerator / denominator
    p_value = f.sf(F_statistic, p2 - p1, numNeurons - p2)


# plotting population PSTH
popPSTH = np.array(popPSTH)
popPSTHBase = np.array(popPSTHBase)
avgPSTH = np.nanmean(popPSTH, axis=0)
avgPSTHBase = np.nanmean(popPSTHBase, axis=0)

fig, ax = plt.subplots()
ax.plot(avgPSTH)
ax.plot(avgPSTHBase, linestyle='--', color='grey')
ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
               2 * histPrePostMS + trueStimDurMS])
ax.set_xticklabels([])
ax.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
           color='grey', alpha=0.1)
ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
ax.set_xlabel('Stimulus Duration (ms)', fontsize=7)
ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                    trueStimDurMS + histPrePostMS],
                   fontsize=7)
ax.set_title('population PSTH for all conditions, hist pre/post = 200ms')
plt.show()


# plotting population PSTH of 7x7 matrix
PSTH7by7arrays = np.array(PSTH7by7arrays)
PSTH7by7Avg = np.mean(PSTH7by7arrays, axis=1)
maxYVal = np.max(PSTH7by7Avg)

fig, axes = plt.subplots(7, 7, figsize=(10, 10))
# Iterate over each element in the 7x7 array
count = 0
for i in range(7):
    for j in range(7):
        ax = axes[i, j]
        # smoothPlot = gaussian_filter1d(PSTH7by7Avg[count], 5)
        # ax.plot(smoothPlot)
        ax.plot(PSTH7by7Avg[count])
        ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                       2 * histPrePostMS + trueStimDurMS])
        ax.set_xticklabels([])
        ax.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
        ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                   color='grey', alpha=0.1)
        ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                            trueStimDurMS + histPrePostMS],
                           fontsize=7)
        ax.set_ylim([0, maxYVal])
        if i == 6 and j == 0:
            ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
            ax.set_xlabel('Stimulus Duration (ms)', fontsize=7)
            # ax.set_title('population PSTH for all conditions, hist pre/post = 200ms')
        count += 1

plt.tight_layout()
plt.show()


popSpikeCountMat7by7 = np.array(popSpikeCountMat7by7)
popSpikeCountAvg = np.mean(popSpikeCountMat7by7, axis=1)
popSpikeCountAvg = popSpikeCountAvg.reshape(7, 7)
fig, ax = plt.subplots()
tickLabels = np.array(['0', '60', '120', '180', '240', '300', 'blank'])

ax = sns.heatmap(popSpikeCountAvg, square=True, linewidths=0.2, vmin=0,
                 annot=True, annot_kws={'fontsize': 7}, cbar=False)
ax.set_xticks(np.arange(7) + 0.5)
ax.set_title(f'Raw Data', y=-0.1, fontsize=7)
ax.set_xlabel('Location 0 Stimulus Direction', fontsize=7)
ax.xaxis.set_label_position('top')
ax.set_xticklabels(tickLabels, rotation=45, fontsize=7)
ax.set_ylabel('Location 1 Stimulus Direction', fontsize=7)
ax.xaxis.set_ticks_position("top")
ax.set_yticks(np.arange(7) + 0.5)
ax.set_yticklabels(tickLabels, rotation=0, fontsize=7)

plt.show()



## one tailed t-test and then wilcoxon signed rank test (stats.wilcoxin)

fig, ax = plt.subplots()
ax.plot(PSTH7by7Avg[0], color='grey')
ax.plot(PSTH7by7Avg[3], color='pink')
ax.plot(PSTH7by7Avg[21], color='pink')
ax.plot(PSTH7by7Avg[24], color='black')
plt.show()

