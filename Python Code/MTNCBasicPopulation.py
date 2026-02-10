# imports
import matplotlib.pyplot as plt
from usefulFns import *
from normalizationFunctions import *
from itertools import combinations
from scipy import stats
import time
import numpy as np
import pandas as pd

################################################ MEETZ #################################################################
# good sessions: 221110, 221115, 221117, 221124, 221128, 221208, 221229, 230123, 230126
# okay sessions: 221010, 221013, 221108, 221206
# bad sessions: 230124
#
# # # for loop to run through all good files
# fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
#             'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
#             'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
#             'Meetz_230126']

# # individual file
# fileList = ['Meetz_221010']

################################################ AKSHAN ################################################################

# good sessions: 240927, 240930, 241016, 241017, 241023, 241025
# okay sessions: 240826, 241002 (13 blocks), 241021 (17 blocks)
# bad sessions: 240827, 240828
#
# fileList = ['Akshan_240826', 'Akshan_240927', 'Akshan_240930', 'Akshan_241002',
#             'Akshan_241016', 'Akshan_241017', 'Akshan_241021', 'Akshan_241023',
#             'Akshan_241025']

# fileList = ['Akshan_240927', 'Akshan_241016', 'Akshan_241017', 'Akshan_241025']


################################################# Both #################################################################
# # # for loop to run through all files
# fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
#             'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
#             'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
#             'Meetz_230126', 'Akshan_240826', 'Akshan_240927', 'Akshan_240930',
#             'Akshan_241002', 'Akshan_241016', 'Akshan_241017', 'Akshan_241021',
#             'Akshan_241023', 'Akshan_241025']

# # for loop to run through all files but excluding some of akshans
fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
            'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
            'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
            'Meetz_230126', 'Akshan_240927', 'Akshan_241016', 'Akshan_241017',
            'Akshan_241025']

t0 = time.time()

popPrefFR = []
gaborSigmaSep = []
electrodeArr = np.array(np.arange(0, 32)).reshape(16, 2)
pairDistance = []
masterGoodUnits = []
masterAllUnits = []
corr_records_RFWeight = []
allUnitsLocationTuning = []


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
        masterAllUnits.append(f'{seshDate}_{unit}')

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
    # spikeCountMat = np.zeros((len(units), blocksDone+1, 49))
    spikeCountMat = np.full((len(units), blocksDone + 1, 49), np.nan)
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
    meanSpike = np.nanmean(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSD = np.nanstd(spikeCountMat[:, :blocksDone, :], axis=1)
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
        baselineResp = meanSpike[unitCount][48]
        combResp = (loc0Resp + loc1Resp) / 2
        extRespMat = np.concatenate((combResp[3:], combResp, combResp[:3]), axis=0)
        maxIndex = np.where(combResp == np.max(combResp))[0][0] + 3
        x = angleMat[maxIndex-3:maxIndex+4]
        y = extRespMat[maxIndex-3:maxIndex+4]
        params = gaussFit(x, y)
        yPred = gauss(x, *params)
        r2 = r2_score(y, yPred)
        # if r2 > 0.90:   # working version
        if r2 > 0.60:
            filterUnits.append(unit)
        unitGaussMean[unitCount] = params[2] % 360
        unitGaussSig[unitCount] = params[3]

        # fit von Mises for each location and save to allUnitsLocationTuning ---
        angles_loc0 = np.arange(0, 360, 60)  # degrees for loc 0
        angles_loc1 = np.arange(0, 360, 60)  # degrees for loc 1

        vm_loc0 = fit_von_mises_params_fixed_baseline(angles_loc0, loc0Resp, baselineResp)
        vm_loc1 = fit_von_mises_params_fixed_baseline(angles_loc1, loc1Resp, baselineResp)

        allUnitsLocationTuning.append({
            "unit_id": f"{seshDate}_{unit}",
            "monkey": monkeyName,
            "loc0_params": vm_loc0,  # dict: A, mu_deg, kappa, B
            "loc1_params": vm_loc1,  # dict: A, mu_deg, kappa, B
        })
    ####################################################################################################################
    ######################################## get units whose preferred response is #####################################
    ############################### statistically reliable above baseline in both locations ############################
    ####################################################################################################################
    def safe_mwu(a, b):
        a_clean = a[np.isfinite(a)]
        b_clean = b[np.isfinite(b)]
        if len(a_clean) < 3 or len(b_clean) < 3:
            return np.nan, np.nan
        return mannwhitneyu(a_clean, b_clean, alternative='greater')

    criterionPassedUnits = []
    # for unit in filterUnits:
    #     unitCount = np.where(unit == units)[0][0]
    #     # find direction tested that is closest to the pref dir
    #     # and reindex around this so that it is in the middle of the grid
    #     prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
    #
    #     nullDirIndex = np.where(dirArray == nullDir)[0][0]
    #     reIndex = (np.array([0, 1, 2, 3, 4, 5]) + nullDirIndex) % 6
    #
    #     # matrix of indices for b and bReIndexed
    #     stimMat, stimMatReIndex = reIndexedStimMat(reIndex)
    #     prefLoc1Indx = int(stimMatReIndex[3, 6])
    #     npLoc1Indx = int(stimMatReIndex[0, 6])
    #     prefLoc0Indx = int(stimMatReIndex[6, 3])
    #     npLoc0Indx = int(stimMatReIndex[6, 0])
    #
    #     sponSpikes = spikeCountMat[unitCount][:blocksDone, 48]
    #     prefLoc0Spikes = spikeCountMat[unitCount, :blocksDone, prefLoc0Indx]
    #     npLoc0Spikes = spikeCountMat[unitCount, :blocksDone, npLoc0Indx]
    #     prefLoc1Spikes = spikeCountMat[unitCount, :blocksDone, prefLoc1Indx]
    #     npLoc1Spikes = spikeCountMat[unitCount, :blocksDone, npLoc1Indx]
    #
    #     stat0, p_value0 = mannwhitneyu(prefLoc0Spikes, sponSpikes,
    #                                    alternative='greater')
    #     stat1, p_value1 = mannwhitneyu(prefLoc1Spikes, sponSpikes,
    #                                    alternative='greater')
    #     stat2, p_value2 = mannwhitneyu(prefLoc0Spikes, npLoc0Spikes,
    #                                    alternative='greater')
    #     stat3, p_value3 = mannwhitneyu(prefLoc1Spikes, npLoc1Spikes,
    #                                    alternative='greater')
    #     if p_value0 < 0.05 and p_value1 < 0.05 and p_value2 < 0.05 and p_value3 < 0.05:
    #         criterionPassedUnits.append(unit)
    #     # if p_value0 < 0.05 and p_value1 < 0.05:
    #     #     criterionPassedUnits.append(unit)
    for unit in filterUnits:
        unitCount = np.where(unit == units)[0][0]
        prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
        nullDirIndex = np.where(dirArray == nullDir)[0][0]
        reIndex = (np.array([0, 1, 2, 3, 4, 5]) + nullDirIndex) % 6

        stimMat, stimMatReIndex = reIndexedStimMat(reIndex)
        prefLoc1Indx = int(stimMatReIndex[3, 6])
        npLoc1Indx = int(stimMatReIndex[0, 6])
        prefLoc0Indx = int(stimMatReIndex[6, 3])
        npLoc0Indx = int(stimMatReIndex[6, 0])

        sponSpikes = spikeCountMat[unitCount][:blocksDone, 48]
        prefLoc0Spikes = spikeCountMat[unitCount, :blocksDone, prefLoc0Indx]
        npLoc0Spikes = spikeCountMat[unitCount, :blocksDone, npLoc0Indx]
        prefLoc1Spikes = spikeCountMat[unitCount, :blocksDone, prefLoc1Indx]
        npLoc1Spikes = spikeCountMat[unitCount, :blocksDone, npLoc1Indx]

        stat0, p_value0 = safe_mwu(prefLoc0Spikes, sponSpikes)
        stat1, p_value1 = safe_mwu(prefLoc1Spikes, sponSpikes)
        stat2, p_value2 = safe_mwu(prefLoc0Spikes, npLoc0Spikes)
        stat3, p_value3 = safe_mwu(prefLoc1Spikes, npLoc1Spikes)

        if all(p < 0.05 for p in [p_value0, p_value1, p_value2, p_value3]):
            criterionPassedUnits.append(unit)
            masterGoodUnits.append(f'{seshDate}_{unit}')

    ####################################################################################################################
    ############################################# different combinations of units ######################################
    ####################################################################################################################

    # combinations of units that pass inclusion criteria
    criterionPassedCombs = [i for i in combinations(criterionPassedUnits, 2)]
    # combinations of units separated by more than 250ums
    distCombs = []
    for i in criterionPassedCombs:
        n1 = unitsChannel[np.where(units == i[0])[0][0]]
        n2 = unitsChannel[np.where(units == i[1])[0][0]]
        n1ElectrodePos = np.where(electrodeArr == n1)[0][0]
        n2ElectrodePos = np.where(electrodeArr == n2)[0][0]
        pairDistOnElectrode = abs(n1ElectrodePos - n2ElectrodePos)
        pairDistance.append(pairDistOnElectrode * 50)
        # if pairDistOnElectrode >= 5:  # working version
        if pairDistOnElectrode >= 5:
            distCombs.append(i)

    # population analysis for reduced stim set (how to fit normalziation equations)
    # split the 6x6 matrix into all possible 2x2 matrices and run the normalization fit on
    # all the sub-matrices (technically 7x7 matrix with blank condition)
    rfWeightR2 = []
    bramR2 = []
    # initialize lists for paired Gabors
    pairPairedCorr = []
    pairSelectivityIndex = []
    pairNonPrefSuppIndex = []
    # initialize lists for single Gabors
    pairSingleCorr = []
    pairSingleSelectivityIndex = []
    pairSingleNonPrefSuppIndex = []

    pairs_180_apart = [pair for pair in combinations(range(6), 2) if (abs(pair[0] - pair[1]) % 6 == 3)]
    subsampleMatrix = [i for i in combinations([0, 1, 2, 3, 4, 5], 2)]
    for subpair in subsampleMatrix:
    # for subpair in pairs_180_apart:
        skipSubpair = False
        unitNormFitEstimate = [[] for i in range(len(units))]
        sli = [subpair[0], subpair[1], 6]
        dir1 = dirArray[sli[0]]
        dir2 = dirArray[sli[1]]
        condArr = np.array([dir1, dir2, dir1, dir2])

        for unit in units:
            unitCount = np.where(unit == units)[0][0]
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

            fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9),
                                                    stimIndexDict, dir1, dir2)
            resp = bCondensed.reshape(9)
            # rf weight normalization
            # parameters for fitting
            resp_loc0 = np.mean([bCondensed[2, 0], bCondensed[2, 1]])
            resp_loc1 = np.mean([bCondensed[0, 2], bCondensed[1, 2]])
            eps = 1e-6
            w1_guess = resp_loc1/(resp_loc0 + eps)  # epsilon to prevent guess from being division by 0
            log_w1_guess = np.log(w1_guess)

            # # with baseline
            guess0 = np.concatenate((bCondensed[2, :-1]+eps,
                                     [log_w1_guess],
                                     [0.1], [bCondensed[2, 2]+eps]), axis=0)
            lb = [0, 0, np.log(0.01), 0, 0]
            ub = [np.inf, np.inf, np.log(100), np.inf, np.inf]

            # # # with sigma and no baseline
            # guess0 = np.concatenate((bCondensed[2, :-1]+eps,
            #                          [log_w1_guess],
            #                          [0.1]), axis=0)
            # lb = [0, 0, np.log(0.01), 0]
            # ub = [np.inf, np.inf, np.log(100), np.inf]

            # # without sigma and no baseline
            # guess0 = np.concatenate((bCondensed[2, :-1]+eps,
            #                          [log_w1_guess]), axis=0)
            # lb = [0, 0, np.log(0.01)]
            # ub = [150, 150, np.log(100)]

            # pOpt, pCov = multi_start_curve_fit(rfWeightCondensed, fixedVals, resp.squeeze(),
            #                                    n_starts=10, seed=42)
            try:
                # Fit model
                pOpt, pCov = curve_fit(rfWeightCondensed, fixedVals, resp.squeeze(),
                                       maxfev=10000000, p0=guess0, bounds=[lb, ub])

                y_pred = rfWeightCondensed(fixedVals, *pOpt)
                r2 = r2_score(resp.squeeze(), y_pred)

                # w1 = np.exp(pOpt[2])  # Since you're fitting log(w1)
                # if w1 < 0.001 or w1 > 20:
                #     skipSubpair = True
                #     break

                # Append fit results
                rfWeightR2.append(r2)
                unitNormFitEstimate[unitCount] = pOpt

            except (RuntimeError, ValueError) as e:
                print(f"Fit failed for unit {unit} in subpair {subpair}: {e}, session, {file}")
                skipSubpair = True
                break  # break out of unit loop

        if skipSubpair:
            continue  # skip rest of subpair loop,

        # for pairCount, pair in enumerate(distCombs):
        for pairCount, pair in enumerate(criterionPassedCombs):
            n1 = np.where(units == pair[0])[0][0]
            n2 = np.where(units == pair[1])[0][0]

            u1 = f"{seshDate}_{pair[0]}"
            u2 = f"{seshDate}_{pair[1]}"

            # find predicted responses for the 3x3 matrix for each neuron:
            fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9),
                                                    stimIndexDict, dir1, dir2)
            # Neuron 1
            y_predN1 = rfWeightCondensed(fixedVals, *unitNormFitEstimate[n1]).reshape(3, 3)
            # Neuron 2
            y_predN2 = rfWeightCondensed(fixedVals, *unitNormFitEstimate[n2]).reshape(3, 3)

            # indices and correlations for paired Gabor stimuli
            for count, i in enumerate(np.int_(stimMatCond[:2, :2].reshape(4))):
                # correlation for that Gabor pair b/w 2 units excluding trials where
                # spike counts exceeded 3 SD from mean
                n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
                pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat,
                                                                     method='shrinkage')

                # Breakpoint logic: skip this count if correlation is not computable
                if np.isnan(pairStimCorr):  # or abs(pairStimCorr) > 0.4:
                    continue  # skip to next i in the loop

                # extract directions of the gabor pair
                loc0Dir = stimIndexDict[i][0]['direction']
                loc1Dir = stimIndexDict[i][1]['direction']

                # extract row, col of paired condition
                row = count // 2
                col = count % 2

                # bram selectivity
                # n1 selectivity and suppression index
                loc0Resp = y_predN1[2, col]
                loc1Resp = y_predN1[row, 2]
                # loc0Resp = unitNormFitEstimate[n1][np.where(condArr == loc0Dir)[0][0]]
                # loc1Resp = unitNormFitEstimate[n1][np.where(condArr == loc1Dir)[0][0]] * np.exp(
                #     unitNormFitEstimate[n1][2])
                w0 = 1
                w1 = np.exp(unitNormFitEstimate[n1][2])
                sig1 = unitNormFitEstimate[n1][3]
                n1Selectivity, n1NonPrefSupp = getSelAndSuppIndxWithSig(loc0Resp, loc1Resp, w0, w1, sig1)

                # n2 selectivity and suppression index
                loc0Resp = y_predN2[2, col]
                loc1Resp = y_predN2[row, 2]
                # loc0Resp = unitNormFitEstimate[n2][np.where(condArr == loc0Dir)[0][0]]
                # loc1Resp = unitNormFitEstimate[n2][np.where(condArr == loc1Dir)[0][0]] * np.exp(
                #     unitNormFitEstimate[n2][2])
                w0 = 1
                w1 = np.exp(unitNormFitEstimate[n2][2])
                sig2 = unitNormFitEstimate[n2][3]
                n2Selectivity, n2NonPrefSupp = getSelAndSuppIndxWithSig(loc0Resp, loc1Resp, w0, w1, sig2)

                # pair selectivity and suppression index
                pairSelectivity = (np.sign(n1Selectivity) * np.sign(n2Selectivity) *
                                   np.sqrt(abs(n1Selectivity) * abs(n2Selectivity)))
                # eps = 1e-6  # to prevent atanh(±1)
                # # Clip values to avoid infinite z-transform
                # n1_clipped = np.clip(abs(n1Selectivity), -1 + eps, 1 - eps)
                # n2_clipped = np.clip(abs(n2Selectivity), -1 + eps, 1 - eps)
                # # Fisher z-transform (atanh), average, then back-transform (tanh)
                # z1 = np.arctanh(n1_clipped)
                # z2 = np.arctanh(n2_clipped)
                # z_mean = 0.5 * (z1 + z2)
                # pairSelectivity = np.sign(n1Selectivity) * np.sign(n2Selectivity) * np.tanh(z_mean)
                pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)

                # append to list
                pairPairedCorr.append(pairStimCorr)
                pairSelectivityIndex.append(pairSelectivity)
                pairNonPrefSuppIndex.append(pairSuppression)

                # loc 0 single Gabor corr, excluding trials where spike counts > 3 SD
                stimIndex = np.int_(stimMatCond[2, :][np.where(condArr == loc0Dir)[0][0]])

                n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
                singleStimLoc0Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat,
                                                                           method='shrinkage')
                # # n1
                # w0 = 1
                # sig = unitNormFitEstimate[n1][3]
                # n1supp = w0/(w0+sig)
                # # n2
                # w0 = 1
                # sig = unitNormFitEstimate[n2][3]
                # n2supp = w0/(w0+sig)
                #
                # pairSuppression = np.sqrt(n1supp * n2supp)

                pairSingleCorr.append(singleStimLoc0Corr)
                pairSingleSelectivityIndex.append(pairSelectivity)
                pairSingleNonPrefSuppIndex.append(pairSuppression)

                # loc 1 single Gabor corr, excluding trials where spike counts > 3 SD
                stimIndex = np.int_(stimMatCond[:, 2][np.where(condArr == loc1Dir)[0][0]])

                n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
                singleStimLoc1Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat,
                                                                           method='shrinkage')

                # # n1
                # w1 = np.exp(unitNormFitEstimate[n1][2])
                # sig = unitNormFitEstimate[n1][3]
                # n1supp = w1 / (w1 + sig)
                # # n2
                # w1 = np.exp(unitNormFitEstimate[n2][2])
                # sig = unitNormFitEstimate[n2][3]
                # n2supp = w1 / (w1 + sig)
                #
                # pairSuppression = np.sqrt(n1supp * n2supp)

                pairSingleCorr.append(singleStimLoc1Corr)
                pairSingleSelectivityIndex.append(pairSelectivity)
                pairSingleNonPrefSuppIndex.append(pairSuppression)

                corr_records_RFWeight.append({
                    "date": seshDate,
                    "monkey": monkeyName,
                    "unit_i": u1,
                    "unit_j": u2,
                    "pair_key": tuple(sorted((u1, u2))),  # handy for grouping/joining later
                    "cond_index": int(i),  # stimulus condition index you used
                    "cond_order": int(count),  # 0..3 within the 2x2 paired block (optional)
                    "loc0_dir": int(loc0Dir),
                    "loc1_dir": int(loc1Dir),

                    # correlations
                    "corr_paired": float(pairStimCorr),
                    "corr_single_loc0": float(singleStimLoc0Corr),
                    "corr_single_loc1": float(singleStimLoc1Corr),

                    # pair metrics
                    "neuron 1 preferred stim": n1Selectivity,
                    "neuron 2 preferred stim": n2Selectivity,
                    "neuron 1 nonpref suppression": n1NonPrefSupp,
                    "neuron 2 nonpref suppression": n2NonPrefSupp,
                    "neuron 1 w0": 1,
                    "neuron 2 w0": 1,
                    "neuron 1 w1": np.exp(unitNormFitEstimate[n1][2]),
                    "neuron 2 w1": np.exp(unitNormFitEstimate[n2][2]),
                    "neuron 1 sig": unitNormFitEstimate[n1][3],
                    "neuron 2 sig": unitNormFitEstimate[n2][3],
                    "pair_selectivity": float(pairSelectivity),
                    "pair_suppression": float(pairSuppression),
                })

    pairPairedCorr = np.array(pairPairedCorr)
    pairSelectivityIndex = np.array(pairSelectivityIndex)
    pairNonPrefSuppIndex = np.array(pairNonPrefSuppIndex)

    pairSingleCorr = np.array(pairSingleCorr)
    pairSingleSelectivityIndex = np.array(pairSingleSelectivityIndex)
    pairSingleNonPrefSuppIndex = np.array(pairSingleNonPrefSuppIndex)

    # compile similarity scores and correlation matrix into one matrix
    gaborPairCompiledArr = np.array([pairPairedCorr,
                                     pairSelectivityIndex,
                                     pairNonPrefSuppIndex])
    gaborSingleCompiledArr = np.array([pairSingleCorr,
                                       pairSingleSelectivityIndex,
                                       pairSingleNonPrefSuppIndex])
    np.save(f'../../condRFWSingleCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborSingleCompiledArr)
    np.save(f'../../condRFWPairedCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborPairCompiledArr)
    os.chdir('../../../Python Code')

print(time.time()-t0)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # plot # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import pandas as pd

eps = 1e-9  # numerical guard

def delta_angle(deg1, deg2):
    """Return absolute difference between angles in degrees, folded to [0,180]."""
    diff = abs(deg1 - deg2) % 360
    if diff > 180:
        diff = 360 - diff
    return diff

# Map: unit_id -> avg mu_deg (from your location tuning struct/list)
mu_map = {}
for u in allUnitsLocationTuning:
    unit_id = u["unit_id"]
    mu0 = u["loc0_params"]["mu_deg"]
    mu1 = u["loc1_params"]["mu_deg"]
    mu_map[unit_id] = np.mean([mu0, mu1])  # average preferred direction across locations

new_records = []
for rec in corr_records_RFWeight:
    # Ensure we work on a copy
    rec = dict(rec)

    u1, u2 = rec["unit_i"], rec["unit_j"]

    # --- preferred direction difference ---
    mu1 = mu_map.get(u1, np.nan)
    mu2 = mu_map.get(u2, np.nan)
    delta_mu = delta_angle(mu1, mu2) if np.isfinite(mu1) and np.isfinite(mu2) else np.nan

    # --- shares computed from fields in corr_records_RFWeight ---
    # neuron 1
    w0_1  = float(rec.get("neuron 1 w0", np.nan))
    w1_1  = float(rec.get("neuron 1 w1", np.nan))
    sig_1 = float(rec.get("neuron 1 sig", np.nan))

    # neuron 2
    w0_2  = float(rec.get("neuron 2 w0", np.nan))
    w1_2  = float(rec.get("neuron 2 w1", np.nan))
    sig_2 = float(rec.get("neuron 2 sig", np.nan))

    # s = w / (w0 + w + sigma), with guards for NaNs
    def safe_share(w0, w, sig):
        if not (np.isfinite(w0) and np.isfinite(w) and np.isfinite(sig)):
            return np.nan
        denom = w0 + w + sig + eps
        if denom <= 0:
            return np.nan
        return w / denom

    s1 = safe_share(w0_1, w1_1, sig_1)
    s2 = safe_share(w0_2, w1_2, sig_2)

    # Geometric mean for the pair
    if np.isfinite(s1) and np.isfinite(s2) and s1 >= 0 and s2 >= 0:
        pair_share = np.sqrt(s1 * s2)
    else:
        pair_share = np.nan

    # --- write out fields ---
    rec["delta_pref_dir"]   = delta_mu
    rec["unit1_avg_weight"] = s1
    rec["unit2_avg_weight"] = s2
    rec["pair_avg_weight"]  = pair_share

    new_records.append(rec)

df_corr_ext = pd.DataFrame(new_records)

# Quick check of the new columns
print(
    df_corr_ext[
        ["unit_i", "unit_j", "delta_pref_dir",
         "unit1_avg_weight", "unit2_avg_weight", "pair_avg_weight"]
    ].head()
)

for filler in range(1):
    savemat("corr_data.mat", {
        "nonPrefSupp": df_corr_ext['pair_avg_weight'].values,
        # "selectivity": df_corr_ext['delta_pref_dir'].values,
        "selectivity": df_corr_ext['pair_selectivity'].values,
        "pairedCorr": df_corr_ext['corr_paired'].values,
        "singleCorrLoc0": df_corr_ext['corr_single_loc0'].values,
        "singleCorrLoc1": df_corr_ext['corr_single_loc1'].values,
        "monkey": df_corr_ext['monkey'].values.astype('U')  # save as strings
    })

    x = df_corr_ext['pair_avg_weight'].values
    y = df_corr_ext['delta_pref_dir'].values
    c = df_corr_ext['corr_paired'].values

    plt.figure(figsize=(7, 6))
    sc = plt.scatter(
        x, y,
        c=c,
        cmap='bwr',  # red = positive, blue = negative, white ~ 0
        vmin=-1, vmax=1,  # fix scale from -1 to 1
        alpha=0.7,
        edgecolor='k', linewidth=0.3
    )
    plt.colorbar(sc, label='Paired correlation')
    plt.xlabel('Pair geometric mean suppression')
    plt.ylabel('Pair selectivity')
    plt.title('Pair selectivity vs suppression\nColored by correlation')
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.show()



