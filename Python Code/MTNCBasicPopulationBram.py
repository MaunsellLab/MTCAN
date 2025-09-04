# imports
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

# fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
#             'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
#             'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
#             'Meetz_230126']

########################################### AKSHAN ###########################################

# good sessions: 240927, 240930, 241016, 241017, 241023, 241025
# okay sessions: 240826, 241002 (13 blocks), 241021 (17 blocks)
# bad sessions: 240827, 240828

# fileList = ['Akshan_240826', 'Akshan_240927', 'Akshan_240930', 'Akshan_241002',
#             'Akshan_241016', 'Akshan_241017', 'Akshan_241021', 'Akshan_241023',
#             'Akshan_241025']
# fileList = ['Akshan_240927', 'Akshan_241016', 'Akshan_241017', 'Akshan_241025']


########################################### Both ###########################################
# # # for loop to run through all files
# fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
#             'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
#             'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
#             'Meetz_230126', 'Akshan_240826', 'Akshan_240927', 'Akshan_240930',
#             'Akshan_241002', 'Akshan_241016', 'Akshan_241017', 'Akshan_241021',
#             'Akshan_241023', 'Akshan_241025']

# some of akshan's sessions removed
fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
            'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
            'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
            'Meetz_230126', 'Akshan_240927', 'Akshan_241016', 'Akshan_241017',
            'Akshan_241025']


t0 = time.time()
meetzGlobalNMI = []
meetzGlobalA1 = []
akshanGlobalNMI = []
akshanGlobalA1 = []
electrodeArr = np.array(np.arange(0, 32)).reshape(16, 2)
masterGoodUnits = []
masterAllUnits = []
corr_records = []
allUnitsLocationTuning = []

for fileIterator in fileList:
    print(fileIterator)
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
    if 'blockStatus' in allTrials[corrTrials[-1]]:
        blocksDone = allTrials[corrTrials[-1]]['blockStatus']['data']['blocksDone']
    else:
        blocksDone = allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone']
    highContrast, zeroContrast = max(stimIndexDF['loc0 Contrast'].unique()), \
                                 min(stimIndexDF['loc0 Contrast'].unique())
    zeroDir = 0
    dirArray = np.array([0, 60, 120, 180, 240, 300])
    loc0IndexArray = np.array([36, 37, 38, 39, 40, 41])
    loc1IndexArray = np.array([42, 43, 44, 45, 46, 47])
    angleMat = np.arange(180, 900, 60)
    spikeCountMat = np.full((len(units), blocksDone+1, 49), np.nan)
    spikeCountLong = []
    onLatency = 50/1000  # time in MS for counting window latency after stim on
    offLatency = 50/1000  # time in MS for counting window latency after stim off
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

        # --- NEW: fit von Mises for each location and save to allUnitsLocationTuning ---
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
    def safe_mwu(a, _):
        a_clean = a[np.isfinite(a)]
        b_clean = _[np.isfinite(_)]
        if len(a_clean) < 3 or len(b_clean) < 3:
            return np.nan, np.nan
        return mannwhitneyu(a_clean, b_clean, alternative='greater')

    criterionPassedUnits = []
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
    ################################################# Different Combs ##################################################
    ####################################################################################################################
    # the different combinations of neuron pairs from total units
    combs = [i for i in combinations(units, 2)]
    # combinations of units that pass inclusion criteria
    criterionPassedCombs = [i for i in combinations(criterionPassedUnits, 2)]
    # combinations of units separated by more than 250ums
    distCombs = []
    # for i in combs:  #### working version
    for i in criterionPassedCombs:
        n1 = unitsChannel[np.where(units == i[0])[0][0]]
        n2 = unitsChannel[np.where(units == i[1])[0][0]]
        n1ElectrodePos = np.where(electrodeArr == n1)[0][0]
        n2ElectrodePos = np.where(electrodeArr == n2)[0][0]
        pairDistOnElectrode = abs(n1ElectrodePos - n2ElectrodePos)
        if pairDistOnElectrode >= 5:
            distCombs.append(i)

    ####################################################################################################################
    ############################################## Compute Correlations ################################################
    ####################################################################################################################
    # initialize lists for paired Gabors
    pairPairedCorr = []
    pairSelectivityIndex = []
    pairNonPrefSuppIndex = []
    # initialize lists for single Gabor presentations
    pairSingleCorr = []
    pairSingleSelectivityIndex = []
    pairSingleNonPrefSuppIndex = []
    # initialize lists for blank Gabor presentations
    pairBlankCorr = []

    # # subsampling the 6x6 matrix of directions into unique 2x2 matrices
    subsampleMatrix = [i for i in combinations([0, 1, 2, 3, 4, 5], 2)]
    for subpair in subsampleMatrix:
        skipSubPair = False
        unitPairedNormFit = []
        unitPairedNormR2 = np.zeros(len(units))
        unitNormFitEstimate = [[] for i in range(len(units))]
        sli = [subpair[0], subpair[1], 6]
        dir1 = dirArray[sli[0]]
        dir2 = dirArray[sli[1]]
        condArr = np.array([dir1, dir2, dir1, dir2])

        for unitCount, unit in enumerate(units):
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

            # (working version has initial guess for alpha = 0.2)
            # Generic Normalization (L1+L2)/(1+al2) w.o scalar
            guess0 = np.concatenate((bCondensed[2, :-1],
                                     bCondensed[:-1, 2],
                                     [0.2]), axis=0)

            resp = bCondensed.reshape(9)[:-1]  # all responses except for baseline
            fixedVals = fixedValsForEMSGenCondensed(stimMatCond.reshape(9)[:-1],
                                                    stimIndexDict, dir1, dir2)

            try:
                # fit model
                pOpt, pCov, = curve_fit(genNormCondensed, fixedVals, resp.squeeze(),
                                        maxfev=10000000)

                y_pred = genNormCondensed(fixedVals, *pOpt)
                r2 = r2_score(resp.squeeze(), y_pred)

                # bram trick to keep values > 0
                pOpt = pOpt ** 2

                # Append fit parameters for condensed matrix
                unitPairedNormR2[unitCount] = r2
                unitPairedNormFit.append(pOpt)
                unitNormFitEstimate[unitCount] = pOpt
                if monkeyName == 'Meetz':
                    meetzGlobalA1.append(pOpt[4])
                else:
                    akshanGlobalA1.append(pOpt[4])

            except (RuntimeError, ValueError) as e:
                print(f'Fit failed for unit {unit} in subpair {subpair}: {e}, session, {fileIterator}')
                skipSubPair = True
                break

        if skipSubPair:
            continue

        # generate paired stimulus correlation, selectivity index,
        # suppression, and NMI index for each gabor pair for each pair
        # of neurons (similar to Bram fig 2c)
        tempMat = stats.zscore(spikeCountMat[:, :blocksDone, :], axis=1, nan_policy='omit')
        stimIndexMat = np.arange(36).reshape(6, 6)
        upperTriangle = upperTriMasking(stimIndexMat)

        for pairCount, pair in enumerate(criterionPassedCombs):
        # for pairCount, pair in enumerate(distCombs):  # working version
            n1 = np.where(units == pair[0])[0][0]
            n2 = np.where(units == pair[1])[0][0]

            u1 = f"{seshDate}_{pair[0]}"
            u2 = f"{seshDate}_{pair[1]}"

            # unit's preferred directions
            n1PrefDir = unitGaussMean[n1]
            n2PrefDir = unitGaussMean[n2]

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
                pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)

                pairPairedCorr.append(pairStimCorr)
                pairSelectivityIndex.append(pairSelectivity)
                pairNonPrefSuppIndex.append(pairSuppression)

                # loc 0 single Gabor corr, excluding trials where spike counts > 3 SD
                stimIndex = np.int_(stimMatCond[2, :][np.where(condArr == loc0Dir)[0][0]])
                n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
                singleStimLoc0Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                pairSingleCorr.append(singleStimLoc0Corr)
                pairSingleSelectivityIndex.append(pairSelectivity)
                pairSingleNonPrefSuppIndex.append(pairSuppression)

                # loc 1 single Gabor corr, excluding trials where spike counts > 3 SD
                stimIndex = np.int_(stimMatCond[:, 2][np.where(condArr == loc1Dir)[0][0]])
                n1SpikeMat = spikeCountMat[n1, :blocksDone, stimIndex]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, stimIndex]
                singleStimLoc1Corr, pairSCov, pairSSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                pairSingleCorr.append(singleStimLoc1Corr)
                pairSingleSelectivityIndex.append(pairSelectivity)
                pairSingleNonPrefSuppIndex.append(pairSuppression)

                # pair blank corr
                blankIndex = 48
                skipTrials = []
                n1SpikeMat = spikeCountMat[n1, :blocksDone, blankIndex]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, blankIndex]
                blankCorr, blankCov, blankSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                pairBlankCorr.append(blankCorr)

                corr_records.append({
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
                    "neuron 1 a0": 1,
                    "neuron 2 a0": 1,
                    "neuron 1 a1": unitNormFitEstimate[n1][4],
                    "neuron 2 a1": unitNormFitEstimate[n2][4],
                    "pair_selectivity": float(pairSelectivity),
                    "pair_suppression": float(pairSuppression),

                    # optional provenance / extras:
                    "n1_pref_dir": float(n1PrefDir),
                    "n2_pref_dir": float(n2PrefDir),
                })

    ####################################################################################################################
    ################################################### Convert to Arrays ##############################################
    ################################################ append to master list #############################################
    ####################################################################################################################
    pairPairedCorr = np.array(pairPairedCorr)
    pairSelectivityIndex = np.array(pairSelectivityIndex)
    pairNonPrefSuppIndex = np.array(pairNonPrefSuppIndex)
    pairSingleCorr = np.array(pairSingleCorr)
    pairSingleSelectivityIndex = np.array(pairSingleSelectivityIndex)
    pairSingleNonPrefSuppIndex = np.array(pairSingleNonPrefSuppIndex)
    pairBlankCorr = np.array(pairBlankCorr)

    # compile similarity scores and correlation matrix into one matrix
    gaborPairCompiledArr = np.array([pairPairedCorr,
                                     pairSelectivityIndex,
                                     pairNonPrefSuppIndex,
                                     pairBlankCorr])
    gaborSingleCompiledArr = np.array([pairSingleCorr,
                                       pairSingleSelectivityIndex,
                                       pairSingleNonPrefSuppIndex])
    np.save(f'../../gaborSingleCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborSingleCompiledArr)
    np.save(f'../../gaborPairCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborPairCompiledArr)
    os.chdir('../../../Python Code')

print(time.time()-t0)


########################################################################################################################
######################################## Assign relative weights to each location ######################################
################################################### from RF Tuning data ################################################
########################################################################################################################


def gaussian2d(coords, amp, x0, y0, sx, sy, rho, baseline):
    X, Y = coords
    sx = np.maximum(sx, 1e-6)
    sy = np.maximum(sy, 1e-6)
    rho = np.clip(rho, -0.95, 0.95)
    Xn = (X - x0) / sx
    Yn = (Y - y0) / sy
    z = (Xn**2 - 2*rho*Xn*Yn + Yn**2) / (2*(1 - rho**2))
    return baseline + amp * np.exp(-z)


def fit_gaussian2d_fixed_baseline(Z, azis, eles, baseline_fixed):
    """
    Fit 2D Gaussian with correlation (rho) BUT with baseline fixed to 'baseline_fixed'.
    Returns (params_dict, Zhat) or (None, None) on failure.
    """
    if not np.isfinite(baseline_fixed):
        return None, None

    A, E = np.meshgrid(azis, eles, indexing='ij')
    X = A.ravel(); Y = E.ravel(); y = Z.ravel()
    mask = np.isfinite(y)
    if mask.sum() < 7:
        return None, None

    Xv, Yv, yv = X[mask], Y[mask], y[mask]

    # Initial guesses: amplitude relative to fixed baseline
    amp0 = float(max(np.nanmax(y) - baseline_fixed, 1e-6)) if np.isfinite(np.nanmax(y)) else 1.0

    # Weighted COM (above baseline) for center
    w = np.clip(y - baseline_fixed, 0, None); w[~np.isfinite(w)] = 0
    if w.sum() > 0:
        x0 = float(np.sum(A.ravel() * w) / w.sum())
        y0 = float(np.sum(E.ravel() * w) / w.sum())
    else:
        x0 = float(np.median(azis)); y0 = float(np.median(eles))

    sx0 = max(np.std(azis) / 2, 1e-3)
    sy0 = max(np.std(eles) / 2, 1e-3)
    rho0 = 0.0

    # Params to fit: amp, x0, y0, sx, sy, rho  (baseline fixed)
    p0 = [amp0, x0, y0, sx0, sy0, rho0]
    lb = [0.0,  min(azis)-10, min(eles)-10, 1e-3, 1e-3, -0.95]
    ub = [np.inf, max(azis)+10, max(eles)+10, np.ptp(azis)+50, np.ptp(eles)+50, 0.95]

    # Model wrapper with fixed baseline
    def _g2d_fixedB(coords, amp, x0, y0, sx, sy, rho):
        Xc, Yc = coords
        sx_ = np.maximum(sx, 1e-6)
        sy_ = np.maximum(sy, 1e-6)
        rho_ = np.clip(rho, -0.95, 0.95)
        Xn = (Xc - x0) / sx_
        Yn = (Yc - y0) / sy_
        z = (Xn**2 - 2*rho_*Xn*Yn + Yn**2) / (2*(1 - rho_**2))
        return baseline_fixed + amp * np.exp(-z)

    try:
        popt, _ = curve_fit(_g2d_fixedB, np.vstack([Xv, Yv]), yv,
                            p0=p0, bounds=(lb, ub), maxfev=20000)
        # success with rho
        amp, x0, y0, sx, sy, rho = popt
    except Exception:
        # Fallback: axis-aligned (rho = 0) with fixed baseline
        def _g2d_axis_fixedB(coords, amp, x0, y0, sx, sy):
            Xc, Yc = coords
            sx_ = np.maximum(sx, 1e-6)
            sy_ = np.maximum(sy, 1e-6)
            z = ((Xc - x0)**2) / (2*sx_**2) + ((Yc - y0)**2) / (2*sy_**2)
            return baseline_fixed + amp * np.exp(-z)

        p0b = [amp0, x0, y0, sx0, sy0]
        lbb = [0.0,  min(azis)-10, min(eles)-10, 1e-3, 1e-3]
        ubb = [np.inf, max(azis)+10, max(eles)+10, np.ptp(azis)+50, np.ptp(eles)+50]
        try:
            popt_b, _ = curve_fit(_g2d_axis_fixedB, np.vstack([Xv, Yv]), yv,
                                  p0=p0b, bounds=(lbb, ubb), maxfev=20000)
            amp, x0, y0, sx, sy = popt_b
            rho = 0.0
        except Exception:
            return None, None

    # Build outputs
    Zhat = (baseline_fixed + amp * np.exp(
        -(((A - x0)/max(sx,1e-6))**2
          - 2*np.clip(rho, -0.95, 0.95)*((A - x0)/max(sx,1e-6))*((E - y0)/max(sy,1e-6))
          + ((E - y0)/max(sy,1e-6))**2
        ) / (2*(1 - np.clip(rho, -0.95, 0.95)**2))
    )).astype(float)

    params = {
        "amp": float(amp), "x0": float(x0), "y0": float(y0),
        "sx": float(sx), "sy": float(sy), "rho": float(rho),
        "baseline": float(baseline_fixed)  # <- fixed to measured spontaneous rate
    }
    return params, Zhat


def eval_weights_at_stims(params, stim_locs):
    if params is None:
        return np.array([np.nan, np.nan]), np.nan
    peak = gaussian2d(
        np.vstack([np.array([params["x0"]]), np.array([params["y0"]])]),
        params["amp"], params["x0"], params["y0"], params["sx"], params["sy"], params["rho"], params["baseline"]
    )[0]
    vals = []
    for (az, el) in stim_locs:
        vals.append(gaussian2d(
            np.vstack([np.array([az]), np.array([el])]),
            params["amp"], params["x0"], params["y0"], params["sx"], params["sy"], params["rho"], params["baseline"]
        )[0])
    vals = np.array(vals, dtype=float)
    weights = np.full_like(vals, np.nan, dtype=float) if (not np.isfinite(peak) or abs(peak)<1e-9) else vals/peak
    return weights, peak


# ------------ main ------------
fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
            'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
            'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
            'Meetz_230126', 'Akshan_240927', 'Akshan_241016', 'Akshan_241017',
            'Akshan_241025']

base_dir = "/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN"

results_all = []
for date in fileList:
    monkeyName, sessionDate = date.split('_')
    sesh_dir = os.path.join(base_dir, monkeyName, sessionDate, "RFLoc Tuning")
    os.chdir(sesh_dir)

    unitsRFRespData = np.load('unitsRFLocMat.npy')   # (n_units, n_azi, n_ele)
    unitsBaselineRespData = np.load('unitsBaseline.npy')  # (n_units)
    sessionUnitIDs = np.asarray(np.load('seshUnitsID.npy', allow_pickle=True)).astype(str).tolist()
    seshAziEle = np.load('seshAziEle.npy')
    seshStimLocs = np.load('seshStimLocs.npy')

    if unitsRFRespData.shape[0] != len(sessionUnitIDs):
        print(f"[WARN] {sessionDate}: rows ({unitsRFRespData.shape[0]}) != sessionUnitIDs ({len(sessionUnitIDs)})")

    # derive grid axes
    if seshAziEle.ndim == 2 and seshAziEle.shape[1] == 2:
        azis_unique = np.unique(seshAziEle[:, 0])
        eles_unique = np.unique(seshAziEle[:, 1])
    else:
        raise ValueError(f"seshAziEle unexpected shape {seshAziEle.shape}; expected (N,2).")

    n_azi, n_ele = len(azis_unique), len(eles_unique)
    if (n_azi, n_ele) != unitsRFRespData.shape[1:]:
        if (n_ele, n_azi) == unitsRFRespData.shape[1:]:
            azis_unique, eles_unique = eles_unique, azis_unique
            n_azi, n_ele = n_ele, n_azi
        else:
            print(f"[ERROR] {sessionDate}: grid mismatch: data {unitsRFRespData.shape[1:]}, grid ({n_azi},{n_ele}). Skipping.")
            continue

    # Good units for this day (strings like '221010_77')
    goodUnits_today = [u for u in masterGoodUnits if u.startswith(sessionDate + "_")]

    for gu in goodUnits_today:
        # Resolve index in sessionUnitIDs
        try:
            _, unit_suffix = gu.split('_', 1)
        except ValueError:
            unit_suffix = gu

        matches = [i for i, s in enumerate(sessionUnitIDs) if s == gu]
        if not matches:
            matches = [i for i, s in enumerate(sessionUnitIDs)
                       if s.endswith("_" + unit_suffix) or s == unit_suffix]
        if not matches:
            continue

        uidx = matches[0]
        Z = unitsRFRespData[uidx]  # (n_azi, n_ele)

        baseline_fixed = float(unitsBaselineRespData[uidx])  # measured spontaneous rate for this unit
        params, Zhat = fit_gaussian2d_fixed_baseline(Z, azis_unique, eles_unique, baseline_fixed)
        if params is None:
            results_all.append((
                f"{sessionDate}_{unit_suffix}",  # unit_key
                sessionDate,                     # seshDate
                monkeyName,                      # monkey name
                sessionUnitIDs[uidx],            # unitID_raw
                False,                           # fit_ok
                None,                            # params
                np.array([np.nan, np.nan]),      # weights
                np.nan                           # peak_val
            ))
            continue

        weights, peak_val = eval_weights_at_stims(params, seshStimLocs)

        results_all.append((
            f"{sessionDate}_{unit_suffix}",   # unit_key
            sessionDate,                      # seshDate
            sessionUnitIDs[uidx],             # unitID_raw (as stored on disk)
            True,                             # fit_ok
            params,                           # dict of fit params
            weights.astype(float),            # (2,)
            float(peak_val)
        ))

        os.chdir(base_dir)

# If you want a NumPy "array", convert the Python list -> object array:
results_array = np.array(results_all, dtype=object)

# plot weight distributions
for filler in range(1):
    # Extract w0 and w1 from results_array
    w0_vals = []
    w1_vals = []

    for row in results_array:
        weights = row[-2]  # weights are second-to-last
        if isinstance(weights, np.ndarray) and weights.size == 2:
            w0_vals.append(weights[0])
            w1_vals.append(weights[1])

    w0_vals = np.array(w0_vals, dtype=float)
    w1_vals = np.array(w1_vals, dtype=float)

    # Plot histograms
    plt.figure(figsize=(7, 5))
    plt.hist(w0_vals, bins=30, alpha=0.6, label='w0 (stim loc 0)')
    plt.hist(w1_vals, bins=30, alpha=0.6, label='w1 (stim loc 1)')
    plt.xlabel('Weight value')
    plt.ylabel('Count')
    plt.title('Distribution of w0 and w1 across neurons')
    plt.legend()
    plt.show()


########################################################################################################################
################################### Find difference in preferred direction and avg weight ##############################
######################################### for a pair to make a heatmap of their corr####################################
########################################################################################################################
def delta_angle(deg1, deg2):
    """Return absolute difference between angles in degrees, folded to [0,180]."""
    diff = abs(deg1 - deg2) % 360
    if diff > 180:
        diff = 360 - diff
    return diff


# Map: unit_id -> avg mu_deg
mu_map = {}
for u in allUnitsLocationTuning:
    unit_id = u["unit_id"]
    mu0 = u["loc0_params"]["mu_deg"]
    mu1 = u["loc1_params"]["mu_deg"]
    mu_map[unit_id] = np.mean([mu0, mu1])  # average pref direction

# Map: unit_id -> average weight
weight_map = {}
for row in results_array:
    if isinstance(row, np.ndarray) and row.ndim == 0:
        row = row.item()
    unit_id = str(row[0])
    weights = row[-2]
    if isinstance(weights, np.ndarray) and weights.size == 2:
        weight_map[unit_id] = float(np.mean(weights))

new_records = []
for rec in corr_records:
    u1, u2 = rec["unit_i"], rec["unit_j"]

    # --- preferred direction difference ---
    mu1 = mu_map.get(u1, np.nan)
    mu2 = mu_map.get(u2, np.nan)
    delta_mu = delta_angle(mu1, mu2) if np.isfinite(mu1) and np.isfinite(mu2) else np.nan

    # --- weights ---
    w1 = weight_map.get(u1, np.nan)
    w2 = weight_map.get(u2, np.nan)
    avg_w = np.nanmean([w1, w2])

    rec = rec.copy()
    rec["delta_pref_dir"] = delta_mu
    rec["unit1_avg_weight"] = w1
    rec["unit2_avg_weight"] = w2
    rec["pair_avg_weight"] = avg_w

    new_records.append(rec)

df_corr_ext = pd.DataFrame(new_records)

print(df_corr_ext[["unit_i", "unit_j", "delta_pref_dir",
                   "unit1_avg_weight", "unit2_avg_weight", "pair_avg_weight"]].head())

for filler in range(1):
    savemat("corr_data.mat", {
        "nonPrefSupp": df_corr_ext['pair_avg_weight'].values,
        "selectivity": df_corr_ext['delta_pref_dir'].values,
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






########################################################################################################################
################################### Recompute Bram's Metrics using Weight based off of #################################
################################################## mapping for each pair ###############################################
########################################################################################################################
# Column names present in corr_records
UNIT1_COL = 'unit_i'
UNIT2_COL = 'unit_j'
SEL1_COL  = 'neuron 1 preferred stim'   # >0 => prefers site 0; <0 => prefers site 1
SEL2_COL  = 'neuron 2 preferred stim'

A1_A0_COL = 'neuron 1 a0'
A1_A1_COL = 'neuron 1 a1'
A2_A0_COL = 'neuron 2 a0'
A2_A1_COL = 'neuron 2 a1'

# Output column names
N1_SUPP_WRF_COL = 'neuron1_supp_wRF_VM'
N2_SUPP_WRF_COL = 'neuron2_supp_wRF_VM'
PAIR_SUPP_WRF_COL = 'pair_supp_wRF_VM'

EPS = 1e-12


def build_weights_map(results_array):
    """
    Map: unit_key -> np.array([w0, w1]).
    Your rows look like:
      [unit_key, seshDate, unitID_raw, fit_ok, params_dict, weights, peak_val]
    """
    wm = {}
    for row in results_array:
        if isinstance(row, np.ndarray) and row.ndim == 0:
            row = row.item()
        try:
            unit_key = str(row[0])
            fit_ok   = bool(row[3])
            weights  = row[-2]
            if fit_ok and isinstance(weights, np.ndarray) and weights.size == 2 and np.all(np.isfinite(weights)):
                wm[unit_key] = weights.astype(float)
        except Exception:
            continue
    return wm


def vm_rf_weighted_suppression(w0, w1, a0, a1, sel_sign,
                               gamma_w=1.0, gamma_a=1.0, eps=1e-12):
    """
    Verhoef & Maunsell suppression with RF weighting applied upstream.
    Two tuning knobs:
      - gamma_w: exponent on RF weights
      - gamma_a: exponent on alphas
    """
    if not (np.isfinite(w0) and np.isfinite(w1) and np.isfinite(a0) and np.isfinite(a1)):
        return np.nan

    a0_eff = (a0 ** gamma_a) * (w0 ** gamma_w)
    a1_eff = (a1 ** gamma_a) * (w1 ** gamma_w)
    denom = a0_eff + a1_eff + eps

    if sel_sign > 0:      # prefers site 0
        return a1_eff / denom
    elif sel_sign < 0:    # prefers site 1
        return a0_eff / denom
    else:
        return np.nan


def geom_mean_nonneg(a, b):
    if np.isfinite(a) and np.isfinite(b) and a >= 0 and b >= 0:
        return float(np.sqrt(a*b))
    return np.nan


weights_map = build_weights_map(results_array)

# If corr_records is not already a DataFrame:
df_corr = pd.DataFrame(corr_records)

gamma_val_for_weights = 1
gamma_val_for_alpha = 1
# Neuron 1
df_corr[N1_SUPP_WRF_COL] = [
    vm_rf_weighted_suppression(
        *weights_map.get(str(u1), (np.nan, np.nan)),     # w0, w1
        float(row[A1_A0_COL]), float(row[A1_A1_COL]),    # a0, a1
        float(row[SEL1_COL]),
        gamma_w=gamma_val_for_weights,
        gamma_a=gamma_val_for_alpha
    )
    for u1, (_, row) in zip(df_corr[UNIT1_COL], df_corr.iterrows())
]

# Neuron 2
df_corr[N2_SUPP_WRF_COL] = [
    vm_rf_weighted_suppression(
        *weights_map.get(str(u2), (np.nan, np.nan)),     # w0, w1
        float(row[A2_A0_COL]), float(row[A2_A1_COL]),    # a0, a1
        float(row[SEL2_COL]),
        gamma_w=gamma_val_for_weights,
        gamma_a=gamma_val_for_alpha
    )
    for u2, (_, row) in zip(df_corr[UNIT2_COL], df_corr.iterrows())
]

# Pair (geometric mean)
df_corr[PAIR_SUPP_WRF_COL] = [
    geom_mean_nonneg(a, b) for a, b in zip(df_corr[N1_SUPP_WRF_COL], df_corr[N2_SUPP_WRF_COL])
]

print(df_corr[[N1_SUPP_WRF_COL, N2_SUPP_WRF_COL, PAIR_SUPP_WRF_COL]].describe())
print("Any NaNs in pair_supp_wRF_VM?", df_corr[PAIR_SUPP_WRF_COL].isna().mean())

# plot and save correlations metrics
for filler in range(1):
    # Extract columns
    # x = df_corr['pair_suppression'].values
    x = df_corr['pair_supp_wRF_VM'].values
    y = df_corr['pair_selectivity'].values
    c = df_corr['corr_paired'].values  # correlation values
    c_singleLoc0 = df_corr['corr_single_loc0'].values
    c_singleLoc1 = df_corr['corr_single_loc1'].values

    savemat("corr_data.mat", {
        "nonPrefSupp": df_corr['pair_supp_wRF_VM'].values,
        "selectivity": df_corr['pair_selectivity'].values,
        "pairedCorr": df_corr['corr_paired'].values,
        "singleCorrLoc0": df_corr['corr_single_loc0'].values,
        "singleCorrLoc1": df_corr['corr_single_loc1'].values,
        "monkey": df_corr['monkey'].values.astype('U')  # save as strings
    })

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


########################################################################################################################
######################################## Compute weight suppression from RF location ###################################
################################################## mapping for each pair ###############################################
########################################################################################################################
# ---- CONFIG: correct column names for corr_records ----
UNIT1_COL = 'unit_i'                        # neuron 1 ID
UNIT2_COL = 'unit_j'                        # neuron 2 ID
PREF1_COL = 'neuron 1 preferred stim'       # >0 => prefers loc 0 ; <0 => prefers loc 1
PREF2_COL = 'neuron 2 preferred stim'       # >0 => prefers loc 0 ; <0 => prefers loc 1

# Non-preferred suppression columns (alpha) from corr_records
# If your actual column names use a space ("non pref"), change these two lines.
N1_NONPREF_COL = 'neuron 1 nonpref suppression'
N2_NONPREF_COL = 'neuron 2 nonpref suppression'

# Previously added outputs
OUT_COL_N1 = 'neuron1_weight_supp'
OUT_COL_N2 = 'neuron2_weight_supp'
OUT_COL_PAIR = 'pair_geom_mean_supp'

# New outputs you requested
OUT_COL_N1_AW = 'neuron1_alpha_weight_supp'
OUT_COL_N2_AW = 'neuron2_alpha_weight_supp'
OUT_COL_PAIR_AW = 'pair_geom_mean_alpha_weight'

# --- NEW output column names (weak-weight fractions) ---
OUT_COL_N1_WEAK = 'neuron1_weak_weight_frac'   # min(w0,w1)/(w0+w1)
OUT_COL_N2_WEAK = 'neuron2_weak_weight_frac'   # min(w0,w1)/(w0+w1)
OUT_COL_PAIR_WEAK = 'pair_geom_mean_weak'      # geometric mean of the two columns


# ---- helpers ----
def to_dataframe(maybe_array_or_df):
    """Return a pandas DataFrame given a pandas DataFrame or a NumPy structured/record array."""
    if isinstance(maybe_array_or_df, pd.DataFrame):
        return maybe_array_or_df.copy()
    return pd.DataFrame(maybe_array_or_df)


def build_weights_map(results_array):
    """
    Build dict: unit_key -> np.array([w0, w1]) from results_array rows.
    Expected row layout (two common variants):
      With monkeyName (len=8):
        (unit_key, seshDate, monkeyName, unitID_raw, fit_ok, params, weights, peak_val)
      Without monkeyName (len=7):
        (unit_key, seshDate, unitID_raw, fit_ok, params, weights, peak_val)
    """
    wm = {}
    for row in results_array:
        if isinstance(row, np.ndarray) and row.ndim == 0:
            row = row.item()
        try:
            if len(row) >= 8:
                unit_key, _, _, _, fit_ok, _, weights, _ = row
            elif len(row) == 7:
                unit_key, _, _, fit_ok, _, weights, _ = row
            else:
                unit_key = str(row[0])
                weights = row[-2]
                fit_ok = next((x for x in row if isinstance(x, (bool, np.bool_))), False)
            if fit_ok and isinstance(weights, np.ndarray) and weights.size == 2:
                wm[str(unit_key)] = weights.astype(float)
        except Exception:
            continue
    return wm


def neuron_weight_supp(weights_map, unit_key, pref_value):
    """
    Weight suppression for one neuron:
      if prefers loc 0 -> w1 / (w0 + w1)
      if prefers loc 1 -> w0 / (w0 + w1)
      else -> np.nan
    """
    w = weights_map.get(str(unit_key))
    if w is None or not np.all(np.isfinite(w)) or w.size != 2:
        return np.nan
    w0, w1 = float(w[0]), float(w[1])
    denom = w0 + w1
    if not np.isfinite(denom) or denom <= 0:
        return np.nan
    if pref_value > 0:
        return w1 / denom
    elif pref_value < 0:
        return w0 / denom
    else:
        return np.nan


def neuron_alpha_weight_supp(weights_map, unit_key, pref_value, nonpref_alpha):
    """
    New suppression (alpha x weight ratio):
      if prefers loc 0 -> alpha * (w1 / w0)
      if prefers loc 1 -> alpha * (w0 / w1)
      else -> np.nan
    """
    if not np.isfinite(nonpref_alpha):
        return np.nan
    w = weights_map.get(str(unit_key))
    if w is None or not np.all(np.isfinite(w)) or w.size != 2:
        return np.nan
    w0, w1 = float(w[0]), float(w[1])

    if pref_value > 0:      # prefers loc 0
        ratio_denom = w0
        ratio_num   = w1
    elif pref_value < 0:    # prefers loc 1
        ratio_denom = w1
        ratio_num   = w0
    else:                   # ambiguous preference
        return np.nan

    if not np.isfinite(ratio_denom) or ratio_denom == 0:
        return np.nan
    ratio = ratio_num / ratio_denom
    return nonpref_alpha * ratio


def geom_mean_nonneg(a, b):
    """Geometric mean for nonnegative finite a, b; else NaN."""
    if np.isfinite(a) and np.isfinite(b) and a >= 0 and b >= 0:
        return np.sqrt(a * b)
    return np.nan


# 1) Build the weights lookup from your per-neuron results_array
weights_map = build_weights_map(results_array)

# 2) Convert corr_records to a DataFrame (if it isn't one already)
df_corr = to_dataframe(corr_records)

# 3) (Already in your code) neuron-level weight suppression and pair geom mean
df_corr[OUT_COL_N1] = [
    neuron_weight_supp(weights_map, u1, p1)
    for u1, p1 in zip(df_corr[UNIT1_COL], df_corr[PREF1_COL])
]
df_corr[OUT_COL_N2] = [
    neuron_weight_supp(weights_map, u2, p2)
    for u2, p2 in zip(df_corr[UNIT2_COL], df_corr[PREF2_COL])
]
df_corr[OUT_COL_PAIR] = [
    geom_mean_nonneg(a, b)
    for a, b in zip(df_corr[OUT_COL_N1], df_corr[OUT_COL_N2])
]

# 4) NEW: alpha-weighted suppression per neuron
# If your corr_records columns are named with a space "non pref", change N1_NONPREF_COL / N2_NONPREF_COL above.
df_corr[OUT_COL_N1_AW] = [
    neuron_alpha_weight_supp(weights_map, u1, p1, a1)
    for u1, p1, a1 in zip(df_corr[UNIT1_COL], df_corr[PREF1_COL], df_corr[N1_NONPREF_COL])
]
df_corr[OUT_COL_N2_AW] = [
    neuron_alpha_weight_supp(weights_map, u2, p2, a2)
    for u2, p2, a2 in zip(df_corr[UNIT2_COL], df_corr[PREF2_COL], df_corr[N2_NONPREF_COL])
]

# 5) NEW: geometric mean of the alpha-weighted suppressions
df_corr[OUT_COL_PAIR_AW] = [
    geom_mean_nonneg(a, b)
    for a, b in zip(df_corr[OUT_COL_N1_AW], df_corr[OUT_COL_N2_AW])
]

# plot scatter
for filler in range(1):
    # Extract columns
    # x = df_corr['pair_geom_mean_supp'].values
    # x = df_corr['pair_suppression'].values
    x = df_corr['pair_geom_mean_alpha_weight'].values
    y = df_corr['pair_selectivity'].values
    c = df_corr['corr_paired'].values  # correlation values

    # Save to .mat file with desired variable names
    savemat("corr_data.mat", {
        "nonPrefSupp": x,
        "selectivity": y,
        "pairedCorr": c
    })

    plt.figure(figsize=(7, 6))
    sc = plt.scatter(
        x, y,
        c=c,
        cmap='bwr',          # red = positive, blue = negative, white ~ 0
        vmin=-1, vmax=1,     # fix scale from -1 to 1
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



########################################################################################################################
######################################## Look at tuning similarity based off of ########################################
############################################  responses at both locations ##############################################
########################################################################################################################


def vm_eval_on_grid(params, deg_grid):
    """
    params: dict with keys A, mu_deg, kappa, B
    deg_grid: array of angles in degrees
    returns f(deg_grid)
    """
    A = params["A"]; mu_deg = params["mu_deg"]; kappa = params["kappa"]; B = params["B"]
    theta = np.deg2rad(deg_grid)
    mu = np.deg2rad(mu_deg % 360.0)
    return A * np.exp(kappa * np.cos(theta - mu)) + B


def circ_diff_deg(a_deg, b_deg):
    """Smallest absolute circular difference in degrees, in [0,180]."""
    d = (a_deg - b_deg + 180.0) % 360.0 - 180.0
    return abs(d)


def overlap_coefficient(p, q, deg_grid):
    """
    p, q nonnegative arrays over degrees. We first normalize to integrate to 1.
    Returns ∫ min(p,q) dθ (Riemann sum in degrees), in [0,1].
    """
    p = np.clip(p, 0, np.inf)
    q = np.clip(q, 0, np.inf)
    # normalize to area 1; if flat/zero, return 0 overlap
    dp = np.trapz(p, deg_grid)
    dq = np.trapz(q, deg_grid)
    if dp <= 0 or dq <= 0:
        return 0.0
    p = p / dp
    q = q / dq
    return float(np.trapz(np.minimum(p, q), deg_grid))


def safe_pearson(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if np.all(~np.isfinite(x)) or np.all(~np.isfinite(y)):
        return np.nan
    # if no variance, define as nan
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan
    r, _ = stats.pearsonr(x, y)
    return float(r)


# --- Build fast lookup from your allUnitsLocationTuning list ---
# Expected entries like:
# {"unit_id": "230607_77", "monkey": "Meetz", "loc0_params": {...}, "loc1_params": {...}}
lut = {d["unit_id"]: d for d in allUnitsLocationTuning}

# --- Main pairwise similarity computation ---
deg_grid = np.arange(0.0, 360.0, 1.0)  # 1-degree grid; make finer if you like
similarity_records = []

for date in fileList:
    monkeyName, sessionDate = date.split('_')

    # Good units for this day (you already have this)
    goodUnits_today = [u for u in masterGoodUnits if u.startswith(sessionDate + "_")]

    # Unique pairs of 2 units
    for u1, u2 in combinations(goodUnits_today, 2):
        rec = {
            "pair": (u1, u2),
            "monkey": monkeyName,
            "sessionDate": sessionDate,
            "loc0": {},
            "loc1": {}
        }
        missing = False
        for uid in (u1, u2):
            if uid not in lut:
                missing = True
                break
        if missing:
            # Skip or record as incomplete
            continue

        # For each location, evaluate curves and compute metrics
        for loc_key in ("loc0_params", "loc1_params"):
            # Fetch params
            p1 = lut[u1][loc_key]
            p2 = lut[u2][loc_key]

            # If a fit failed earlier (NaNs), handle gracefully
            if any([not np.isfinite(p1[k]) for k in ("A","mu_deg","kappa","B")]) or \
               any([not np.isfinite(p2[k]) for k in ("A","mu_deg","kappa","B")]):
                loc_res = {"dmu_deg": np.nan, "dkappa": np.nan, "overlap": np.nan, "corr": np.nan}
            else:
                f1 = vm_eval_on_grid(p1, deg_grid)
                f2 = vm_eval_on_grid(p2, deg_grid)

                # baseline-remove for shape similarity
                f1_nr = np.clip(f1 - p1["B"], 0, np.inf)
                f2_nr = np.clip(f2 - p2["B"], 0, np.inf)

                dmu = circ_diff_deg(p1["mu_deg"], p2["mu_deg"])
                dk = abs(p1["kappa"] - p2["kappa"])
                ovl = overlap_coefficient(f1_nr, f2_nr, deg_grid)
                rr = safe_pearson(f1_nr, f2_nr)

                loc_res = {"dmu_deg": dmu, "dkappa": dk, "overlap": ovl, "corr": rr}

            if loc_key == "loc0_params":
                rec["loc0"] = loc_res
            else:
                rec["loc1"] = loc_res

        # Aggregate a few summaries across locations
        loc0, loc1 = rec["loc0"], rec["loc1"]
        rec["summary"] = {
            "overlap_mean": np.nanmean([loc0["overlap"], loc1["overlap"]]),
            "overlap_min": np.nanmin([loc0["overlap"], loc1["overlap"]]),
            "corr_mean": np.nanmean([loc0["corr"], loc1["corr"]]),
            "corr_min": np.nanmin([loc0["corr"], loc1["corr"]]),
            "dmu_mean_deg": np.nanmean([loc0["dmu_deg"], loc1["dmu_deg"]]),
            "dkappa_mean": np.nanmean([loc0["dkappa"], loc1["dkappa"]]),
        }

        similarity_records.append(rec)

### plot distribution of delta alignment of tuning curves for each monkey
# Extract Δμ values
for filler in range(1):
    dmu_meetz = [rec["summary"]["dmu_mean_deg"] for rec in similarity_records
                 if rec["monkey"].lower() == "meetz"]
    dmu_akshan = [rec["summary"]["dmu_mean_deg"] for rec in similarity_records
                  if rec["monkey"].lower() == "akshan"]

    plt.figure(figsize=(7,5))

    # Histogram, bins in 10° steps from 0–180
    bins = np.arange(0, 185, 10)

    plt.hist(dmu_meetz, bins=bins, alpha=0.6, label="Meetz", color="royalblue")
    plt.hist(dmu_akshan, bins=bins, alpha=0.6, label="Akshan", color="darkorange")

    plt.xlabel("Δμ (mean preferred direction difference, deg)")
    plt.ylabel("Number of unit pairs")
    plt.title("Distribution of alignment differences (Δμ) by monkey")
    plt.legend()
    plt.tight_layout()
    plt.show()

########################################################################################################################
######################################## Compare correlations with differences in ######################################
################################################## direction tuning ####################################################
########################################################################################################################
# 1) Build a lookup from similarity_records
sim_lut = {}
for rec in similarity_records:
    pair = rec["pair"]  # e.g., ("230607_77", "230607_81")
    key = frozenset(pair)
    sim_lut[key] = {
        "monkey": rec["monkey"],
        "dmu_mean_deg": rec["summary"]["dmu_mean_deg"],
        "overlap_mean": rec["summary"]["overlap_mean"],
    }

# 2) Choose which correlation field to use
# Options seen in your example: 'corr_paired', 'corr_single_loc0', 'corr_single_loc1'
corr_field = "corr_paired"   # <- change to "corr_single_loc0" or "corr_single_loc1" if desired

# 3) Collect per-record points (every corr_records entry)
X_dmu_meetz, Y_corr_meetz = [], []
X_dmu_akshan, Y_corr_akshan = [], []
X_ovl_meetz,  Y2_corr_meetz = [], []
X_ovl_akshan, Y2_corr_akshan = [], []

missing_pairs = 0

for rec in corr_records:
    # Build unordered key for the pair
    if "pair_key" in rec and isinstance(rec["pair_key"], (tuple, list)):
        key = frozenset(rec["pair_key"])
    else:
        key = frozenset((rec["unit_i"], rec["unit_j"]))

    # Lookup similarity metrics
    sim = sim_lut.get(key, None)
    if sim is None:
        missing_pairs += 1
        continue

    dmu = sim["dmu_mean_deg"]
    ovl = sim["overlap_mean"]
    corr_val = rec.get(corr_field, np.nan)

    if not np.isfinite(corr_val) or not np.isfinite(dmu) or not np.isfinite(ovl):
        continue

    if sim["monkey"].lower() == "meetz":
        X_dmu_meetz.append(dmu)
        Y_corr_meetz.append(corr_val)
        X_ovl_meetz.append(ovl)
        Y2_corr_meetz.append(corr_val)
    else:
        X_dmu_akshan.append(dmu);
        Y_corr_akshan.append(corr_val)
        X_ovl_akshan.append(ovl)
        Y2_corr_akshan.append(corr_val)

print(f"Skipped {missing_pairs} corr_records with no similarity match")


########################################################################################################################
############################ BASED OFF OF DIRECTION TUNING FOR DIRECTION TUNING MAPPING (GRF2) #########################
########################################################################################################################

########################################################################################################################
######################################## Look at tuning similarity based off of ########################################
################################################## direction tuning ####################################################
########################################################################################################################

# Von Mises model (with baseline and amplitude)
def von_mises_func(theta, baseline, amplitude, pref_dir, kappa):
    # theta, pref_dir in radians
    # normalized von Mises bump + baseline
    return baseline + amplitude * np.exp(kappa * np.cos(theta - pref_dir)) / (2 * np.pi * i0(kappa))


# Fit one neuron helper ---
def fit_von_mises_single(responses, theta_rad):
    # Initial guesses
    baseline_guess = np.nanmin(responses)
    amp_guess = np.nanmax(responses) - baseline_guess
    pref_guess = theta_rad[np.nanargmax(responses)]
    kappa_guess = 2.0

    p0 = [baseline_guess, amp_guess, pref_guess, kappa_guess]

    # Bounds: baseline free, amplitude >= 0, pref_dir in [0, 2π], kappa >= 1e-6
    bounds = (
        [-np.inf, 0.0, 0.0, 1e-6],
        [np.inf, np.inf, 2*np.pi, 1000.0]
    )

    # Fit
    params, _ = curve_fit(
        von_mises_func, theta_rad, responses, p0=p0,
        bounds=bounds, maxfev=20000
    )
    baseline, amplitude, pref_dir, kappa = params

    # Normalize preferred direction to [0, 2π)
    pref_dir = np.mod(pref_dir, 2*np.pi)
    return baseline, amplitude, pref_dir, kappa


def ecdf(y):
    y = np.sort(np.asarray(y))
    if y.size == 0:
        return np.array([]), np.array([])
    x = np.arange(1, y.size + 1) / y.size
    return y, x


# go through each session with multiple units and extract the direction tuning from the cells of interest
similarity_records = []  # will accumulate dicts across all sessions
for date in fileList:
    monkeyName, sessionDate = date.split('_')
    os.chdir(f'{monkeyName}/{sessionDate}/Direction Tuning')
    dirData = np.load('unitsDirTuningMat.npy')  # shape: (n_units_today, n_dirs)
    sessionUnitIDs = np.load('sessionUnitIDs.npy', allow_pickle=True)
    sessionUnitIDs = np.asarray(sessionUnitIDs).astype(str).tolist()  # ensure list[str]

    # Sanity: rows of dirData should match sessionUnitIDs length
    if dirData.shape[0] != len(sessionUnitIDs):
        print(f"[WARN] {sessionDate}: dirData rows ({dirData.shape[0]}) "
              f"!= sessionUnitIDs length ({len(sessionUnitIDs)})")

    # Good units for this day (strings like '230607_77')
    goodUnits_today = [u for u in masterGoodUnits if u.startswith(sessionDate + "_")]

    # Build mapping using sessionUnitIDs (NOT masterAllUnits)
    unit_to_row = {u: i for i, u in enumerate(sessionUnitIDs)}

    # Keep only good units that exist in this session’s ID list
    goodUnits_today = [u for u in goodUnits_today if u in unit_to_row]

    # If fewer than 2 good units, nothing to compare
    if len(goodUnits_today) < 2:
        os.chdir('../../../')
        continue

    # --- Angle grid for tuning fits ---
    n_rows, n_dirs = dirData.shape
    angles_deg = np.linspace(0, 360, n_dirs, endpoint=False)
    angles_rad = np.deg2rad(angles_deg)

    # --- Fit von Mises for good units using sessionUnitIDs indexing ---
    pref_dirs_rad = np.full(len(goodUnits_today), np.nan)  # start as NaN; fill on success

    for gi, u in enumerate(goodUnits_today):
        idx = unit_to_row[u]  # row in dirData per sessionUnitIDs
        if idx < 0 or idx >= n_rows:
            print(f"[SKIP] {sessionDate}: unit {u} row idx {idx} out of bounds")
            continue

        yi = dirData[idx, :]
        try:
            _, _, prefD, _ = fit_von_mises_single(yi, angles_rad)
            pref_dirs_rad[gi] = prefD
        except RuntimeError as e:
            print(f"[FIT FAIL] {sessionDate}: von Mises fit failed for {u} (row {idx}): {e}")

    # Only proceed with units that got a valid PD
    valid_idx = [k for k in range(len(goodUnits_today)) if np.isfinite(pref_dirs_rad[k])]
    if len(valid_idx) < 2:
        os.chdir('../../../')
        continue

    # Compute similarity among valid good units only
    th = pref_dirs_rad[valid_idx]
    c = np.cos(th).reshape(-1, 1)
    s = np.sin(th).reshape(-1, 1)
    SI = (c @ c.T) + (s @ s.T)          # [-1, 1]
    SI01 = 0.5 * (SI + 1.0)             # [0, 1]

    pref_dirs_deg = np.rad2deg(pref_dirs_rad)

    # Record only unique pairs among valid indices
    for a, b in combinations(range(len(valid_idx)), 2):
        i = valid_idx[a]
        j = valid_idx[b]
        similarity_records.append({
            "date": sessionDate,
            "monkey": monkeyName,
            "unit_i": goodUnits_today[i],
            "unit_j": goodUnits_today[j],
            "PD_i_deg": float(pref_dirs_deg[i]),
            "PD_j_deg": float(pref_dirs_deg[j]),
            "SI": float(SI[a, b]),
            "SI01": float(SI01[a, b]),
        })

    # go back to root
    os.chdir('../../../')

# plot tuning similarity distribution
# ---------- Build DataFrame ----------
sim_df = pd.DataFrame(similarity_records)
if "SI" not in sim_df or sim_df.empty:
    raise ValueError("similarity_records is empty or missing 'SI'.")
for filler in range(1):
    sim_df = sim_df[np.isfinite(sim_df["SI"])].copy()
    if sim_df.empty:
        raise ValueError("All 'SI' values are NaN/inf. Check how similarities are attached.")

    colors = {"Meetz": "tab:blue", "Akshan": "tab:orange"}

    # ---------- Plot ----------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # (A) Overall histogram (counts)
    bins_all = np.linspace(-1, 1, 21)
    axes[0].hist(sim_df["SI"], bins=bins_all, edgecolor="k", alpha=0.85)
    axes[0].set_xlabel("Similarity (cos Δθ, -1 to 1)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Overall distribution (all sessions)")
    axes[0].axvline(0, color="k", linestyle="--", alpha=0.5)

    # (B) Per-monkey density histograms
    bins_monkey = np.linspace(-1, 1, 21)
    for mk, g in sim_df.groupby("monkey"):
        vals = g["SI"].dropna().to_numpy()
        if vals.size == 0:
            continue
        axes[1].hist(
            vals, bins=bins_monkey, density=True, histtype="stepfilled",
            alpha=0.40, color=colors.get(mk, "gray"), edgecolor="k", linewidth=1.1,
            label=f"{mk} (n={vals.size})"
        )
    axes[1].set_xlabel("Similarity (cos Δθ, -1 to 1)")
    axes[1].set_ylabel("Density (area = 1)")
    axes[1].set_title("Distribution by monkey (density)")
    axes[1].axvline(0, color="k", linestyle="--", alpha=0.5)
    axes[1].legend(fontsize=9)

    # (C) Per-session ECDFs, colored by monkey (two colors only)
    # Each session/date is a curve; color chosen by that session's monkey
    for date, g in sim_df.groupby("date"):
        vals = g["SI"].dropna().to_numpy()
        if vals.size < 2:
            continue
        # Identify monkey for this session (assumes one monkey per date)
        mk = g["monkey"].iloc[0]
        y, x = ecdf(vals)
        axes[2].plot(y, x, alpha=0.55, linewidth=1.5, color=colors.get(mk, "gray"))

    axes[2].set_xlabel("Similarity (cos Δθ)")
    axes[2].set_ylabel("ECDF")
    axes[2].set_title("Per-session ECDF (colored by monkey)")
    axes[2].axvline(0, color="k", linestyle="--", alpha=0.4)

    # Two-color legend (no per-session clutter)
    handles = [
        Line2D([0], [0], color=colors["Meetz"], lw=2, label="Meetz sessions"),
        Line2D([0], [0], color=colors["Akshan"], lw=2, label="Akshan sessions"),
    ]
    axes[2].legend(handles=handles, loc="lower right", fontsize=9)

    plt.tight_layout()
    plt.show()

########################################################################################################################
######################################## Correlations and Similarity in Tuning #########################################
########################################################################################################################
# Make a DataFrame from similarity_records
sim_df = pd.DataFrame(similarity_records)
# Create a key for each pair (order-independent)
sim_df["pair_key"] = sim_df.apply(
    lambda r: tuple(sorted((r["unit_i"], r["unit_j"]))), axis=1
)
# Reduce to one row per pair (SI, SI01 are same for all)
sim_df = sim_df.drop_duplicates("pair_key")[["pair_key", "SI", "SI01"]]

# Turn corr_records into DataFrame
corr_df = pd.DataFrame(corr_records)
# Merge on pair_key
merged_df = corr_df.merge(sim_df, on="pair_key", how="left")

# plots
# Bin by similarity
bins = np.linspace(-1, 1, 9)  # 8 bins
pair_means["bin"] = pd.cut(pair_means["SI"], bins)

# Compute both mean and median per bin
bin_stats = pair_means.groupby("bin")["mean_corr"].agg(["mean", "median"])
bin_centers = [interval.mid for interval in bin_stats.index]

plt.figure(figsize=(6, 5))

# scatter all pairs
plt.scatter(pair_means["SI"], pair_means["mean_corr"], alpha=0.3)

# plot binned mean
plt.plot(bin_centers, bin_stats["mean"].values, "o-", color="red", label="Mean")

# plot binned median
plt.plot(bin_centers, bin_stats["median"].values, "s--", color="blue", label="Median")

plt.xlabel("Tuning similarity (cos Δθ, -1 to 1)")
plt.ylabel("Correlation")
plt.title("Binned correlation vs similarity")
plt.legend()
plt.show()

########## MEETZ ONLY
# keep only Meetz
merged_df = corr_df.merge(sim_df, on="pair_key", how="left")
merged_dfMeetz = merged_df[merged_df["monkey"] == "Meetz"]

# Now compute pair_means on Meetz-only data
pair_means = merged_dfMeetz.groupby("pair_key").agg(
    mean_corr=("corr_paired", "mean"),
    SI=("SI", "first"),
    SI01=("SI01", "first")
).reset_index()

# Bin by similarity
bins = np.linspace(-1, 1, 9)
pair_means["bin"] = pd.cut(pair_means["SI"], bins)

# Compute mean & median per bin
bin_stats = pair_means.groupby("bin")["mean_corr"].agg(["mean", "median"])
bin_centers = [interval.mid for interval in bin_stats.index]

# Plot
plt.figure(figsize=(6,5))
plt.scatter(pair_means["SI"], pair_means["mean_corr"], alpha=0.3, label="Pairs")
plt.plot(bin_centers, bin_stats["mean"].values, "o-", color="red", label="Mean")
plt.plot(bin_centers, bin_stats["median"].values, "s--", color="blue", label="Median")

plt.xlabel("Tuning similarity (cos Δθ, -1 to 1)")
plt.ylabel("Correlation")
plt.title("Binned correlation vs similarity (Meetz only)")
plt.legend()
plt.show()

########## Akshan ONLY
merged_df = corr_df.merge(sim_df, on="pair_key", how="left")
merged_dfAkshan = merged_df[merged_df["monkey"] == "Akshan"]

# Now compute pair_means on Meetz-only data
pair_means = merged_dfAkshan.groupby("pair_key").agg(
    mean_corr=("corr_paired", "mean"),
    SI=("SI", "first"),
    SI01=("SI01", "first")
).reset_index()

# Bin by similarity
bins = np.linspace(-1, 1, 9)
pair_means["bin"] = pd.cut(pair_means["SI"], bins)

# Compute mean & median per bin
bin_stats = pair_means.groupby("bin")["mean_corr"].agg(["mean", "median"])
bin_centers = [interval.mid for interval in bin_stats.index]

# Plot
plt.figure(figsize=(6,5))
plt.scatter(pair_means["SI"], pair_means["mean_corr"], alpha=0.3, label="Pairs")
plt.plot(bin_centers, bin_stats["mean"].values, "o-", color="red", label="Mean")
plt.plot(bin_centers, bin_stats["median"].values, "s--", color="blue", label="Median")

plt.xlabel("Tuning similarity (cos Δθ, -1 to 1)")
plt.ylabel("Correlation")
plt.title("Binned correlation vs similarity (Akshan only)")
plt.legend()
plt.show()

