# Imports
from usefulFns import *
from normalizationFunctions import *

p0 = time.time()

# # Both monkeys
fileList = ['Meetz_230607', 'Meetz_230608', 'Meetz_230613', 'Meetz_230615',
            'Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628',
            'Meetz_230630', 'Meetz_230707', 'Meetz_230710', 'Meetz_230711',
            'Meetz_230718', 'Meetz_230719', 'Meetz_230720', 'Akshan_240530',
            'Akshan_240603', 'Akshan_240606', 'Akshan_240607', 'Akshan_240610',
            'Akshan_240611', 'Akshan_240628', 'Akshan_240703', 'Akshan_240704',
            'Akshan_240705', 'Akshan_240709', 'Akshan_241014', 'Akshan_241024']

unitList = ['230607_149', '230607_147', '230607_146', '230607_144', '230607_126',
            '230607_125', '230607_80', '230607_77', '230608_117', '230608_140',
            '230608_143', '230608_149', '230613_187', '230615_166', '230620_58',
            '230626_58', '230626_63', '230626_64', '230626_65', '230626_131',
            '230627_14', '230627_70', '230627_71', '230627_79', '230627_93',
            '230627_106', '230627_112', '230627_141', '230627_147', '230627_148',
            '230627_152', '230628_42', '230628_143', '230630_50', '230630_138',
            '230707_140', '230707_143', '230707_144', '230710_137', '230711_45',
            '230711_50', '230718_178', '230719_112', '230719_147', '230719_155',
            '230719_156', '230720_149', '240530_39', '240603_127', '240606_74',
            '240606_84', '240606_90', '240606_94', '240606_103', '240607_31',
            '240610_29', '240610_32', '240610_37', '240610_92', '240611_22',
            '240611_31', '240611_37', '240611_42', '240611_52', '240611_93',
            '240611_95', '240611_98', '240611_104', '240611_112', '240611_116',
            '240611_125', '240611_127', '240611_129', '240628_49', '240703_7',
            '240703_57', '240704_24', '240705_56', '240705_60', '240705_64',
            '240705_65', '240709_72', '240709_108', '241014_3', '241014_15',
            '241014_61', '241024_2', '241024_21', '241024_23', '241024_25']

# # just meetz
# fileList = ['Meetz_230607', 'Meetz_230608', 'Meetz_230613', 'Meetz_230615',
#             'Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628', 'Meetz_230630',
#             'Meetz_230707', 'Meetz_230710', 'Meetz_230711', 'Meetz_230718', 'Meetz_230719',
#             'Meetz_230720']
#
# unitList = ['230607_149', '230607_147', '230607_146', '230607_144', '230607_126',
#             '230607_125', '230607_80', '230607_77', '230608_117', '230608_140',
#             '230608_143', '230608_149', '230613_187', '230615_166', '230620_58',
#             '230626_58', '230626_63', '230626_64', '230626_65', '230626_131',
#             '230627_14', '230627_70', '230627_71', '230627_79', '230627_93',
#             '230627_106', '230627_112', '230627_141', '230627_147', '230627_148',
#             '230627_152', '230628_42', '230628_143', '230630_50', '230630_138',
#             '230707_140', '230707_143', '230707_144', '230710_137', '230711_45',
#             '230711_50', '230718_178', '230719_112', '230719_147', '230719_155',
#             '230719_156', '230720_149']

# # just akshan
# fileList = ['Akshan_240603', 'Akshan_240606', 'Akshan_240607', 'Akshan_240610',
#             'Akshan_240611', 'Akshan_240628', 'Akshan_240703', 'Akshan_240704',
#             'Akshan_240705', 'Akshan_240709', 'Akshan_241014', 'Akshan_241024']
#
# unitList = ['240530_39', '240603_127', '240606_74',
#             '240606_84', '240606_90', '240606_94', '240606_103', '240607_31',
#             '240610_29', '240610_32', '240610_37', '240610_92', '240611_22',
#             '240611_31', '240611_37', '240611_42', '240611_52', '240611_93',
#             '240611_95', '240611_98', '240611_104', '240611_112', '240611_116',
#             '240611_125', '240611_127', '240611_129', '240628_49', '240703_7',
#             '240703_57', '240704_24', '240705_56', '240705_60', '240705_64',
#             '240705_65', '240709_72', '240709_108', '241014_3', '241014_15',
#             '241014_61', '241024_2', '241024_21', '241024_23', '241024_25']

# fileList = ['Meetz_230719']
# unitList = ['230719_147']

masterPupilDelta = []
masterImDelta = []
trialPupilAreaList = []
trialPupilAreaTime = []
trialPAAlignedStim1 = []
trialPAAlignedStim1Time = []
trialPAAlignedStim1Neuron = []
trialPAAlignedStim1TimeNeuron = []
pupilAreaAcrossBlockPop = []
masterAboveMedianSplit = []
masterBelowMedianSplit = []
popBlockMeanZscore = []
popTrialStim1 = []
popTrialStim2 = []
popTrialStim3 = []
popTrialStim4 = []
popTrialStim1PSTH = []
popTrialStim2PSTH = []
popTrialStim3PSTH = []
popTrialStim4PSTH = []
globalTrialIM = []

for file in fileList:
    # Load relevant file here with pyMat reader
    monkeyName, seshDate = file.split('_')

    fileName = f'{monkeyName}_{seshDate}_MTND_Spikes.mat'
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

    # list of indices of correctTrials (non-instruct, valid trialCertify)
    corrTrials = correctTrialsMTX(allTrials)

    # generate list of unique active units, and their channel
    units = activeUnits('spikeData', allTrials)
    unitCluster = allTrials[corrTrials[0]]['spikeTempInfo']['cgs']
    unitsChannel = unitsInfo(units, corrTrials, allTrials)

    # create list of Good Units from this session (compiled from goodUnits list)
    sessionGoodUnits = []
    for unit in unitList:
        date, unitID = unit.split('_')
        if date == seshDate:
            sessionGoodUnits.append(int(unitID))
    sessionGoodUnits = np.array(sessionGoodUnits)

    # change stimDesc to be list of dictionaries
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        nStim = len(currTrial['stimDesc']['data']['listType'])
        currTrial['stimDesc']['data'] = [{k: v[i] for k, v in currTrial['stimDesc']['data'].items()}
                                         for i in range(nStim)]

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
    trueStimDurMS = 250
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
        trueStimDurMS = np.int32(np.around(1000 / frameRateHz * stimDurFrame[0]))

    # get mapping stimulus azi/ele positions and index
    stepIndex = []
    loc1AziEle = []
    offsetDegSep = []
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        stimDesc = currTrial['stimDesc']['data']
        for stim in stimDesc:
            if stim['stimLoc'] == 1 and stim['stepIndex'] not in stepIndex:
                stepIndex.append(stim['stepIndex'])
                loc1AziEle.append((stim['azimuthDeg'], stim['elevationDeg']))
                offsetDegSep.append(stim['offsetDeg'])

    stepIndex = np.array(stepIndex)
    loc1AziEle = np.array(loc1AziEle)
    offsetDegSep = np.array(offsetDegSep)
    offsetAziEle = loc1AziEle[np.argsort(stepIndex)]
    offsetDegSep = offsetDegSep[np.argsort(stepIndex)]
    eccentricity = np.sqrt(offsetAziEle[0][0]**2 + offsetAziEle[0][1]**2)
    rfGaborSigma = allTrials[0]['rfGabor']['data']['sigmaDeg'] / eccentricity
    offsetDegSepNorm = (offsetDegSep[np.argsort(offsetDegSep)]) / eccentricity

    # initialize lists/arrays/dataframes for counting spikeCounts and for analysis
    if 'prefDirDeg' in header['blockStatus']['data']:
        prefDir = header['blockStatus']['data']['prefDirDeg'][0]
        nonPrefDir = header['blockStatus']['data']['nonPrefDirDeg'][0]
    else:
        prefDir = header['mapSetting']['data']['directionDeg'][0]
        nonPrefDir = (prefDir + 180) % 360
    blocksDone = allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone'][0]
    numOffsets = header['blockStatus']['data']['numOffsets'][0]
    numSteps = int((numOffsets - 1) / 2)
    spikeCountMat = np.zeros((len(units), blocksDone+1, numSteps*7+3))
    pupilAreaMat = np.zeros((blocksDone+1, numSteps*7+3))
    trialStimPosition = np.zeros((blocksDone+1, numSteps*7+3))
    stimCount = np.zeros((3, numOffsets), dtype=int)
    stimCountIndex = np.arange((numSteps*7)+3-numSteps)
    stimCountIndex = stimCountIndex.reshape(3, int(len(stimCountIndex)/3))
    transectStimCount = np.zeros(numSteps, dtype=int)
    transectCountIndex = np.arange((numSteps*7)+3-numSteps, numSteps*7+3)
    onLatency = 50 / 1000  # time in MS for counting window latency after stim on
    offLatency = 50 / 1000  # time in MS for counting window latency after stim off
    histPrePostMS = 100  # 100ms window pre/post stimulus on/off
    spikeHists = np.zeros((len(units), numSteps*7+3, trueStimDurMS + (2*histPrePostMS+1)))

    # insert spikes from valid stimulus presentations into spike count matrices
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        if 'spikeData' in currTrial:
            stimDesc = currTrial['stimDesc']['data']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            trialStartS = currTrial['taskEvents']['trialStart']['timeS']
            pupilAreaTimeS = trialStartS + 2/1000 * np.arange(len(currTrial['eyeLPData']['data']))
            if currTrial['taskEvents']['catchTrial']['data'] != 1:
                targetOnTime = currTrial['taskEvents']['targetOn']['time']
                tempIndx = np.where(pupilAreaTimeS < targetOnTime)
                trialPupilAreaList.append(currTrial['eyeLPData']['data'][tempIndx])
                tempPupilTime = 2/1000 * np.arange(len(currTrial['eyeLPData']['data']))
                trialPupilAreaTime.append(tempPupilTime[tempIndx])
                tempIndx = np.where((pupilAreaTimeS > stim1TimeS) & (pupilAreaTimeS < targetOnTime))
                trialPAAlignedStim1.append(currTrial['eyeLPData']['data'][tempIndx])
                trialPAAlignedStim1Time.append(2/1000 * np.arange(len(currTrial['eyeLPData']['data'][tempIndx])))
            trialStimNumber = 0
            for stim in stimDesc:
                if stim['listType'] == 1 and stim['stimLoc'] == 0:
                    centIndex = stim['centerIndex']
                    offIdx = stim['offsetIndex']
                    step = stim['stepIndex']
                    stimOnTimeS = ((1000 / frameRateHz * stim['stimOnFrame'])
                                   / 1000) + stim1TimeS
                    stimOffTimeS = ((1000 / frameRateHz * stim['stimOffFrame'])
                                    / 1000) + stim1TimeS
                    if step > numSteps:
                        indx = step - numSteps - 1
                        stCount = int(transectStimCount[indx])
                        transectStimCount[indx] += 1
                        stimIndex = transectCountIndex[indx]
                    else:
                        col = (int(offIdx != 0) * step * 2 - int(offIdx == 1))
                        # ^ basically indexing into stimCount stepIndex * 2 - 1
                        # for every column after the first, which is the blank
                        # condition
                        stCount = int(stimCount[centIndex, col])
                        stimCount[centIndex, col] += 1
                        stimIndex = stimCountIndex[centIndex, col]

                    # pupil Area
                    pupilIndx = np.where((pupilAreaTimeS >= stimOnTimeS) &
                                         (pupilAreaTimeS <= stimOffTimeS))[0]
                    stimMeanPupilArea = np.mean(currTrial['eyeLPData']['data'][pupilIndx])
                    pupilAreaMat[stCount][stimIndex] = stimMeanPupilArea
                    trialStimPosition[stCount][stimIndex] = trialStimNumber

                    # only units that are considered good for population analysis
                    # for unit in sessionGoodUnits:
                    for unit in units:
                        unitCount = np.where(units == unit)[0][0]
                        if unit in currTrial['spikeData']['unit']:
                            unitIndex = np.where(currTrial['spikeData']['unit'] == unit)[0]
                            unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                            stimSpikes = np.where((unitTimeStamps >= (stimOnTimeS + onLatency)) &
                                                  (unitTimeStamps <= (stimOffTimeS + offLatency)))[0]
                            spikeCountMat[unitCount][stCount][stimIndex] \
                                = len(stimSpikes)

                            # PSTHs
                            stimOnPreSNEV = stimOnTimeS - (histPrePostMS / 1000)
                            stimOffPostSNEV = stimOffTimeS + (histPrePostMS / 1000)
                            histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                                             & (unitTimeStamps <= stimOffPostSNEV)
                                                             )] - stimOnPreSNEV
                            histStimSpikes = np.int32(histStimSpikes * 1000)
                            spikeHists[unitCount, stimIndex, histStimSpikes] += 1

                    trialStimNumber += 1

    # mean, SEM, and reshaping of spikeCount matrices
    meanSpike = np.mean(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSD = np.std(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSEM = spikeCountSD/np.sqrt(blocksDone)
    meanSpikeReshaped = np.zeros((len(units), 3, int(((numSteps*7)+3-numSteps)/3)))
    SEMReshaped = np.zeros((len(units), 3, int(((numSteps*7)+3-numSteps)/3)))
    for count, i in enumerate(meanSpikeReshaped):
        meanSpikeReshaped[count] = (meanSpike[count, :(numSteps*7)+3-numSteps].
                                    reshape(3, int(((numSteps*7)+3-numSteps)/3)) *
                                    1000/trueStimDurMS)
        SEMReshaped[count] = (spikeCountSEM[count, :(numSteps*7)+3-numSteps].
                              reshape(3, int(((numSteps*7)+3-numSteps)/3)) *
                              1000/trueStimDurMS)
    transectMeanSpike = meanSpike[:, (numSteps*7)+3-numSteps:] * 1000 / trueStimDurMS
    transectSEM = spikeCountSEM[:, (numSteps*7)+3-numSteps:] * 1000 / trueStimDurMS

    # add avg pupil area for each block to master list
    avgPupilAreaAcrossBlock = np.mean(pupilAreaMat[:blocksDone, :], axis=1)
    pupilAreaAcrossBlockPop.append(avgPupilAreaAcrossBlock)

    # known variables for fitting model that includes the baseline/spontaneous term
    contrast_centerSpon = [1, 0, 0, 0, 0,
                           1, 0, 0, 0, 0,
                           1, 1, 1, 1,
                           1, 1, 1, 1,
                           1, 1, 1, 1,
                           1, 1, 1, 1,
                           0]  # First response has contrast in the center, second in periphery
    contrast_peripherySpon = [0, 1, 1, 1, 1,
                              0, 1, 1, 1, 1,
                              1, 1, 1, 1,
                              1, 1, 1, 1,
                              1, 1, 1, 1,
                              1, 1, 1, 1,
                              0]  # First has no periphery, second has periphery with contrast
    locationsSpon = np.array([-1, 0, 1, 2, 3,
                              -1, 0, 1, 2, 3,
                              0, 1, 2, 3,
                              0, 1, 2, 3,
                              0, 1, 2, 3,
                              0, 1, 2, 3,
                              -1])  # First has no peripheral stimulus, second is at location 1
    stim_type_centerSpon = np.array([1, -1, -1, -1, -1,
                                     0, -1, -1, -1, -1,
                                     1, 1, 1, 1,
                                     0, 0, 0, 0,
                                     1, 1, 1, 1,
                                     0, 0, 0, 0,
                                     1])
    stim_type_peripherySpon = np.array([-1, 1, 1, 1, 1,
                                        -1, 0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        1, 1, 1, 1,
                                        1, 1, 1, 1,
                                        0, 0, 0, 0,
                                        1])

    # fit mean responses to model, then fit each block to model but only allow Im to vary
    # then see if changes in Im from mean Im varies with changes in pupil area from mean pupil area
    # 1. normalize pupil area mat

    pupilAreaMat_normalized = np.zeros_like(pupilAreaMat)
    n_trials, n_conds = pupilAreaMat.shape
    for trial in range(n_trials):
        for cond in range(n_conds):
            # find the position of this condition in this trial
            pos = trialStimPosition[trial, cond]

            # find all conditions in this trial that have the same position
            same_pos_mask = (trialStimPosition[trial, :] == pos)

            # mean pupil area at this position in this trial
            mean_pupil_at_pos = np.mean(pupilAreaMat[trial, same_pos_mask])

            # normalize this condition's pupil area
            pupilAreaMat_normalized[trial, cond] = pupilAreaMat[trial, cond] / mean_pupil_at_pos

    for unit in sessionGoodUnits:
        i = np.where(unit == units)[0][0]

        # ====== Indices for condition ordering ======
        pref_idx = [18, 2, 4, 6, 8]
        nonpref_idx = [9, 1, 3, 5, 7]
        prefNonpref_idx = [19, 21, 23, 25]
        nonprefPref_idx = [11, 13, 15, 17]
        pp_idx = [20, 22, 24, 26]
        nn_idx = [10, 12, 14, 16]
        spon_idx = [0]

        # ====== Full fit ======
        meanSpikes = meanSpikeReshaped[i].reshape(27)
        # normalize to pref at center
        meanSpikesNormalized = meanSpikes / meanSpikes[18]

        # properly ordered mean response
        pref = meanSpikesNormalized[pref_idx]
        nonpref = meanSpikesNormalized[nonpref_idx]
        prefNonpref = meanSpikesNormalized[prefNonpref_idx]
        nonprefPref = meanSpikesNormalized[nonprefPref_idx]
        pp = meanSpikesNormalized[pp_idx]
        nn = meanSpikesNormalized[nn_idx]
        spon = meanSpikesNormalized[spon_idx]

        resp = np.concatenate((
            pref, nonpref, prefNonpref, nonprefPref, pp, nn, spon
        ), axis=0)

        initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10, float(spon)]
        popt, pcov = curve_fit(
            model_wrapperHeuristic,
            (contrast_centerSpon, contrast_peripherySpon, locationsSpon, stim_type_centerSpon, stim_type_peripherySpon),
            resp,
            p0=initial_guess
        )

        L0_full = popt[6]

        # ====== Normalize all spike counts ======
        # normalized all spike counts, normalized to mean pref center activity
        normalizedSpikeCount = (spikeCountMat[i, :blocksDone, :27] * 1000/trueStimDurMS) / meanSpikes[18]
        # mean_pupil_area = pupilAreaMat[:blocksDone, :27].mean()
        mean_pupil_area = pupilAreaMat_normalized[:blocksDone, :27].mean()

        n_blocks = blocksDone
        trials_per_block = blocksDone // n_blocks

        L0_blocks = []
        pupil_blocks = []

        for b in range(n_blocks):
            start = b * trials_per_block
            end = (b + 1) * trials_per_block if b < n_blocks - 1 else blocksDone

            block_mean_full = normalizedSpikeCount[start:end, :].mean(axis=0)

            block_pref = block_mean_full[pref_idx]
            block_nonpref = block_mean_full[nonpref_idx]
            block_prefNonpref = block_mean_full[prefNonpref_idx]
            block_nonprefPref = block_mean_full[nonprefPref_idx]
            block_pp = block_mean_full[pp_idx]
            block_nn = block_mean_full[nn_idx]
            block_spon = block_mean_full[spon_idx]

            block_resp = np.concatenate((
                block_pref,
                block_nonpref,
                block_prefNonpref,
                block_nonprefPref,
                block_pp,
                block_nn,
                block_spon
            ), axis=0)

            # block_pupil = pupilAreaMat[start:end, :27].mean()
            block_pupil = pupilAreaMat_normalized[start:end, :27].mean()

            fixed_params = popt.copy()


            def model_L0only(L0):
                params = fixed_params.copy()
                params[6] = L0
                return model_wrapperHeuristic(
                    (contrast_centerSpon, contrast_peripherySpon,
                     locationsSpon, stim_type_centerSpon, stim_type_peripherySpon),
                    *params
                )


            res = minimize_scalar(
                lambda L0: np.sum((model_L0only(L0) - block_resp) ** 2),
                bounds=(0.001, 10), method='bounded'
            )

            L0_blocks.append(res.x)
            pupil_blocks.append(block_pupil)

        L0_blocks = np.array(L0_blocks)
        pupil_blocks = np.array(pupil_blocks)

        # ====== Compute deltas ======
        delta_L0_blocks = L0_blocks - L0_full
        delta_L0_blocks_z = (delta_L0_blocks - np.mean(delta_L0_blocks)) / np.std(delta_L0_blocks)  # z-score
        # delta_pupil_blocks = (pupil_blocks - mean_pupil_area) / np.std(pupilAreaMat[:blocksDone, :27])
        delta_pupil_blocks = (pupil_blocks - mean_pupil_area) / np.std(pupilAreaMat_normalized[:blocksDone, :27])


        # masterImDelta.extend(delta_L0_blocks)
        masterImDelta.extend(delta_L0_blocks_z)
        masterPupilDelta.extend(delta_pupil_blocks)

    # =============================================================================
    # # median split of pupil to see performance changes
    # 1. compute mean pupil area for every trial excluding ignored trials
    pupilAreaBeforeSaccade = []
    pupilAreaTrialIndx = []
    for count, currTrial in enumerate(allTrials):
        if currTrial['extendedEOT']['data'] != 4:
            trialStartMS = currTrial['trialStart']['timeMS']
            if currTrial['trial']['data']['catchTrial'] == 1:
                rewardTimeMS = currTrial['reward']['timeMS']
                if type(rewardTimeMS) != int:
                    rewardTimeMS = rewardTimeMS[0]
                timeDiff = (rewardTimeMS - trialStartMS) / 1000  # convert to seconds
            else:
                if 'saccade' not in currTrial:
                    continue
                else:
                    saccadeTimeMS = currTrial['saccade']['timeMS']
                    timeDiff = (saccadeTimeMS - trialStartMS) / 1000
            pupilAreaTimeMS = 2/1000 * np.arange(len(currTrial['eyeLPData']['data']))
            indx = np.where(pupilAreaTimeMS < timeDiff)[0]
            avgPupilAreaBeforeSaccade = np.mean(currTrial['eyeLPData']['data'])
            pupilAreaBeforeSaccade.append(avgPupilAreaBeforeSaccade)
            pupilAreaTrialIndx.append(count)

    # 2. perform median split
    pupilAreaBeforeSaccade = np.array(pupilAreaBeforeSaccade)
    pupilAreaTrialIndx = np.array(pupilAreaTrialIndx)
    median_val = np.median(pupilAreaBeforeSaccade)
    low_indx = np.where(pupilAreaBeforeSaccade <= median_val)[0]
    low_indices = pupilAreaTrialIndx[low_indx]
    high_indx = np.where(pupilAreaBeforeSaccade > median_val)[0]
    high_indices = pupilAreaTrialIndx[high_indx]

    # 3. compute performance for median split trials
    corrTrialsHighSplit = 0
    corrTrialsLowSplit = 0
    for count, currTrial in enumerate(allTrials):
        # skip if ignored trial (data == 4)
        if currTrial['extendedEOT']['data'] == 4:
            continue
        if 'saccade' not in currTrial:
            continue
        # correct trial
        if currTrial['extendedEOT']['data'] == 0:
            if count in high_indices:
                corrTrialsHighSplit += 1
            elif count in low_indices:
                corrTrialsLowSplit += 1

    print(f'pupil area above median split performance: {corrTrialsHighSplit/len(high_indices)}')
    print(f'pupil area below median split performance: {corrTrialsLowSplit/len(low_indices)}')
    masterAboveMedianSplit.append(corrTrialsHighSplit/len(high_indices))
    masterBelowMedianSplit.append(corrTrialsLowSplit/len(low_indices))

    # =============================================================================
    # z-score the spikeCountMat to see if mean rate changes across blocks
    for unit in sessionGoodUnits:
        i = np.where(unit == units)[0][0]
        X = spikeCountMat[i, :blocksDone, :]

        # z-score along rows (i.e., across trials) **for each column**
        X_zscored = stats.zscore(X, axis=0, ddof=1)  # ddof=1 for sample std

        # mean z-score for each block
        meanZscore = np.mean(X_zscored, axis=1)
        popBlockMeanZscore.append(meanZscore)

    # =============================================================================
    # generate responses for each stim to their avg response for avg trial
    allUnitsZscoreSpikeCount = stats.zscore(spikeCountMat[:, :blocksDone, :], axis=1, ddof=1)

    stimCount = np.zeros((3, numOffsets), dtype=int)
    stimCountIndex = np.arange((numSteps*7)+3-numSteps)
    stimCountIndex = stimCountIndex.reshape(3, int(len(stimCountIndex)/3))
    transectStimCount = np.zeros(numSteps, dtype=int)
    transectCountIndex = np.arange((numSteps*7)+3-numSteps, numSteps*7+3)
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        if currTrial['blockStatus']['data']['blocksDone'][0] == blocksDone:
            break
        if 'spikeData' in currTrial:
            stimDesc = currTrial['stimDesc']['data']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            trialStartS = currTrial['taskEvents']['trialStart']['timeS']
            trialStimNumber = 0
            for stim in stimDesc:
                if stim['listType'] == 1 and stim['stimLoc'] == 0:

                    centIndex = stim['centerIndex']
                    offIdx = stim['offsetIndex']
                    step = stim['stepIndex']
                    stimOnTimeS = ((1000 / frameRateHz * stim['stimOnFrame'])
                                   / 1000) + stim1TimeS
                    stimOffTimeS = ((1000 / frameRateHz * stim['stimOffFrame'])
                                    / 1000) + stim1TimeS
                    if step > numSteps:
                        indx = step - numSteps - 1
                        stCount = int(transectStimCount[indx])
                        transectStimCount[indx] += 1
                        stimIndex = transectCountIndex[indx]
                    else:
                        col = (int(offIdx != 0) * step * 2 - int(offIdx == 1))
                        # ^ basically indexing into stimCount stepIndex * 2 - 1
                        # for every column after the first, which is the blank
                        # condition
                        stCount = int(stimCount[centIndex, col])
                        stimCount[centIndex, col] += 1
                        stimIndex = stimCountIndex[centIndex, col]

                    # only units that are considered good for population analysis
                    for unit in sessionGoodUnits:
                        unitCount = np.where(units == unit)[0][0]
                        if unit in currTrial['spikeData']['unit']:
                            zScoreResp = allUnitsZscoreSpikeCount[unitCount, stCount, stimIndex]
                            if trialStimNumber == 0:
                                popTrialStim1.append(zScoreResp)
                            if trialStimNumber == 1:
                                popTrialStim2.append(zScoreResp)
                            if trialStimNumber == 2:
                                popTrialStim3.append(zScoreResp)
                            if trialStimNumber == 3:
                                popTrialStim4.append(zScoreResp)

                    trialStimNumber += 1

    # =============================================================================
    # find IM for every stim in a trial across trials
    # 1. fit model to data and store fit params to array
    paramFits = np.empty(len(units), dtype=object)
    for unit in sessionGoodUnits:
        i = np.where(unit == units)[0][0]
        meanSpikes = meanSpikeReshaped[i].reshape(27)
        # normalize to pref at center
        meanSpikesNormalized = meanSpikes / meanSpikes[18]

        pref = meanSpikesNormalized[[18, 2, 4, 6, 8]]
        nonpref = meanSpikesNormalized[[9, 1, 3, 5, 7]]
        prefNonpref = meanSpikesNormalized[[19, 21, 23, 25]]
        nonprefPref = meanSpikesNormalized[[11, 13, 15, 17]]
        pp = meanSpikesNormalized[[20, 22, 24, 26]]
        nn = meanSpikesNormalized[[10, 12, 14, 16]]
        spon = meanSpikesNormalized[0]

        # fit heuristic model to average data from each condition to get parameters
        resp = np.concatenate((pref, nonpref, prefNonpref, nonprefPref, pp, nn, [spon]), axis=0)
        initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                         float(spon)]
        popt, pcov = curve_fit(model_wrapperHeuristic, (contrast_centerSpon,
                                                        contrast_peripherySpon,
                                                        locationsSpon,
                                                        stim_type_centerSpon,
                                                        stim_type_peripherySpon), resp,
                               p0=initial_guess)
        paramFits[i] = popt

    stimCount = np.zeros((3, numOffsets), dtype=int)
    stimCountIndex = np.arange((numSteps*7)+3-numSteps)
    stimCountIndex = stimCountIndex.reshape(3, int(len(stimCountIndex)/3))
    transectStimCount = np.zeros(numSteps, dtype=int)
    transectCountIndex = np.arange((numSteps*7)+3-numSteps, numSteps*7+3)
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        if 'spikeData' in currTrial:
            stimDesc = currTrial['stimDesc']['data']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            trialStartS = currTrial['taskEvents']['trialStart']['timeS']
            pupilAreaTimeS = trialStartS + 2/1000 * np.arange(len(currTrial['eyeLPData']['data']))
            catchTrial = currTrial['taskEvents']['catchTrial']['data']

            trialIM_per_neuron = [[] for _ in sessionGoodUnits]
            for stim in stimDesc:
                if stim['listType'] == 3 and stim['stimLoc'] == 0:
                    centIndex = stim['centerIndex']
                    offIdx = stim['offsetIndex']
                    step = stim['stepIndex']
                    stimOnTimeS = ((1000 / frameRateHz * stim['stimOnFrame'])
                                   / 1000) + stim1TimeS
                    stimOffTimeS = ((1000 / frameRateHz * stim['stimOffFrame'])
                                    / 1000) + stim1TimeS
                    if step > numSteps:
                        indx = step - numSteps - 1
                        stimIndex = transectCountIndex[indx]
                    else:
                        col = (int(offIdx != 0) * step * 2 - int(offIdx == 1))
                        # ^ basically indexing into stimCount stepIndex * 2 - 1
                        # for every column after the first, which is the blank
                        # condition
                        stimIndex = stimCountIndex[centIndex, col]

                    if stimIndex <= 26 and catchTrial != 1:
                        # only units that are considered good for population analysis
                        for idx, unit in enumerate(sessionGoodUnits):
                            unitCount = np.where(units == unit)[0][0]
                            meanSpikes = meanSpikeReshaped[unitCount].reshape(27)
                            # normalized to pref center spike count

                            if unit in currTrial['spikeData']['unit']:
                                unitIndex = np.where(currTrial['spikeData']['unit'] == unit)[0]
                                unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                                stimSpikes = np.where((unitTimeStamps >= (stimOnTimeS + onLatency)) &
                                                      (unitTimeStamps <= (stimOffTimeS + offLatency)))[0]
                                trialSpikes = len(stimSpikes)
                            else:
                                trialSpikes = 0

                            trialSpikes = (trialSpikes * 1000/trueStimDurMS) / meanSpikes[18]
                            # find Im
                            popt = paramFits[unitCount]
                            im = compute_im(stimIndex, trialSpikes, popt)
                            if abs(im) > 10:
                                trialIM_per_neuron[idx].append([])
                            else:
                                trialIM_per_neuron[idx].append(im)

                if stim['listType'] == 1 and stim['stimLoc'] == 0:
                    centIndex = stim['centerIndex']
                    offIdx = stim['offsetIndex']
                    step = stim['stepIndex']
                    stimOnTimeS = ((1000 / frameRateHz * stim['stimOnFrame'])
                                   / 1000) + stim1TimeS
                    stimOffTimeS = ((1000 / frameRateHz * stim['stimOffFrame'])
                                    / 1000) + stim1TimeS
                    if step > numSteps:
                        indx = step - numSteps - 1
                        stCount = int(transectStimCount[indx])
                        transectStimCount[indx] += 1
                        stimIndex = transectCountIndex[indx]
                    else:
                        col = (int(offIdx != 0) * step * 2 - int(offIdx == 1))
                        # ^ basically indexing into stimCount stepIndex * 2 - 1
                        # for every column after the first, which is the blank
                        # condition
                        stCount = int(stimCount[centIndex, col])
                        stimCount[centIndex, col] += 1
                        stimIndex = stimCountIndex[centIndex, col]
                    if stCount > blocksDone:
                        break
                    if stimIndex <= 26 and catchTrial != 1:
                        # only units that are considered good for population analysis
                        for idx, unit in enumerate(sessionGoodUnits):
                            unitCount = np.where(units == unit)[0][0]
                            meanSpikes = meanSpikeReshaped[unitCount].reshape(27)
                            # normalized to pref center spike count
                            trialSpikes = ((spikeCountMat[unitCount, stCount, stimIndex] * 1000/trueStimDurMS)
                                           / meanSpikes[18])
                            # find Im
                            popt = paramFits[unitCount]
                            im = compute_im(stimIndex, trialSpikes, popt)
                            if abs(im) > 10:
                                trialIM_per_neuron[idx].append([])
                            else:
                                trialIM_per_neuron[idx].append(im)
            if catchTrial != 1:
                for im_list in trialIM_per_neuron:
                    globalTrialIM.append(im_list)
                targetOnTime = currTrial['taskEvents']['targetOn']['time']
                tempIndx = np.where((pupilAreaTimeS > stim1TimeS) & (pupilAreaTimeS < targetOnTime))
                trace = currTrial['eyeLPData']['data'][tempIndx]
                traceTime = 2 / 1000 * np.arange(len(trace))
                for _ in range(len(sessionGoodUnits)):
                    trialPAAlignedStim1Neuron.append(trace)
                    trialPAAlignedStim1TimeNeuron.append(traceTime)

    # to open another file in the loop
    os.chdir('../../Python Code')

print(time.time()-p0)

# ======================================================================================================================
# ================================================ Plotting Code =======================================================
# ======================================================================================================================

# ====== Initialize Lists to Arrays ======
masterImDelta = np.array(masterImDelta)
masterPupilDelta = np.array(masterPupilDelta)
popTrialStim1 = np.array(popTrialStim1)
popTrialStim2 = np.array(popTrialStim2)
popTrialStim3 = np.array(popTrialStim3)
popTrialStim4 = np.array(popTrialStim4)

# ====== Plot ======
for filler in range(1):
    x = np.array(masterImDelta)
    y = np.array(masterPupilDelta)
    # create a mask for valid (non-NaN) pairs
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_valid = x[mask]
    y_valid = y[mask]
    rho, pval = stats.spearmanr(x_valid, y_valid)

    plt.figure(figsize=(6, 5))
    # plt.scatter(delta_L0_blocks, delta_pupil_blocks, c='b')
    plt.scatter(masterImDelta, masterPupilDelta, c='b')
    plt.xlabel(r"$\Delta L_0$ (block - full)")
    plt.ylabel(r"$\Delta$ Pupil Area (block - mean)")
    plt.title(f"Block-wise ΔL0 vs Δ Pupil Area, spearman r = {rho:.3f}, p = {pval:.3e}")
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# to filter for outlier data points
for filler in range(1):
    mask = np.abs(masterImDelta) < 5
    filtered_L0 = masterImDelta[mask]
    filtered_pupil = masterPupilDelta[mask]
    # compute Spearman correlation
    rho, pval = stats.spearmanr(filtered_L0, filtered_pupil)

    plt.figure(figsize=(6, 5))
    plt.scatter(filtered_L0, filtered_pupil, c='b')
    plt.xlabel(r"$\Delta L_0$ (block - full)")
    plt.ylabel(r"$\Delta$ Pupil Area (block - mean)")
    plt.title(r"Block-wise $\Delta L_0$ vs $\Delta$ Pupil Area (|$\Delta L_0$| < 5)")
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.grid(True)
    # add Spearman rho & p-value as text box
    textstr = '\n'.join((
        r'Spearman $\rho$ = {:.3f}'.format(rho),
        r'$p$ = {:.3e}'.format(pval)))
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)

    plt.gca().text(0.05, 0.10, textstr, transform=plt.gca().transAxes,
                   fontsize=10, verticalalignment='top', bbox=props)
    plt.tight_layout()
    plt.show()


# visualize average pupil area trace (normalized x-axis)
for filler in range(1):
    # Example: set number of bins for normalized time
    n_bins = 100

    resampled_traces = []  # will hold resampled trials

    for trace, times in zip(trialPupilAreaList, trialPupilAreaTime):
        # Normalize time to [0,1]
        norm_times = (times - times[0]) / (times[-1] - times[0])
        # Define common normalized time base
        norm_time_bins = np.linspace(0, 1, n_bins)
        # Interpolate pupil trace onto the normalized time bins
        resampled = np.interp(norm_time_bins, norm_times, trace)
        resampled_traces.append(resampled)

    # Convert to a 2D array: (n_trials, n_bins)
    resampled_traces = np.array(resampled_traces)

    # Compute mean and SEM across trials
    mean_pupil = resampled_traces.mean(axis=0)
    sem_pupil = resampled_traces.std(axis=0) / np.sqrt(resampled_traces.shape[0])

    # Plot
    plt.figure(figsize=(6, 4))
    plt.plot(np.linspace(0, 1, n_bins), mean_pupil, label='Mean Pupil Area')
    plt.fill_between(np.linspace(0, 1, n_bins),
                     mean_pupil - sem_pupil,
                     mean_pupil + sem_pupil,
                     alpha=0.3, label='SEM')
    plt.xlabel('Normalized Trial Time')
    plt.ylabel('Pupil Area')
    plt.title('Average Pupil Area Across Trials')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# visualize the average pupil area trace with raw time x-axis & per-trial normalization
for filler in range(1):
    all_times = []
    all_pupils = []

    for trace, times in zip(trialPupilAreaList, trialPupilAreaTime):
        trace = np.array(trace)
        times = np.array(times)

        # already starts at 0, but convert to ms
        rel_times = times * 1000  # now in ms

        # normalize this trial’s pupil trace to its own max
        if np.max(trace) > 0:
            trace_norm = trace / np.max(trace)
        else:
            trace_norm = trace  # in case of all zeros (unlikely)

        all_times.append(rel_times)
        all_pupils.append(trace_norm)

    all_times = np.hstack(all_times)
    all_pupils = np.hstack(all_pupils)

    # Define bins
    bin_width = 20  # ms
    time_min = 0
    time_max = np.percentile(all_times, 99)  # ignore extreme long trials
    time_bins = np.arange(time_min, time_max + bin_width, bin_width)
    bin_centers = (time_bins[:-1] + time_bins[1:]) / 2

    mean_pupil = []
    sem_pupil = []

    for start, end in zip(time_bins[:-1], time_bins[1:]):
        mask = (all_times >= start) & (all_times < end)
        values = all_pupils[mask]
        if len(values) > 0:
            mean_pupil.append(values.mean())
            sem_pupil.append(values.std() / np.sqrt(len(values)))
        else:
            mean_pupil.append(np.nan)
            sem_pupil.append(np.nan)

    mean_pupil = np.array(mean_pupil)
    sem_pupil = np.array(sem_pupil)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(bin_centers / 1000, mean_pupil, label='Mean Pupil Area (normalized per trial)', color='b')
    plt.fill_between(bin_centers / 1000,
                     mean_pupil - sem_pupil,
                     mean_pupil + sem_pupil,
                     color='b', alpha=0.3, label='SEM')
    plt.xlabel('Time within Trial (s)')
    plt.ylabel('Normalized Pupil Area (max=1 per trial)')
    plt.title('Average Normalized Pupil Area Across Trials (Raw Time)')
    plt.xlim([0, time_max / 1000])  # limit x-axis
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# visualize the average pupil area trace with raw time x-axis & per-trial normalization
# aligned to stim 1
for filler in range(1):
    all_times = []
    all_pupils = []

    for trace, times in zip(trialPAAlignedStim1, trialPAAlignedStim1Time):
        trace = np.array(trace)
        times = np.array(times)

        # already starts at 0, but convert to ms
        rel_times = times * 1000  # now in ms

        # normalize this trial’s pupil trace to its own max
        if np.max(trace) > 0:
            trace_norm = trace / np.max(trace)
        else:
            trace_norm = trace  # in case of all zeros (unlikely)

        all_times.append(rel_times)
        all_pupils.append(trace_norm)

    all_times = np.hstack(all_times)
    all_pupils = np.hstack(all_pupils)

    # Define bins
    bin_width = 20  # ms
    time_min = 0
    time_max = np.percentile(all_times, 99)  # ignore extreme long trials
    time_bins = np.arange(time_min, time_max + bin_width, bin_width)
    bin_centers = (time_bins[:-1] + time_bins[1:]) / 2

    mean_pupil = []
    sem_pupil = []

    for start, end in zip(time_bins[:-1], time_bins[1:]):
        mask = (all_times >= start) & (all_times < end)
        values = all_pupils[mask]
        if len(values) > 0:
            mean_pupil.append(values.mean())
            sem_pupil.append(values.std() / np.sqrt(len(values)))
        else:
            mean_pupil.append(np.nan)
            sem_pupil.append(np.nan)

    mean_pupil = np.array(mean_pupil)
    sem_pupil = np.array(sem_pupil)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(bin_centers / 1000, mean_pupil, label='Mean Pupil Area (normalized per trial)', color='b')
    plt.fill_between(bin_centers / 1000,
                     mean_pupil - sem_pupil,
                     mean_pupil + sem_pupil,
                     color='b', alpha=0.3, label='SEM')
    plt.xlabel('Time within Trial (s)')
    plt.ylabel('Normalized Pupil Area (max=1 per trial)')
    plt.title('Average Normalized Pupil Area Across Trials Aligned to Stim 1 (Raw Time)')
    plt.xlim([0, time_max / 1000])  # limit x-axis
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

########################################################################################################################
# visualize the average pupil area trace with raw time x-axis & per-trial normalization
# aligned to stim 1 - to compare with how IM changes across a trial
for filler in range(1):
    all_times = []
    all_pupils = []

    for trace, times in zip(trialPAAlignedStim1Neuron, trialPAAlignedStim1TimeNeuron):
        trace = np.array(trace)
        times = np.array(times)

        # already starts at 0, but convert to ms
        rel_times = times * 1000  # now in ms

        # normalize this trial’s pupil trace to its own max
        if np.max(trace) > 0:
            trace_norm = trace / np.max(trace)
        else:
            trace_norm = trace  # in case of all zeros (unlikely)

        all_times.append(rel_times)
        all_pupils.append(trace_norm)

    all_times = np.hstack(all_times)
    all_pupils = np.hstack(all_pupils)

    # Define bins
    bin_width = 20  # ms
    time_min = 0
    time_max = np.percentile(all_times, 99.5)  # ignore extreme long trials
    time_bins = np.arange(time_min, time_max + bin_width, bin_width)
    bin_centers = (time_bins[:-1] + time_bins[1:]) / 2

    mean_pupil = []
    sem_pupil = []

    for start, end in zip(time_bins[:-1], time_bins[1:]):
        mask = (all_times >= start) & (all_times < end)
        values = all_pupils[mask]
        if len(values) > 0:
            mean_pupil.append(values.mean())
            sem_pupil.append(values.std() / np.sqrt(len(values)))
        else:
            mean_pupil.append(np.nan)
            sem_pupil.append(np.nan)

    mean_pupil = np.array(mean_pupil)
    sem_pupil = np.array(sem_pupil)

    # Plot
    plt.figure(figsize=(8, 5))

    plt.plot(bin_centers / 1000, mean_pupil, label='Mean Pupil Area (normalized per trial)', color='b')
    plt.fill_between(bin_centers / 1000,
                     mean_pupil - sem_pupil,
                     mean_pupil + sem_pupil,
                     color='b', alpha=0.3, label='SEM')

    plt.xlabel('Time within Trial (s)')
    plt.ylabel('Normalized Pupil Area (max=1 per trial)')
    plt.title('Average Normalized Pupil Area Across Trials Aligned to Stim 1 (Raw Time)')
    plt.xlim([0, time_max / 1000])  # limit x-axis
    plt.grid(True)
    plt.legend()

    # Add black bars at specified intervals
    bar_duration_ms = 250
    gap_duration_ms = 200
    total_cycle_ms = bar_duration_ms + gap_duration_ms

    # convert to seconds for plotting
    bar_duration_s = bar_duration_ms / 1000
    gap_duration_s = gap_duration_ms / 1000
    total_cycle_s = total_cycle_ms / 1000

    y_min, y_max = plt.ylim()
    bar_height = y_min - 0.05 * (y_max - y_min)  # slightly below min y
    bar_thickness = 0.02 * (y_max - y_min)  # height of the bar

    t = 0
    while t < (time_max / 1000):
        plt.axvspan(t, t + bar_duration_s,
                    ymin=0, ymax=0.02,  # small height above x-axis
                    color='k')
        t += total_cycle_s

    plt.tight_layout()
    plt.show()

# including position 1 Im
for filler in range(1):
    # step 1: gather values at each position
    position_dict = defaultdict(list)

    for sublist in globalTrialIM:
        for idx, val in enumerate(sublist):
            # skip if val is empty or not a number
            if isinstance(val, (int, float, np.integer, np.floating)):
                position_dict[idx].append(val)

    # step 2: compute stats
    positions = sorted(position_dict.keys())
    means, sems, ns = [], [], []

    for pos in positions:
        vals = np.array(position_dict[pos])
        if len(vals) == 0:
            means.append(np.nan)
            sems.append(np.nan)
            ns.append(0)
            continue

        means.append(np.mean(vals))
        sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)))  # SEM
        ns.append(len(vals))

    means = np.array(means)
    sems = np.array(sems)
    ns = np.array(ns)

    # step 3: print table
    print(f"{'Position':>10} {'N':>5} {'Mean':>10} {'SEM':>10}")
    for pos, m, s, n in zip(positions, means, sems, ns):
        print(f"{pos+1:>10} {n:>5} {m:>10.3f} {s:>10.3f}")

    # step 4: plot
    plt.figure(figsize=(8, 5))
    plt.errorbar(
        np.array(positions) + 1,
        means, yerr=sems,
        fmt='o-', capsize=4, color='blue', ecolor='gray', label='Mean ± SEM'
    )
    plt.xlabel('Position in Trial')
    plt.ylabel('IM')
    plt.title('Average IM at each Position (ignoring empty entries)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()

# excluding position 1 Im
for filler in range(1):
    # step 1: gather values at each position
    position_dict = defaultdict(list)

    for sublist in globalTrialIM:
        for idx, val in enumerate(sublist):
            if isinstance(val, (int, float, np.integer, np.floating)):
                position_dict[idx].append(val)

    # step 2: compute stats
    positions = sorted(position_dict.keys())  # starts at 0
    means, sems, ns = [], [], []

    for pos in positions:
        vals = np.array(position_dict[pos])
        if len(vals) == 0:
            means.append(np.nan)
            sems.append(np.nan)
            ns.append(0)
            continue

        means.append(np.median(vals))
        sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)))
        ns.append(len(vals))

    means = np.array(means)
    sems = np.array(sems)
    ns = np.array(ns)

    # 🔷 Insert empty at position 1 → shift everything right
    means = np.insert(means, 0, np.nan)
    sems = np.insert(sems, 0, np.nan)
    ns = np.insert(ns, 0, 0)

    positions = np.arange(1, len(means) + 1)

    # Print table
    print(f"{'Position':>10} {'N':>5} {'Mean':>10} {'SEM':>10}")
    for pos, m, s, n in zip(positions, means, sems, ns):
        print(f"{pos:>10} {n:>5} {m:>10.3f} {s:>10.3f}")

    # Plot
    plt.figure(figsize=(8, 5))
    plt.errorbar(
        positions,
        means,
        yerr=sems,
        fmt='o-', capsize=4, color='blue', ecolor='gray', label='Mean ± SEM'
    )

    plt.xticks(positions)  # ensure all positions incl. 1 are labeled
    plt.xlabel('Position in Trial')
    plt.ylabel('IM')
    plt.title('Average IM at each Position (Position 1 left empty)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()

# plot scatter of z-scored avg pupil area for each stimulus vs im for that stimulus

########################################################################################################################

# visualize the avg pupil area across blocks for each session
PAMeetz = pupilAreaAcrossBlockPop[:15]
PAAkshan = pupilAreaAcrossBlockPop[15:]
PABoth = pupilAreaAcrossBlockPop
for filler in range(1):
    max_blocks = max(len(session) for session in PABoth)
    # initialize matrix
    data_matrix = np.full((len(PABoth), max_blocks), np.nan)

    for i, session in enumerate(PABoth):
        session = np.array(session, dtype=float)
        session_norm = session / np.max(session)  # normalize to max=1
        data_matrix[i, :len(session_norm)] = session_norm

    # compute population mean & SEM ignoring NaNs
    mean_pupil = np.nanmean(data_matrix, axis=0)
    sem_pupil = np.nanstd(data_matrix, axis=0) / np.sqrt(np.sum(~np.isnan(data_matrix), axis=0))

    blocks = np.arange(1, max_blocks + 1)

    plt.figure(figsize=(10, 6))

    # Plot individual normalized sessions (faint lines)
    for i in range(data_matrix.shape[0]):
        plt.plot(blocks, data_matrix[i, :], alpha=0.2, color='gray', linewidth=1)

    # Plot population mean & SEM on top
    plt.plot(blocks, mean_pupil, label='Mean Normalized Pupil Area', color='blue', linewidth=2)
    plt.fill_between(blocks, mean_pupil - sem_pupil, mean_pupil + sem_pupil,
                     alpha=0.3, color='blue', label='SEM')

    plt.xlabel('Block Number')
    plt.ylabel('Normalized Pupil Area (max=1 per session)')
    plt.title('Normalized Pupil Area Across Blocks (Population + Individual Sessions (both))')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
# akshan last session is a bad one for pupil area


# visualize population mean-zscore response per block
for filler in range(1):
    # Find the maximum number of blocks
    max_blocks = max(len(neuron) for neuron in popBlockMeanZscore)

    # Initialize a matrix of nans
    data_matrix = np.full((len(popBlockMeanZscore), max_blocks), np.nan)

    # Fill matrix
    for i, neuron in enumerate(popBlockMeanZscore):
        data_matrix[i, :len(neuron)] = neuron

    # Compute mean & SEM at each block, ignoring NaNs
    mean_block = np.nanmean(data_matrix, axis=0)
    sem_block = np.nanstd(data_matrix, axis=0) / np.sqrt(np.sum(~np.isnan(data_matrix), axis=0))

    blocks = np.arange(1, max_blocks + 1)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(blocks, mean_block, color='b', label='Mean response')
    plt.fill_between(blocks, mean_block - sem_block, mean_block + sem_block,
                     color='b', alpha=0.3, label='SEM')
    plt.xlabel('Block number')
    plt.ylabel('Mean response (z-score)')
    plt.title('Population mean response per block')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# z-scored spike rate as function of stimulus position in trial, averaged across all trials
# will show whether stimulus position in the trial (first or last) changes the spike rate
for filler in range(1):
    # remove NaNs from each
    data = [
        popTrialStim1[~np.isnan(popTrialStim1)],
        popTrialStim2[~np.isnan(popTrialStim2)],
        popTrialStim3[~np.isnan(popTrialStim3)],
        popTrialStim4[~np.isnan(popTrialStim4)]
    ]
    labels = ['Stim 1', 'Stim 2', 'Stim 3', 'Stim 4']

    plt.figure(figsize=(8,6))
    plt.violinplot(data, showmeans=True, showmedians=True)

    plt.xticks(np.arange(1, 5), labels)
    plt.ylabel('Response (z-score)')
    plt.title('Distribution of Responses Across Stimuli (Violin Plots)')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

## attempt at single trial
# for unit in sessionGoodUnits:
#     i = np.where(unit == units)[0]
#     meanSpikes = meanSpikeReshaped[i].reshape(27)
#     # normalize to pref at center
#     meanSpikesNormalized = meanSpikes / meanSpikes[18]
#
#     pref = meanSpikesNormalized[[18, 2, 4, 6, 8]]
#     nonpref = meanSpikesNormalized[[9, 1, 3, 5, 7]]
#     prefNonpref = meanSpikesNormalized[[19, 21, 23, 25]]
#     nonprefPref = meanSpikesNormalized[[11, 13, 15, 17]]
#     pp = meanSpikesNormalized[[20, 22, 24, 26]]
#     nn = meanSpikesNormalized[[10, 12, 14, 16]]
#     spon = meanSpikesNormalized[0]
#
#     # fit heuristic model to average data from each condition to get parameters
#     resp = np.concatenate((pref, nonpref, prefNonpref, nonprefPref, pp, nn, [spon]), axis=0)
#     initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
#                      float(spon)]
#     popt, pcov = curve_fit(model_wrapperHeuristic, (contrast_centerSpon,
#                                                     contrast_peripherySpon,
#                                                     locationsSpon,
#                                                     stim_type_centerSpon,
#                                                     stim_type_peripherySpon), resp,
#                            p0=initial_guess)
#
#     # normalized all spike counts, normalized to mean pref center activity
#     normalizedSpikeCount = (spikeCountMat[i, :blocksDone, :] * 1000/trueStimDurMS) / meanSpikes[18]
#
#     individualIm = []
#     allPupilArea = []
#     meanPupilArea = np.mean(pupilAreaMat[i, :blocksDone, :])
#     # for single stim in periphery conditions
#     for cond in np.arange(1, 9):
#         group = ((cond - 1) // 2) + 1
#
#         # trials with single stim in offset
#         if group == 1:  # first offset
#             if cond % 2 == 1:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = ((popt[1] - trialSpikes) * popt[2]) / (trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 0:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = ((popt[0] - trialSpikes) * popt[2]) / (trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 2:  # second offset
#             if cond % 2 == 1:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = ((popt[1] - trialSpikes) * popt[3]) / (trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 0:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = ((popt[0] - trialSpikes) * popt[3]) / (trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 3:  # third offset
#             if cond % 2 == 1:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = ((popt[1] - trialSpikes) * popt[4]) / (trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 0:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = ((popt[0] - trialSpikes) * popt[4]) / (trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 4:  # fourth offset
#             if cond % 2 == 1:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = ((popt[1] - trialSpikes) * popt[5]) / (trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 0:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = ((popt[0] - trialSpikes) * popt[5]) / (trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#     # for single stim in center
#     for cond in [9, 18]:
#         if cond % 2 == 1:  # non-pref
#             for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                 # Compute z-score
#                 mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                 std_cond = np.std(normalizedSpikeCount[:, cond])
#                 z = (trialSpikes - mean_cond) / std_cond
#                 # Skip if > 3 SD away
#                 if abs(z) > 2:
#                     continue
#                 im = (popt[1] - trialSpikes) / (trialSpikes - popt[7])
#                 individualIm.append(im)
#                 allPupilArea.append(pupilAreaMat[i, count, cond])
#         if cond % 2 == 0:  # pref
#             for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                 # Compute z-score
#                 mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                 std_cond = np.std(normalizedSpikeCount[:, cond])
#                 z = (trialSpikes - mean_cond) / std_cond
#                 # Skip if > 3 SD away
#                 if abs(z) > 2:
#                     continue
#                 im = (popt[0] - trialSpikes) / (trialSpikes - popt[7])
#                 individualIm.append(im)
#                 allPupilArea.append(pupilAreaMat[i, count, cond])
#
#     # for paired stim n in center
#     for cond in np.arange(10, 18):
#         group = ((cond - 10) // 2) + 1
#
#         if group == 1:  # first offset
#             if cond % 2 == 0:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[1] - trialSpikes + (popt[2] * (popt[1] - trialSpikes))) / (
#                         trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 1:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[1] - trialSpikes + (popt[2] * (popt[0] - trialSpikes))) / (
#                             trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 1:  # second offset
#             if cond % 2 == 0:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[1] - trialSpikes + (popt[3] * (popt[1] - trialSpikes))) / (
#                         trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 1:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[1] - trialSpikes + (popt[3] * (popt[0] - trialSpikes))) / (
#                             trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 1:  # third offset
#             if cond % 2 == 0:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[1] - trialSpikes + (popt[4] * (popt[1] - trialSpikes))) / (
#                         trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 1:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[1] - trialSpikes + (popt[4] * (popt[0] - trialSpikes))) / (
#                             trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 1:  # fourth offset
#             if cond % 2 == 0:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[1] - trialSpikes + (popt[5] * (popt[1] - trialSpikes))) / (
#                         trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 1:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[1] - trialSpikes + (popt[5] * (popt[0] - trialSpikes))) / (
#                             trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#     # for paired stim p in center
#     for cond in np.arange(19, 27):
#         group = ((cond - 19) // 2) + 1
#
#         if group == 1:  # first offset
#             if cond % 2 == 1:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[0] - trialSpikes + (popt[2] * (popt[1] - trialSpikes))) / (
#                         trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 0:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[0] - trialSpikes + (popt[2] * (popt[0] - trialSpikes))) / (
#                             trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 1:  # second offset
#             if cond % 2 == 1:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[0] - trialSpikes + (popt[3] * (popt[1] - trialSpikes))) / (
#                         trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 0:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[0] - trialSpikes + (popt[3] * (popt[0] - trialSpikes))) / (
#                             trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 1:  # third offset
#             if cond % 2 == 1:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[0] - trialSpikes + (popt[4] * (popt[1] - trialSpikes))) / (
#                         trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 0:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[0] - trialSpikes + (popt[4] * (popt[0] - trialSpikes))) / (
#                             trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#         if group == 1:  # fourth offset
#             if cond % 2 == 1:  # non-pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[0] - trialSpikes + (popt[5] * (popt[1] - trialSpikes))) / (
#                         trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#             if cond % 2 == 0:  # pref
#                 for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#                     # Compute z-score
#                     mean_cond = np.mean(normalizedSpikeCount[:, cond])
#                     std_cond = np.std(normalizedSpikeCount[:, cond])
#                     z = (trialSpikes - mean_cond) / std_cond
#                     # Skip if > 3 SD away
#                     if abs(z) > 2:
#                         continue
#                     im = (popt[0] - trialSpikes + (popt[5] * (popt[0] - trialSpikes))) / (
#                             trialSpikes - popt[7])
#                     individualIm.append(im)
#                     allPupilArea.append(pupilAreaMat[i, count, cond])
#                     print(im)
#
#     individualIm = np.array(individualIm)
#     allPupilArea = np.array(allPupilArea)
#     deltaIM = individualIm - popt[6]
#     deltaPupilArea = allPupilArea - meanPupilArea
#
#     # get a boolean mask where |deltaIM| < 2
#     mask = np.abs(deltaIM) < 2
#
#     # use the mask to extract the valid points
#     deltaIM_filtered = deltaIM[mask]
#     deltaPupilArea_filtered = deltaPupilArea[mask]
#
#     # plot them
#     import matplotlib.pyplot as plt
#
#     plt.scatter(deltaIM_filtered, deltaPupilArea_filtered)
#     plt.xlabel('Delta IM')
#     plt.ylabel('Delta Pupil Area')
#     plt.title('|Delta IM| < 2')
#     plt.grid(True)
#     plt.show()
#
#     # pupil area plotting
#     plt.plot(np.mean((pupilAreaMat[i, :blocksDone, :]), axis=1))

# # trials with single stim in offset
# if 1 <= stimIndex < 9:
#     group = (stimIndex-1 // 2) + 1
#     if group == 1:  # first offset
#         if stimIndex % 2 == 1:  # non-pref
#             im = ((popt[1] - trialSpikes) * popt[2]) / (trialSpikes - popt[7])
#         if stimIndex % 2 == 0:  # pref
#             im = ((popt[0] - trialSpikes) * popt[2]) / (trialSpikes - popt[7])
#
#     if group == 2:  # second offset
#         if stimIndex % 2 == 1:  # non-pref
#             im = ((popt[1] - trialSpikes) * popt[3]) / (trialSpikes - popt[7])
#         if stimIndex % 2 == 0:  # pref
#             im = ((popt[0] - trialSpikes) * popt[3]) / (trialSpikes - popt[7])
#
#     if group == 3:  # third offset
#             im = ((popt[1] - trialSpikes) * popt[4]) / (trialSpikes - popt[7])
#         if stimIndex % 2 == 0:  # pref
#             im = ((popt[0] - trialSpikes) * popt[4]) / (trialSpikes - popt[7])
#
#     if group == 4:  # fourth offset
#         if stimIndex % 2 == 1:  # non-pref
#             im = ((popt[1] - trialSpikes) * popt[5]) / (trialSpikes - popt[7])
#         if stimIndex % 2 == 0:  # pref
#             im = ((popt[0] - trialSpikes) * popt[5]) / (trialSpikes - popt[7])
#
# # for single stim in center
# if stimIndex == 9:
#     im = (popt[1] - trialSpikes) / (trialSpikes - popt[7])
#
# if stimIndex == 18:
#     im = (popt[0] - trialSpikes) / (trialSpikes - popt[7])
#
# # for paired stim n in center
# if 10 <= stimIndex < 18:
#     group = ((stimIndex - 10) // 2) + 1
#     if group == 1:  # first offset
#         if stimIndex % 2 == 0:  # non-pref
#             im = (popt[1] - trialSpikes + (popt[2] * (popt[1] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#         if stimIndex % 2 == 1:  # pref
#             im = (popt[1] - trialSpikes + (popt[2] * (popt[0] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#
#     if group == 1:  # second offset
#         if stimIndex % 2 == 0:  # non-pref
#             im = (popt[1] - trialSpikes + (popt[3] * (popt[1] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#         if stimIndex % 2 == 1:  # pref
#             im = (popt[1] - trialSpikes + (popt[3] * (popt[0] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#
#     if group == 1:  # third offset
#         if stimIndex % 2 == 0:  # non-pref
#             im = (popt[1] - trialSpikes + (popt[4] * (popt[1] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#         if stimIndex % 2 == 1:  # pref
#             im = (popt[1] - trialSpikes + (popt[4] * (popt[0] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#
#     if group == 1:  # fourth offset
#         if stimIndex % 2 == 0:  # non-pref
#             im = (popt[1] - trialSpikes + (popt[5] * (popt[1] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#         if stimIndex % 2 == 1:  # pref
#             im = (popt[1] - trialSpikes + (popt[5] * (popt[0] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#
# # for paired stim p in center
# if 19 <= stimIndex < 27:
#     group = ((stimIndex - 19) // 2) + 1
#
#     if group == 1:  # first offset
#         if stimIndex % 2 == 1:  # non-pref
#             im = (popt[0] - trialSpikes + (popt[2] * (popt[1] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#         if stimIndex % 2 == 0:  # pref
#             im = (popt[0] - trialSpikes + (popt[2] * (popt[0] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#
#     if group == 1:  # second offset
#         if stimIndex % 2 == 1:  # non-pref
#             im = (popt[0] - trialSpikes + (popt[3] * (popt[1] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#         if stimIndex % 2 == 0:  # pref
#             im = (popt[0] - trialSpikes + (popt[3] * (popt[0] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#
#     if group == 1:  # third offset
#         if stimIndex % 2 == 1:  # non-pref
#             im = (popt[0] - trialSpikes + (popt[4] * (popt[1] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#         if stimIndex % 2 == 0:  # pref
#             for count, trialSpikes in enumerate(normalizedSpikeCount[:, cond]):
#             im = (popt[0] - trialSpikes + (popt[4] * (popt[0] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#
#     if group == 1:  # fourth offset
#         if stimIndex % 2 == 1:  # non-pref
#             im = (popt[0] - trialSpikes + (popt[5] * (popt[1] - trialSpikes))) / (
#                   trialSpikes - popt[7])
#         if stimIndex % 2 == 0:  # pref
#             im = (popt[0] - trialSpikes + (popt[5] * (popt[0] - trialSpikes))) / (
#                   trialSpikes - popt[7])
