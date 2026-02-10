# Imports
from usefulFns import *
from normalizationFunctions import *

p0 = time.time()


# # # Master list (both monkeys)
fileList = ['Meetz_230607', 'Meetz_230608', 'Meetz_230612', 'Meetz_230613', 'Meetz_230615',
            'Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628', 'Meetz_230630',
            'Meetz_230707', 'Meetz_230710', 'Meetz_230711', 'Meetz_230718', 'Meetz_230719',
            'Meetz_230720', 'Akshan_240530', 'Akshan_240603', 'Akshan_240606', 'Akshan_240607',
            'Akshan_240610', 'Akshan_240611', 'Akshan_240613', 'Akshan_240628', 'Akshan_240703',
            'Akshan_240704', 'Akshan_240705', 'Akshan_240709', 'Akshan_241014', 'Akshan_241024']

# # meetz only
# fileList = ['Meetz_230607', 'Meetz_230608', 'Meetz_230612', 'Meetz_230613', 'Meetz_230615',
#             'Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628', 'Meetz_230630',
#             'Meetz_230707', 'Meetz_230710', 'Meetz_230711', 'Meetz_230718', 'Meetz_230719',
#             'Meetz_230720']

# # akshan only
# fileList = ['Akshan_240530', 'Akshan_240603', 'Akshan_240606', 'Akshan_240607', 'Akshan_240610',
#             'Akshan_240611', 'Akshan_240613', 'Akshan_240628', 'Akshan_240703', 'Akshan_240704',
#             'Akshan_240705', 'Akshan_240709', 'Akshan_241014', 'Akshan_241024']


# # updated list with new cell IDs
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

# from old MTND files
# prefNonprefSwitchList = ['230627_89', '230627_110', '230627_114', '230719_54', '230719_124',
#                          '230719_150', '230615_153', '230615_151', '230718_199', '240606_65',
#                          '240606_69', '240606_96', '240606_97', '240606_101', '240607_10',
#                          '240607_76', '240610_22', '240628_47', '240703_63', '240703_69',
#                          '240703_58', '240703_65', '240705_113', '240705_49', '241014_19',
#                          '241014_32', '241014_59', '241014_24', '241014_57', '241024_0']

# from new MTND files with more units added from Phy
prefNonprefSwitchList = ['230627_89', '230627_110', '230627_114',
                         '230719_54', '230719_124', '230719_150', '230719_118', '230719_117',
                         '230615_153', '230615_151', '230615_53', '230615_43', '230615_61',
                         '230620_98',
                         '230626_22',
                         '230718_199',
                         '230607_19',
                         '230608_7', '230608_8', '230608_11', '230608_20', '230608_24',
                         '230608_32', '230608_33', '230608_41', '230608_42', '230608_43',
                         '240603_120', '240603_130',
                         '240606_65', '240606_69', '240606_96', '240606_97', '240606_101',
                         '240607_10', '240607_76'
                         '240610_22',
                         '240628_47',
                         '240703_63', '240703_69', '240703_58', '240703_65',
                         '240705_113', '240705_49',
                         '241014_19', '241014_32', '241014_59', '241014_24', '241014_57',
                         '241024_0', '241024_9', '241024_10']

unitListCopy = np.copy(unitList)
unitListCopy = np.array(unitListCopy)

# population arrays
prefNormalized = []
nonprefNormalized = []
pnNormalized = []
npNormalized = []
ppNormalized = []
nnNormalized = []
transectNormalized = []
sponNormalized = []
offsetDegSepNormPop = []
rfGaborSigmaPop = []
adaptationMat = np.zeros((10000, 4, 100))
adaptationMat[:] = np.nan
adaptC = 0
allBlocksDone = []
masterGoodUnits = []
masterAllUnits = []
totUnits = []
corrLists = [[] for _ in range(4)]
masterCorrUnits = []
masterCorrPairNum = []
masterCorrUnitsSpikeCountMat = []
masterCorrUnitsSpikeCountMatTransient = []
masterCorrUnitsSpikeCountMatSustained = []
popGaborSD = []
popGaborOffsetDeg = []
akshanFixWindow = []
meetzFixWindow = []

for file in fileList:
    # Load relevant file here with pyMat reader
    monkeyName, seshDate = file.split('_')

    fileName = f'{monkeyName}_{seshDate}_MTND_Spikes.mat'
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

    # list of indices of correctTrials (non-instruct, valid trialCertify)
    corrTrials = correctTrialsMTX(allTrials)

    # generate list of unique active units, and their channel
    units = activeUnits('spikeData', allTrials)
    totUnits.append(len(units))
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
    popGaborSD.append(allTrials[0]['rfGabor']['data']['sigmaDeg'])
    popGaborOffsetDeg.append(offsetDegSep)

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
    stimCount = np.zeros((3, numOffsets), dtype=int)
    stimCountIndex = np.arange((numSteps*7)+3-numSteps)
    stimCountIndex = stimCountIndex.reshape(3, int(len(stimCountIndex)/3))
    transectStimCount = np.zeros(numSteps, dtype=int)
    transectCountIndex = np.arange((numSteps*7)+3-numSteps, numSteps*7+3)
    onLatency = 50 / 1000  # time in MS for counting window latency after stim on
    offLatency = 50 / 1000  # time in MS for counting window latency after stim off
    histPrePostMS = 100  # 100ms window pre/post stimulus on/off
    spikeHists = np.zeros((len(units), numSteps*7+3, trueStimDurMS + (2*histPrePostMS+1)))

    # NEW: define the two analysis windows
    latencyWinMS = 50  # 50 ms wide: 50->100 ms after onset
    latencyEnd = onLatency + (latencyWinMS / 1000)  # 100 ms after onset
    sustainedStart = latencyEnd  # 100 ms after onset
    spikeCountMatTransient = np.zeros((len(units), blocksDone + 1, numSteps * 7 + 3))
    spikeCountMatSustained = np.zeros((len(units), blocksDone + 1, numSteps * 7 + 3))

    # insert spikes from valid stimulus presentations into spike count matrices
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        if 'spikeData' in currTrial:
            stimDesc = currTrial['stimDesc']['data']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            fixateTimeS = currTrial['taskEvents']['fixate']['time']
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

                    win_lat_start = stimOnTimeS + onLatency
                    win_lat_end = stimOnTimeS + latencyEnd  # stimOn + 100 ms
                    win_sus_start = stimOnTimeS + sustainedStart  # stimOn + 100 ms
                    win_sus_end = stimOffTimeS + offLatency

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

                            # NEW: 50–100 ms after onset
                            stimSpikes_lat = np.where((unitTimeStamps >= win_lat_start) &
                                                      (unitTimeStamps <= win_lat_end))[0]
                            spikeCountMatTransient[unitCount][stCount][stimIndex] = len(stimSpikes_lat)

                            # NEW: 100 ms after onset → stimOff + offLatency
                            stimSpikes_sus = np.where((unitTimeStamps >= win_sus_start) &
                                                      (unitTimeStamps <= win_sus_end))[0]
                            spikeCountMatSustained[unitCount][stCount][stimIndex] = len(stimSpikes_sus)

                            # PSTHs
                            stimOnPreSNEV = stimOnTimeS - (histPrePostMS / 1000)
                            stimOffPostSNEV = stimOffTimeS + (histPrePostMS / 1000)
                            histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                                             & (unitTimeStamps <= stimOffPostSNEV)
                                                             )] - stimOnPreSNEV
                            histStimSpikes = np.int32(histStimSpikes * 1000)
                            spikeHists[unitCount, stimIndex, histStimSpikes] += 1

    # reindex some spike count for some units because the non-pref direction is actually
    # the pref direction (this is just a.for correlation analysis)
    # start by flipping 2nd and 3rrd row and then flip column pairs for the stim index matrix
    stimIndexMat = np.arange(0, 27).reshape(3, 9)
    stimIndexMat[[1, 2]] = stimIndexMat[[2, 1]]
    for col in range(1, 9, 2):
        stimIndexMat[:, [col, col+1]] = stimIndexMat[:, [col+1, col]]
    order = np.asarray(stimIndexMat).reshape(-1)
    full_order = np.r_[order, np.arange(27, 31)]
    for idx, unit in enumerate(units):
        if f'{seshDate}_{unit}' in prefNonprefSwitchList:
            reordered = spikeCountMat[idx][:, full_order]
            spikeCountMat[idx] = reordered

    # mean, SEM, and reshaping of spikeCount matrices
    meanSpike = np.mean(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSD = np.std(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSEM = spikeCountSD / np.sqrt(blocksDone)
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

    # ================= TRANSIENT (50–100 ms) =================
    meanSpikeTransient = np.mean(spikeCountMatTransient[:, :blocksDone, :], axis=1)
    meanSpikeTransientReshaped = np.zeros((len(units), 3, int(((numSteps * 7) + 3 - numSteps) / 3)))
    for count, i in enumerate(meanSpikeTransientReshaped):
        meanSpikeTransientReshaped[count] = (
                meanSpikeTransient[count, :(numSteps * 7) + 3 - numSteps]
                .reshape(3, int(((numSteps * 7) + 3 - numSteps) / 3))
                * 1000 / 50  # 50 ms window
        )

    # ================= SUSTAINED (100 ms → offLatency) =================
    sustainedWindowMS = trueStimDurMS - 50  # matches your counting logic
    meanSpikeSustained = np.mean(spikeCountMatSustained[:, :blocksDone, :], axis=1)
    meanSpikeSustainedReshaped = np.zeros((len(units), 3, int(((numSteps * 7) + 3 - numSteps) / 3)))
    for count, i in enumerate(meanSpikeSustainedReshaped):
        meanSpikeSustainedReshaped[count] = (
                meanSpikeSustained[count, :(numSteps * 7) + 3 - numSteps]
                .reshape(3, int(((numSteps * 7) + 3 - numSteps) / 3))
                * 1000 / sustainedWindowMS
        )

    # add good units from each session to population arrays
    # for unit in sessionGoodUnits:
    corrUnits = []
    for unit in units:
        count = np.where(units == unit)[0][0]
        masterAllUnits.append(f'{seshDate}_{unit}')

        # transect
        x = np.concatenate((np.arange(-numSteps, 0),
                            np.arange(0, numSteps+1)),
                           axis=0)
        transReverse = transectMeanSpike[count, ::-1]
        transSEMReverse = transectSEM[count, ::-1]
        prefTransect = np.concatenate((transReverse,
                                       [meanSpikeReshaped[count, 2, 0]],
                                       meanSpikeReshaped[count, 0, 2::2]),
                                      axis=0)

        # NMI vs distance (Ni and Maunsell, 2017)
        pCentNMI = []
        nCentNMI = []
        ppNMI = []
        nnNMI = []

        for i in np.arange(1, numSteps+1):
            # gauss fit
            xTransect = np.concatenate((np.arange(-numSteps, 0),
                                       np.arange(0, numSteps + 1)),
                                       axis=0)
            params = gaussFit(xTransect, prefTransect)
            xFull = np.linspace(-numSteps, numSteps, 1000)
            respFull = gauss(xFull, *params)

            # PN
            centerP = meanSpikeReshaped[count, 2, 0]
            offsetN = meanSpikeReshaped[count, 0, i*2-1]
            PN = meanSpikeReshaped[count, 2, i*2-1]
            pnNMI = ((centerP + offsetN) - PN) / ((centerP + offsetN) + PN)
            pCentNMI.append(pnNMI)

            # NP
            centerN = meanSpikeReshaped[count, 1, 0]
            offsetP = meanSpikeReshaped[count, 0, i*2]
            NP = meanSpikeReshaped[count, 1, i*2]
            npNMI = ((centerN + offsetP) - NP) / ((centerN + offsetP) + NP)
            nCentNMI.append(npNMI)

            # PP
            PP = meanSpikeReshaped[count, 2, i*2]
            tempNMI = ((centerP + offsetP) - PP) / ((centerP + offsetP) + PP)
            ppNMI.append(tempNMI)

            # NN
            NN = meanSpikeReshaped[count, 1, i*2-1]
            tempNMI = ((centerN + offsetN) - NN) / ((centerN + offsetN) + NN)
            nnNMI.append(tempNMI)

        # normalize responses and add to population arrays
        # normVal = np.max(meanSpike[count]) * 1000 / trueStimDurMS  # np.max(prefTransect)
        # normVal = prefTransect[4]  # normalize to center response
        normVal = 1  # normalize everything outside the loop
        spon = meanSpikeReshaped[count, 0, 0] / normVal
        nullOnly = (np.concatenate(([meanSpikeReshaped[count, 1, 0]],
                                   meanSpikeReshaped[count, 0, 1::2]), axis=0) /
                    normVal)
        prefOnly = (np.concatenate(([meanSpikeReshaped[count, 2, 0]],
                                   meanSpikeReshaped[count, 0, 2::2]), axis=0) /
                    normVal)
        NP = meanSpikeReshaped[count, 1, 2::2] / normVal
        PN = meanSpikeReshaped[count, 2, 1::2] / normVal
        PP = meanSpikeReshaped[count, 2, 2::2] / normVal
        NN = meanSpikeReshaped[count, 1, 1::2] / normVal
        prefTransect = prefTransect / normVal

        # check to see if unit meets inclusion criterion, if so add to good units
        condition = True

        # check if pref resp is significantly greater than baseline
        prefCenterSpikeCounts = spikeCountMat[count][:blocksDone, 18]
        sponSpikeCounts = spikeCountMat[count][:blocksDone, 0]
        # Perform the one-sided Mann-Whitney U test (pref center > baseline)
        stat, p_value = mannwhitneyu(prefCenterSpikeCounts,
                                     sponSpikeCounts, alternative='greater')
        if p_value > 0.05:
            condition = False

        # check if nonpref resp is significantly greater than baseline
        nonprefCenterSpikeCounts = spikeCountMat[count][:blocksDone, 9]
        stat, p_value = mannwhitneyu(nonprefCenterSpikeCounts,
                                     sponSpikeCounts, alternative='greater')
        if p_value > 0.05:
            condition = False

        # check if pref resp is significantly greater than non-pref at center
        stat, p_value = mannwhitneyu(prefCenterSpikeCounts,
                                     nonprefCenterSpikeCounts, alternative='greater')
        if p_value > 0.05:
            condition = False

        # check if at least one other pref response is significantly greater than baseline
        otherPrefAboveBase = 0
        pref1SpikeCounts = spikeCountMat[count][:blocksDone, 20]
        pref2SpikeCounts = spikeCountMat[count][:blocksDone, 22]
        pref3SpikeCounts = spikeCountMat[count][:blocksDone, 24]
        pref4SpikeCounts = spikeCountMat[count][:blocksDone, 26]
        stat, p_value = mannwhitneyu(pref1SpikeCounts,
                                     sponSpikeCounts, alternative='greater')
        if p_value < 0.05:
            otherPrefAboveBase += 1
        stat, p_value = mannwhitneyu(pref2SpikeCounts,
                                     sponSpikeCounts, alternative='greater')
        if p_value < 0.05:
            otherPrefAboveBase += 1
        stat, p_value = mannwhitneyu(pref3SpikeCounts,
                                     sponSpikeCounts, alternative='greater')
        if p_value < 0.05:
            otherPrefAboveBase += 1
        stat, p_value = mannwhitneyu(pref4SpikeCounts,
                                     sponSpikeCounts, alternative='greater')
        if p_value < 0.05:
            otherPrefAboveBase += 1
        if otherPrefAboveBase < 1:
            condition = False

        # append good units to lists
        if condition:
            masterGoodUnits.append(f'{seshDate}_{unit}')
            sponNormalized.append(spon)
            prefNormalized.append(prefOnly)
            nonprefNormalized.append(nullOnly)
            pnNormalized.append(PN)
            npNormalized.append(NP)
            ppNormalized.append(PP)
            nnNormalized.append(NN)
            transectNormalized.append(prefTransect)
            offsetDegSepNormPop.append(offsetDegSepNorm)
            rfGaborSigmaPop.append(rfGaborSigma)
            allBlocksDone.append(blocksDone)

            # add preferred response at center, null at center, pref at first offset (including transect)
            # responses across all blocks from all valid neurons, to look at adaptation of response.
            p0Indx = stimCountIndex[2, 0]
            adaptationMat[adaptC, 0, :blocksDone] = (spikeCountMat[count, :blocksDone, p0Indx] /
                                                     np.max(spikeCountMat[count, :3, p0Indx]))
            n0Indx = stimCountIndex[1, 0]
            adaptationMat[adaptC, 1, :blocksDone] = (spikeCountMat[count, :blocksDone, n0Indx] /
                                                     np.max(spikeCountMat[count, :3, n0Indx]))
            p1Indx = stimCountIndex[0, 2]
            adaptationMat[adaptC, 2, :blocksDone] = (spikeCountMat[count, :blocksDone, p1Indx] /
                                                     np.max(spikeCountMat[count, :3, p1Indx]))
            pT1Indx = transectCountIndex[0]
            adaptationMat[adaptC, 3, :blocksDone] = (spikeCountMat[count, :blocksDone, pT1Indx] /
                                                     np.max(spikeCountMat[count, :3, pT1Indx]))
            adaptC += 1

            # add unit to be considered for correlation analysis
            corrUnits.append(unit)

    if len(corrUnits) > 1:
        # add units for correlation analysis to master list
        masterCorrPairNum.append(len(corrUnits))
        for unit in corrUnits:
            masterCorrUnits.append(f'{seshDate}_{unit}')

        # add the spike count mat for these units to a master list
        for uid in corrUnits:
            idx = np.where(units == uid)[0][0]
            masterCorrUnitsSpikeCountMat.append(spikeCountMat[idx, :blocksDone, :])
            masterCorrUnitsSpikeCountMatTransient.append(spikeCountMatTransient[idx, :blocksDone, :])
            masterCorrUnitsSpikeCountMatSustained.append(spikeCountMatSustained[idx, :blocksDone, :])

        corrUnits = np.array(corrUnits)
        # the different combinations of neuron pairs from total units
        combs = [i for i in combinations(corrUnits, 2)]
        for j in combs:
            n1 = np.where(units == j[0])[0][0]
            n2 = np.where(units == j[1])[0][0]
            for i in range(10, 18):  # center N offsets P or N
                # correlation for that Gabor pair b/w 2 units excluding trials where
                # spike counts exceeded 3 SD from mean
                n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
                pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                indx = (i - 10) // 2
                corrLists[indx].append(pairStimCorr)

            for i in range(19, 27):  # center P offsets P or N
                n1SpikeMat = spikeCountMat[n1, :blocksDone, i]
                n2SpikeMat = spikeCountMat[n2, :blocksDone, i]
                pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

                indx = (i - 19) // 2
                corrLists[indx].append(pairStimCorr)

    # to open another file in the loop
    os.chdir('../../Python Code')

print(time.time()-p0)

############# FIGURES AND FITTING ##################

# use this for fitting full RF weighted
numSteps = 4
numUnits = len(sponNormalized)
offsetDegSepNormPop = np.array(offsetDegSepNormPop)
rfGaborSigmaPop = np.array(rfGaborSigmaPop)
offsetDegSepNormPop = offsetDegSepNormPop.reshape((numUnits, numSteps*2+1))
sponNormalized = np.array(sponNormalized)
meanSpon = np.mean(sponNormalized)
semSpon = np.std(sponNormalized) / np.sqrt(numUnits)
masterGoodUnits = np.array(masterGoodUnits)
masterCorrUnits = np.array(masterCorrUnits)

prefNormalized = np.array(prefNormalized)
nonprefNormalized = np.array(nonprefNormalized)
pnNormalized = np.array(pnNormalized)
npNormalized = np.array(npNormalized)
ppNormalized = np.array(ppNormalized)
nnNormalized = np.array(nnNormalized)
transectNormalized = np.array(transectNormalized)

# normalize each response matrix by normalizing to
# the pref center response (keep baseline)
normVal = np.copy(prefNormalized[:, 0][:, np.newaxis])
prefNormalizedWithBase = prefNormalized / normVal
nonprefNormalizedWithBase = nonprefNormalized / normVal
pnNormalizedWithBase = pnNormalized / normVal
npNormalizedWithBase = npNormalized / normVal
ppNormalizedWithBase = ppNormalized / normVal
nnNormalizedWithBase = nnNormalized / normVal
sponNormalizedToPrefCenter = sponNormalized[:, np.newaxis] / normVal
normValBase = np.copy(normVal)
transectNormalized = transectNormalized / normValBase

# normalize each response matrix by subtracting baseline
# and then normalizing the pref center response
prefNormalized = prefNormalized - sponNormalized[:, np.newaxis]
normVal = np.copy(prefNormalized[:, 0][:, np.newaxis])
prefNormalized = prefNormalized / normVal
nonprefNormalized = nonprefNormalized - sponNormalized[:, np.newaxis]
nonprefNormalized = nonprefNormalized / normVal
pnNormalized = pnNormalized - sponNormalized[:, np.newaxis]
pnNormalized = pnNormalized / normVal
npNormalized = npNormalized - sponNormalized[:, np.newaxis]
npNormalized = npNormalized / normVal
ppNormalized = ppNormalized - sponNormalized[:, np.newaxis]
ppNormalized = ppNormalized / normVal
nnNormalized = nnNormalized - sponNormalized[:, np.newaxis]
nnNormalized = nnNormalized / normVal
# transectNormalized = transectNormalized - sponNormalized[:, np.newaxis]
# transectNormalized = transectNormalized / normVal

prefMean = np.mean(prefNormalized, axis=0)
nonprefMean = np.mean(nonprefNormalized, axis=0)
pnMean = np.mean(pnNormalized, axis=0)
npMean = np.mean(npNormalized, axis=0)
ppMean = np.mean(ppNormalized, axis=0)
nnMean = np.mean(nnNormalized, axis=0)
transectMean = np.mean(transectNormalized, axis=0)

prefSEM = np.std(prefNormalized, axis=0) / np.sqrt(numUnits)
nonprefSEM = np.std(nonprefNormalized, axis=0) / np.sqrt(numUnits)
pnSEM = np.std(pnNormalized, axis=0) / np.sqrt(numUnits)
npSEM = np.std(npNormalized, axis=0) / np.sqrt(numUnits)
ppSEM = np.std(ppNormalized, axis=0) / np.sqrt(numUnits)
nnSEM = np.std(nnNormalized, axis=0) / np.sqrt(numUnits)
transectSEM = np.std(transectNormalized, axis=0) / np.sqrt(numUnits)

# means of arrays with Baseline
prefMeanWithBase = np.mean(prefNormalizedWithBase, axis=0)
nonprefMeanWithBase = np.mean(nonprefNormalizedWithBase, axis=0)
pnMeanWithBase = np.mean(pnNormalizedWithBase, axis=0)
npMeanWithBase = np.mean(npNormalizedWithBase, axis=0)
ppMeanWithBase = np.mean(ppNormalizedWithBase, axis=0)
nnMeanWithBase = np.mean(nnNormalizedWithBase, axis=0)

prefWithBaseSEM = np.std(prefNormalizedWithBase, axis=0) / np.sqrt(numUnits)
nonprefWithBaseSEM = np.std(nonprefNormalizedWithBase, axis=0) / np.sqrt(numUnits)
pnWithBaseSEM = np.std(pnNormalizedWithBase, axis=0) / np.sqrt(numUnits)
npWithBaseSEM = np.std(npNormalizedWithBase, axis=0) / np.sqrt(numUnits)
ppWithBaseSEM = np.std(ppNormalizedWithBase, axis=0) / np.sqrt(numUnits)
nnWithBaseSEM = np.std(nnNormalizedWithBase, axis=0) / np.sqrt(numUnits)

## key variables
hfont = {'fontname': 'Arial'}
lineWidth = 1.0
mSize = 2
fontSize = 8
x = np.arange(0, numSteps + 1)
xNorm = np.arange(1, numSteps + 1)
transX = np.concatenate((np.arange(-numSteps, 0),
                         np.arange(0, numSteps + 1)),
                        axis=0)

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

########################################################################################################################
####################################### RF Weighted Normalization with baseline ########################################
########################################################################################################################

# individual neurons
r2ScoresRFWeightFullWithBase = []
aicScoresRFWeightFullWithBase = []
rfWeightNeuronFits = []
for i in range(len(prefNormalized)):
    resp = np.concatenate((prefNormalizedWithBase[i], nonprefNormalizedWithBase[i],
                           pnNormalizedWithBase[i], npNormalizedWithBase[i],
                           ppNormalizedWithBase[i], nnNormalizedWithBase[i],
                           sponNormalizedToPrefCenter[i]), axis=0)
    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, 0.10, 0.10,
                     float(sponNormalizedToPrefCenter[i])]  # Lp, Lnp, W1, W2, W3, W4, sigma, b

    # weighted norm with baseline
    popt, pcov = curve_fit(model_wrapperWithBase, (contrast_centerSpon,
                                                   contrast_peripherySpon,
                                                   locationsSpon,
                                                   stim_type_centerSpon,
                                                   stim_type_peripherySpon), resp,
                           p0=initial_guess)
    y_pred1 = apply_fitted_modelWithBase(contrast_centerSpon, contrast_peripherySpon,
                                         locationsSpon, stim_type_centerSpon,
                                         stim_type_peripherySpon, *popt)
    r2 = r2_score(resp, y_pred1)
    r2ScoresRFWeightFullWithBase.append(r2)
    rfWeightNeuronFits.append(popt)

    # compute AIC and append to list
    residuals = resp - y_pred1
    rss = np.sum(residuals**2)
    n = len(resp)         # number of data points
    k = len(popt)         # number of parameters, should be 8 here

    aic = n * np.log(rss / n) + 2 * k
    aicScoresRFWeightFullWithBase.append(aic)

    # single presentation prediction
    pCenterSinglePred = np.array(fullRFWeightedNormWithBase(1, 0, -1, 1, -1, *popt))
    pCenterSinglePred = pCenterSinglePred[np.newaxis]
    pPeriSinglePred = np.array([fullRFWeightedNormWithBase(0, 1, j, -1, 1, *popt) for j in range(4)])
    pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

    npCenterSinglePred = np.array(fullRFWeightedNormWithBase(1, 0, -1, 0, -1, *popt))
    npCenterSinglePred = npCenterSinglePred[np.newaxis]
    npPeriSinglePred = np.array([fullRFWeightedNormWithBase(0, 1, j, -1, 0, *popt) for j in range(4)])
    npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

    # paired presentation prediction
    pnFitPred = [fullRFWeightedNormWithBase(1, 1, j, 1, 0, *popt) for j in range(4)]
    npFitPred = [fullRFWeightedNormWithBase(1, 1, j, 0, 1, *popt) for j in range(4)]
    ppFitPred = [fullRFWeightedNormWithBase(1, 1, j, 1, 1, *popt) for j in range(4)]
    nnFitPred = [fullRFWeightedNormWithBase(1, 1, j, 0, 0, *popt) for j in range(4)]

    predResp = np.concatenate((pSingleFitPred, npSingleFitPred, pnFitPred,
                               npFitPred, ppFitPred, nnFitPred), axis=0)

r2ScoresRFWeightFullWithBase = np.array(r2ScoresRFWeightFullWithBase)
np.median(r2ScoresRFWeightFullWithBase)

########################################################################################################################
################################## RF Weighted Normalization (Intensity Weighted) ######################################
############################################ Heuristic Model with ib ###################################################
########################################################################################################################

# individual neurons
r2ScoresRFWeightHeuristic = []
rfWeightNeuronFits = []

for i in range(len(prefNormalized)):

    resp = np.concatenate((prefNormalizedWithBase[i], nonprefNormalizedWithBase[i],
                           pnNormalizedWithBase[i], npNormalizedWithBase[i],
                           ppNormalizedWithBase[i], nnNormalizedWithBase[i],
                           sponNormalizedToPrefCenter[i]), axis=0)
    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                     float(sponNormalizedToPrefCenter[i])]

    # classic norm with baseline
    popt, pcov = curve_fit(model_wrapperHeuristic, (contrast_centerSpon,
                                                    contrast_peripherySpon,
                                                    locationsSpon,
                                                    stim_type_centerSpon,
                                                    stim_type_peripherySpon), resp,
                           p0=initial_guess)
    y_pred1 = apply_fitted_modelHeuristic(contrast_centerSpon, contrast_peripherySpon,
                                          locationsSpon, stim_type_centerSpon,
                                          stim_type_peripherySpon, *popt)
    r2 = r2_score(resp, y_pred1)
    r2ScoresRFWeightHeuristic.append(r2)
    rfWeightNeuronFits.append(popt)

r2ScoresRFWeightHeuristic = np.array(r2ScoresRFWeightHeuristic)
np.median(r2ScoresRFWeightHeuristic)


########################################################################################################################
###################################### Direction Tuning and Correlation Analysis #######################################
########################################################################################################################


# ========================= Von Mises model & helpers =========================
def von_mises_func(theta, baseline, amplitude, pref_dir, kappa):
    """
    theta, pref_dir in radians
    normalized von Mises bump + baseline
    """
    return baseline + amplitude * np.exp(kappa * np.cos(theta - pref_dir)) / (2 * np.pi * i0(kappa))


def r2_score_safe(y, yhat):
    """Robust R^2 that handles NaNs and near-constant vectors."""
    m = np.isfinite(y) & np.isfinite(yhat)
    if m.sum() < 2:
        return np.nan
    y_ = y[m]; yh = yhat[m]
    ss_res = np.sum((y_ - yh)**2)
    ss_tot = np.sum((y_ - np.mean(y_))**2)
    if ss_tot <= 1e-12:
        # If data are (near) constant, define R^2 sensibly
        return 1.0 if ss_res <= 1e-12 else 0.0
    return 1.0 - ss_res/ss_tot


def fit_von_mises_single(responses, theta_rad):
    """
    Fit a single neuron's direction tuning with baseline + amplitude * vonMises.
    Returns: baseline, amplitude, pref_dir(rad in [0,2π)), kappa
    """
    responses = np.asarray(responses, dtype=float)
    theta_rad = np.asarray(theta_rad, dtype=float)

    # Initial guesses
    baseline_guess = np.nanmin(responses)
    amp_guess = max(np.nanmax(responses) - baseline_guess, 0.0)
    pref_guess = theta_rad[np.nanargmax(responses)]
    kappa_guess = 2.0

    p0 = [baseline_guess, amp_guess, pref_guess, kappa_guess]

    # Bounds: baseline free, amplitude >= 0, pref_dir in [0, 2π], kappa >= 1e-6
    bounds = (
        [-np.inf, 0.0, 0.0, 1e-6],
        [ np.inf, np.inf, 2*np.pi, 1000.0]
    )

    # Fit
    try:
        params, _ = curve_fit(
            von_mises_func, theta_rad, responses, p0=p0,
            bounds=bounds, maxfev=20000
        )
        baseline, amplitude, pref_dir, kappa = params
        pref_dir = np.mod(pref_dir, 2*np.pi)  # normalize to [0, 2π)
        return baseline, amplitude, pref_dir, kappa
    except Exception:
        # Fallback if fitting fails
        return np.nan, np.nan, np.nan, np.nan


# ========================= Build session labels =========================
# Assumes: unitList is defined outside as a list like ['230607_77', ...]
dates = [u.split('_')[0] for u in masterGoodUnits]
date_counts = Counter(dates)
multi_unit_dates = [d for d, count in date_counts.items() if count > 1]

date_labels = []
for d in multi_unit_dates:
    if d.startswith("23"):
        date_labels.append(f"Meetz_{d}")
    elif d.startswith("24"):
        date_labels.append(f"Akshan_{d}")
    else:
        date_labels.append(d)  # fallback if another year shows up

# ========================= Iterate sessions & collect similarity =========================
base_dir = "/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN"

similarity_records = []  # will accumulate dicts across all sessions
for date in date_labels:
    monkeyName, sessionDate = date.split('_')
    sesh_dir = os.path.join(base_dir, monkeyName, sessionDate, "Direction Tuning")
    os.chdir(sesh_dir)

    # shape: (n_units_today, n_dirs)
    dirData = np.load('unitsDirTuningMat.npy')

    # Units present today (strings like '230607_77'), and the good ones
    # Assumes: masterAllUnits and masterGoodUnits exist in your environment
    units_today = [u for u in masterAllUnits if u.startswith(sessionDate + "_")]
    goodUnits_today = [u for u in masterGoodUnits if u.startswith(sessionDate + "_")]

    # Map unit string -> row index in dirData
    # ASSUMPTION: rows of dirData correspond to units_today in order.
    unit_to_row = {u: i for i, u in enumerate(units_today)}

    # Keep only good units that actually exist in today's data
    goodUnits_today = [u for u in goodUnits_today if u in unit_to_row]

    # Need at least a pair
    if len(goodUnits_today) < 2:
        os.chdir(base_dir)
        continue

    # Angle grid
    n_neurons_today, n_dirs = dirData.shape
    angles_deg = np.linspace(0, 360, n_dirs, endpoint=False)
    angles_rad = np.deg2rad(angles_deg)

    # ----- Fit von Mises for good units; store PDs and R^2 aligned to goodUnits_today -----
    pref_dirs_rad = np.zeros(len(goodUnits_today), dtype=float)
    r2_list = np.full(len(goodUnits_today), np.nan, dtype=float)

    for gi, u in enumerate(goodUnits_today):
        idx = unit_to_row[u]
        yi = dirData[idx, :].astype(float)
        b, A, prefD, k = fit_von_mises_single(yi, angles_rad)
        pref_dirs_rad[gi] = prefD

        # Predicted tuning & R^2 (NaN-safe)
        if np.isfinite(b) and np.isfinite(A) and np.isfinite(prefD) and np.isfinite(k):
            yhat = von_mises_func(angles_rad, b, A, prefD, k)
            r2_list[gi] = r2_score_safe(yi, yhat)
        else:
            r2_list[gi] = np.nan

    pref_dirs_deg = np.rad2deg(pref_dirs_rad)

    # ----- Tuning similarity (cos Δθ): SI in [-1,1], SI01 in [0,1] -----
    c = np.cos(pref_dirs_rad).reshape(-1, 1)
    s = np.sin(pref_dirs_rad).reshape(-1, 1)
    SI = (c @ c.T) + (s @ s.T)     # (ngood x ngood), includes self-diagonal=1
    SI01 = 0.5 * (SI + 1.0)

    # ----- Bhattacharyya on raw tuning curves -----
    eps = 1e-12
    pmfs = []
    for u in goodUnits_today:
        idx = unit_to_row[u]
        y = dirData[idx, :].astype(float)
        y = np.clip(y, 0.0, None)                     # ensure nonnegative
        p = (y + eps) / (np.sum(y + eps))             # normalize to sum=1
        pmfs.append(p)
    pmfs = np.vstack(pmfs)

    # ----- Record only unique pairs (i < j) -----
    for i, j in combinations(range(len(goodUnits_today)), 2):
        # Bhattacharyya coefficient & distance
        bc = float(np.sum(np.sqrt(pmfs[i] * pmfs[j])))    # in [0,1]
        bd = float(-np.log(max(bc, eps)))                 # >= 0

        similarity_records.append({
            "date": sessionDate,
            "monkey": monkeyName,
            "unit_i": goodUnits_today[i],
            "unit_j": goodUnits_today[j],
            "PD_i_deg": float(pref_dirs_deg[i]),
            "PD_j_deg": float(pref_dirs_deg[j]),
            "SI": float(SI[i, j]),
            "SI01": float(SI01[i, j]),
            "BC_dir_tuning": bc,      # Bhattacharyya coefficient (0..1)
            "BDir": bd,               # Bhattacharyya distance (-log(BC))
            "R2_i": float(r2_list[i]),  # NEW: von Mises fit R^2 for unit_i
            "R2_j": float(r2_list[j])   # NEW: von Mises fit R^2 for unit_j
        })

    # Return to base directory before next session
    os.chdir(base_dir)

########################################################################################################################
#################################################### Correlation Analysis ##############################################
########################################################################################################################


def build_r2_lookup(similarity_records):
    """
    Map an unordered unit pair -> (R2_i, R2_j) taken from similarity_records.
    Orientation doesn't matter since we only require both >= R2_MIN.
    """
    r2_lookup = {}
    for r in similarity_records:
        a, b = r["unit_i"], r["unit_j"]
        key = tuple(sorted((a, b)))
        r2_i = r.get("R2_i", np.nan)
        r2_j = r.get("R2_j", np.nan)
        r2_lookup[key] = (r2_i, r2_j)
    return r2_lookup


def shared_geom_share_rectified(w1, w2, sigma1, sigma2):
    eps = 1e-9
    w1p = max(w1, 0.0)
    w2p = max(w2, 0.0)
    s1 = (1 + w1p) / (1 + w1p + sigma1 + eps)   # in [0,1)
    s2 = (1 + w2p) / (1 + w2p + sigma2 + eps)
    return np.sqrt(s1 * s2)


# --- circular |ΔPD| helper ---
def circ_delta_abs_deg(a_deg, b_deg):
    """Absolute circular difference between two angles in degrees, in [0, 180]."""
    d = np.abs((a_deg - b_deg) % 360.0)
    return float(np.minimum(d, 360.0 - d))


# Lookups from similarity_records
def build_similarity_lookup(similarity_records):
    sim_SI = {}
    sim_SI01 = {}
    for r in similarity_records:
        a, b = r["unit_i"], r["unit_j"]
        key = tuple(sorted((a, b)))
        sim_SI[key] = r.get("SI")
        sim_SI01[key] = r.get("SI01")
    return sim_SI, sim_SI01


def build_pd_lookup(similarity_records):
    pd_lookup = {}
    for r in similarity_records:
        a, b = r["unit_i"], r["unit_j"]
        key = tuple(sorted((a, b)))
        pd_lookup[key] = (r.get("PD_i_deg"), r.get("PD_j_deg"))
    return pd_lookup


def build_monkey_lookup(similarity_records):
    monkey_lookup = {}
    for r in similarity_records:
        a, b = r["unit_i"], r["unit_j"]
        key = tuple(sorted((a, b)))
        monkey_lookup[key] = r.get("monkey")
    return monkey_lookup


# === NEW: BC/BD lookup from similarity_records (order-agnostic) ===
def build_bcbd_lookup(similarity_records):
    bc_lookup = {}
    bd_lookup = {}
    for r in similarity_records:
        a, b = r["unit_i"], r["unit_j"]
        key = tuple(sorted((a, b)))
        bc_lookup[key] = r.get("BC_dir_tuning")  # may be None if not computed
        bd_lookup[key] = r.get("BDir")           # may be None if not computed
    return bc_lookup, bd_lookup


# === map paired index i to single indices (center/periphery) and site ===
def single_indices_for_pair(i):
    """
    Return (center_single_idx, periph_single_idx, site)
    - For i in 10..17 (center=N paired): center_single_idx=9
    - For i in 19..26 (center=P paired): center_single_idx=18
    - periph_single_idx alternates N/P by within-site parity
      and yields {1..8} depending on site and parity.
    """
    if 10 <= i <= 17:  # center=N paired block
        site = (i - 10) // 2            # 0..3
        within = (i - 10) % 2           # 0 -> periph N, 1 -> periph P
        center_single_idx = 9
        periph_single_idx = 2*site + 1 + within  # e.g., i=10->(9,1), i=11->(9,2)
        return center_single_idx, periph_single_idx, site
    elif 19 <= i <= 26:  # center=P paired block
        site = (i - 19) // 2            # 0..3
        within = (i - 19) % 2           # 0 -> periph N, 1 -> periph P
        center_single_idx = 18
        periph_single_idx = 2*site + 1 + within  # e.g., i=25->(18,7), i=26->(18,8)
        return center_single_idx, periph_single_idx, site
    else:
        raise ValueError(f"i={i} is not in a paired block (10..17 or 19..26).")


# Helper to add a record (now includes each neuron's site-specific weight)
def add_record(block, i, comb, pairStimCorr, pairWeightGeomMean, n1Sel, n2Sel,
               sim_SI=None, sim_SI01=None, delta_PD_abs_deg=None, avg_pair_weight_2w=None,
               monkey=None, singleCorr_center=None, singleCorr_periph=None,
               pairStimCorr_transient=None,
               singleCorr_center_transient=None,
               singleCorr_periph_transient=None,
               pairStimCorr_sustained=None,
               singleCorr_center_sustained=None,
               singleCorr_periph_sustained=None,
               singleWeightGeom_center=None, singleWeightGeom_periph=None,
               BC_dir_tuning=None, BDir=None,
               n1w=None, n2w=None,
               pairMinShareGeom=None, n1MinShare=None, n2MinShare=None,
               n1Sig=None, n2Sig=None, pairCov=None, singleCCov=None,
               singlePCov=None,
               pairMeanFR_paired_hz=None,
               pairMeanFR_single_center_hz=None,
               pairMeanFR_single_periph_hz=None):
    u1, u2 = comb
    date = u1.split("_")[0]  # both units should share the same date
    corr_records.append({
        "date": date,
        "unit_i": u1,
        "unit_j": u2,
        "block": block,              # e.g., "centerN_offsetsPN" or "centerP_offsetsPN"
        "cond_index": i,             # 10..17 or 19..26
        "corr": float(pairStimCorr),
        # NEW: paired corr (transient / sustained)
        "corr_transient": (None if pairStimCorr_transient is None else float(pairStimCorr_transient)),
        "corr_sustained": (None if pairStimCorr_sustained is None else float(pairStimCorr_sustained)),
        # NEW singles (transient / sustained)
        "singleCorr_center_transient": (
            None if singleCorr_center_transient is None else float(singleCorr_center_transient)),
        "singleCorr_periph_transient": (
            None if singleCorr_periph_transient is None else float(singleCorr_periph_transient)),
        "singleCorr_center_sustained": (
            None if singleCorr_center_sustained is None else float(singleCorr_center_sustained)),
        "singleCorr_periph_sustained": (
            None if singleCorr_periph_sustained is None else float(singleCorr_periph_sustained)),
        "pair_weight_geommean": float(pairWeightGeomMean),
        "n1Sel": float(n1Sel),
        "n2Sel": float(n2Sel),
        "pairSel_geom": float(np.sign(n1Sel)*np.sign(n2Sel)*np.sqrt(abs(n1Sel)*abs(n2Sel))),
        "similarity_SI": (None if sim_SI is None else float(sim_SI)),
        "similarity_SI01": (None if sim_SI01 is None else float(sim_SI01)),
        "delta_PD_abs_deg": (None if delta_PD_abs_deg is None else float(delta_PD_abs_deg)),
        "avg_pair_weight_2w": (None if avg_pair_weight_2w is None else float(avg_pair_weight_2w)),
        "monkey": (None if monkey is None else str(monkey)),
        # Singles:
        "singleCorr_center": (None if singleCorr_center is None else float(singleCorr_center)),
        "singleCorr_periph": (None if singleCorr_periph is None else float(singleCorr_periph)),
        "singleWeightGeom_center": (None if singleWeightGeom_center is None else float(singleWeightGeom_center)),
        "singleWeightGeom_periph": (None if singleWeightGeom_periph is None else float(singleWeightGeom_periph)),
        # Bhattacharyya metrics
        "BC_dir_tuning": (None if BC_dir_tuning is None else float(BC_dir_tuning)),
        "BDir": (None if BDir is None else float(BDir)),
        # Per-neuron weight used in this paired condition
        "n1w": (None if n1w is None else float(n1w)),
        "n2w": (None if n2w is None else float(n2w)),
        # NEW: min-share metric
        "pairMinShareGeom": (None if pairMinShareGeom is None else float(pairMinShareGeom)),
        "n1MinShare": (None if n1MinShare is None else float(n1MinShare)),
        "n2MinShare": (None if n2MinShare is None else float(n2MinShare)),
        "n1Sigma": (None if n1Sig is None else float(n1Sig)),
        "n2Sigma": (None if n2Sig is None else float(n2Sig)),
        "pairCov": (None if pairCov is None else float(pairCov)),
        "singleCenterCov": (None if singleCCov is None else float(singleCCov)),
        "singlePeriphCov": (None if singlePCov is None else float(singlePCov)),
        "pairMeanFR_paired_hz": (None if pairMeanFR_paired_hz is None else float(pairMeanFR_paired_hz)),
        "pairMeanFR_single_center_hz": (
            None if pairMeanFR_single_center_hz is None else float(pairMeanFR_single_center_hz)),
        "pairMeanFR_single_periph_hz": (
            None if pairMeanFR_single_periph_hz is None else float(pairMeanFR_single_periph_hz)),
    })

# stim index mat (reference)
stimIndexMat = np.arange(0, 27).reshape(3, 9)

# ---- pairs to evaluate ----
all_combs = []
start = 0
for n in masterCorrPairNum:
    day_units = masterCorrUnits[start:start + n]
    all_combs.extend(combinations(day_units, 2))
    start += n

# index (exclusive) where Monkey 1’s combos end in all_combs
split_idx = sum(n*(n-1)//2 for n in masterCorrPairNum[:9])
monkey1_combs = all_combs[:split_idx]
monkey2_combs = all_combs[split_idx:]

# Lookups
sim_SI_lookup, sim_SI01_lookup = build_similarity_lookup(similarity_records)
pd_lookup = build_pd_lookup(similarity_records)
monkey_lookup = build_monkey_lookup(similarity_records)
bc_lookup, bd_lookup = build_bcbd_lookup(similarity_records)
r2_lookup = build_r2_lookup(similarity_records)

# Faster index maps for master arrays to avoid repeated np.where calls
good_index = {u: i for i, u in enumerate(masterGoodUnits)}
corr_index = {u: i for i, u in enumerate(masterCorrUnits)}

# ===== R^2 gate (toggle this as you like) =====
R2_MIN = 0.0  # e.g., 0.80 means 80%; set to 0.70, 0.90, etc.
# Accumulators
corr_records = []   # list of dict rows
masterCorr = []
masterPairWeight = []
masterPairSel = []
sameSelCorr = []
sameSelWeight = []
oppositeSelCorr = []
oppositeSelWeight = []
HZ_PER_COUNT = 1000.0 / trueStimDurMS

for comb in all_combs:
    if comb[0] not in good_index or comb[1] not in good_index:
        continue

    n1 = good_index[comb[0]]
    n2 = good_index[comb[1]]

    # similarity/PD/BC lookups (order-agnostic)
    key = tuple(sorted((comb[0], comb[1])))
    monkey_name = monkey_lookup.get(key)
    sim_val = sim_SI_lookup.get(key)
    sim01_val = sim_SI01_lookup.get(key)
    bc_val = bc_lookup.get(key)   # may be None if not in similarity_records
    bd_val = bd_lookup.get(key)

    # R² gate for this pair
    r2_i, r2_j = r2_lookup.get(key, (np.nan, np.nan))
    r2_ok = (np.isfinite(r2_i) and np.isfinite(r2_j) and
             (r2_i >= R2_MIN) and (r2_j >= R2_MIN))

    PD_i_deg, PD_j_deg = (None, None)
    if key in pd_lookup:
        PD_i_deg, PD_j_deg = pd_lookup[key]
    deltaPD = None
    if (PD_i_deg is not None) and (PD_j_deg is not None):
        deltaPD = circ_delta_abs_deg(PD_i_deg, PD_j_deg)

    # =========================
    # Block 1: i = 10..17 (center=N offsets P or N)
    # =========================
    for i in range(10, 18):
        n1SpkMatIdx = corr_index[comb[0]]
        n2SpkMatIdx = corr_index[comb[1]]
        n1SpikeMat = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, i]
        n2SpikeMat = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, i]
        pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat,
                                                             method='pearson')

        # --- TRANSIENT (50–100ms) paired correlation ---
        n1SpikeMat_T = masterCorrUnitsSpikeCountMatTransient[n1SpkMatIdx][:, i]
        n2SpikeMat_T = masterCorrUnitsSpikeCountMatTransient[n2SpkMatIdx][:, i]
        pairStimCorr_T, _, _ = pairCorrExclude3SD(n1SpikeMat_T, n2SpikeMat_T, method='pearson')

        # --- SUSTAINED (100ms → offLatency) paired correlation ---
        n1SpikeMat_S = masterCorrUnitsSpikeCountMatSustained[n1SpkMatIdx][:, i]
        n2SpikeMat_S = masterCorrUnitsSpikeCountMatSustained[n2SpkMatIdx][:, i]
        pairStimCorr_S, _, _ = pairCorrExclude3SD(n1SpikeMat_S, n2SpikeMat_S, method='pearson')

        # --- PAIRED avg FR across the two neurons (Hz) ---
        # Each is mean spike count over trials -> convert to Hz -> average across neurons
        n1_mean_spk_paired = float(np.mean(n1SpikeMat))
        n2_mean_spk_paired = float(np.mean(n2SpikeMat))
        n1_fr_paired_hz = n1_mean_spk_paired * HZ_PER_COUNT
        n2_fr_paired_hz = n2_mean_spk_paired * HZ_PER_COUNT
        pairMeanFR_paired_hz = 0.5 * (n1_fr_paired_hz + n2_fr_paired_hz)

        # Correct single indices and site
        center_single_idx, periph_single_idx, site = single_indices_for_pair(i)

        # Weights / sigmas for this site
        w0 = 1.0
        n1w = rfWeightNeuronFits[n1][2 + site]
        sig1 = rfWeightNeuronFits[n1][6]
        n2w = rfWeightNeuronFits[n2][2 + site]
        sig2 = rfWeightNeuronFits[n2][6]
        if n1w <= 0 or n2w <= 0:
            continue

        # --- NEW: min-share per neuron and geometric mean across pair ---
        # For each neuron: share = min(w0, w1) / (w0 + w1 + sigma); here w0=1, w1 = that neuron's site weight
        n1MinShare = min(1.0, float(n1w)) / (1.0 + float(n1w) + float(sig1))
        n2MinShare = min(1.0, float(n2w)) / (1.0 + float(n2w) + float(sig2))
        pairMinShareGeom = np.sqrt(n1MinShare * n2MinShare)

        # Paired weight geom-mean (existing)
        pairWeightGeomMean = shared_geom_share_rectified(w1=n1w, w2=n2w, sigma1=sig1, sigma2=sig2)
        masterPairWeight.append(pairWeightGeomMean)
        masterCorr.append(pairStimCorr)

        # --- single-stim correlations ---
        n1_single_center = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, center_single_idx]
        n2_single_center = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, center_single_idx]
        singleCorr_center, singleCenterCov, _ = pairCorrExclude3SD(n1_single_center, n2_single_center,
                                                                   method='pearson')

        # --- SINGLE (center) avg FR across the two neurons (Hz) ---
        n1_mean_spk_single_c = float(np.mean(n1_single_center))
        n2_mean_spk_single_c = float(np.mean(n2_single_center))
        n1_fr_single_center_hz = n1_mean_spk_single_c * HZ_PER_COUNT
        n2_fr_single_center_hz = n2_mean_spk_single_c * HZ_PER_COUNT
        pairMeanFR_single_center_hz = 0.5 * (n1_fr_single_center_hz + n2_fr_single_center_hz)

        n1_single_periph = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, periph_single_idx]
        n2_single_periph = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, periph_single_idx]
        singleCorr_periph, singlePeriCov, _ = pairCorrExclude3SD(n1_single_periph, n2_single_periph,
                                                                 method='pearson')

        # --- TRANSIENT singles ---
        n1_single_center_T = masterCorrUnitsSpikeCountMatTransient[n1SpkMatIdx][:, center_single_idx]
        n2_single_center_T = masterCorrUnitsSpikeCountMatTransient[n2SpkMatIdx][:, center_single_idx]
        singleCorr_center_T, _, _ = pairCorrExclude3SD(n1_single_center_T, n2_single_center_T, method='pearson')

        n1_single_periph_T = masterCorrUnitsSpikeCountMatTransient[n1SpkMatIdx][:, periph_single_idx]
        n2_single_periph_T = masterCorrUnitsSpikeCountMatTransient[n2SpkMatIdx][:, periph_single_idx]
        singleCorr_periph_T, _, _ = pairCorrExclude3SD(n1_single_periph_T, n2_single_periph_T, method='pearson')
        # --- SUSTAINED singles ---
        n1_single_center_S = masterCorrUnitsSpikeCountMatSustained[n1SpkMatIdx][:, center_single_idx]
        n2_single_center_S = masterCorrUnitsSpikeCountMatSustained[n2SpkMatIdx][:, center_single_idx]
        singleCorr_center_S, _, _ = pairCorrExclude3SD(n1_single_center_S, n2_single_center_S, method='pearson')

        n1_single_periph_S = masterCorrUnitsSpikeCountMatSustained[n1SpkMatIdx][:, periph_single_idx]
        n2_single_periph_S = masterCorrUnitsSpikeCountMatSustained[n2SpkMatIdx][:, periph_single_idx]
        singleCorr_periph_S, _, _ = pairCorrExclude3SD(n1_single_periph_S, n2_single_periph_S, method='pearson')

        # --- SINGLE (periph) avg FR across the two neurons (Hz) ---
        n1_mean_spk_single_p = float(np.mean(n1_single_periph))
        n2_mean_spk_single_p = float(np.mean(n2_single_periph))
        n1_fr_single_periph_hz = n1_mean_spk_single_p * HZ_PER_COUNT
        n2_fr_single_periph_hz = n2_mean_spk_single_p * HZ_PER_COUNT
        pairMeanFR_single_periph_hz = 0.5 * (n1_fr_single_periph_hz + n2_fr_single_periph_hz)

        # --- single-stim weight geom-means ---
        # center: w0=1 -> share = 1/(1+sigma)
        s1_center = 1.0 / (1.0 + sig1)
        s2_center = 1.0 / (1.0 + sig2)
        singleWeightGeom_center = np.sqrt(s1_center * s2_center)

        # periphery: share = w/(w+sigma)
        s1_periph = n1w / (n1w + sig1)
        s2_periph = n2w / (n2w + sig2)
        singleWeightGeom_periph = np.sqrt(s1_periph * s2_periph)

        # Existing avg of each neuron's (nW, w0=1), then across neurons
        avg_n1_2w = 0.5 * (float(n1w) + w0)
        avg_n2_2w = 0.5 * (float(n2w) + w0)
        avg_pair_weight_2w = 0.5 * (avg_n1_2w + avg_n2_2w)

        # ----- Selectivity (center=N) -----
        loc0Resp = rfWeightNeuronFits[n1][1] / (1 + rfWeightNeuronFits[n1][6])
        loc0Resp = loc0Resp + rfWeightNeuronFits[n1][7]
        if i % 2 == 0:
            loc1Resp = (rfWeightNeuronFits[n1][1] * rfWeightNeuronFits[n1][2 + site]) / \
                       (rfWeightNeuronFits[n1][2 + site] + rfWeightNeuronFits[n1][6])
        else:
            loc1Resp = (rfWeightNeuronFits[n1][0] * rfWeightNeuronFits[n1][2 + site]) / \
                       (rfWeightNeuronFits[n1][2 + site] + rfWeightNeuronFits[n1][6])
        loc1Resp = loc1Resp + rfWeightNeuronFits[n1][7]
        n1Sel = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)

        loc0Resp = rfWeightNeuronFits[n2][1] / (1 + rfWeightNeuronFits[n2][6])
        loc0Resp = loc0Resp + rfWeightNeuronFits[n2][7]
        if i % 2 == 0:
            loc1Resp = (rfWeightNeuronFits[n2][1] * rfWeightNeuronFits[n2][2 + site]) / \
                       (rfWeightNeuronFits[n2][2 + site] + rfWeightNeuronFits[n2][6])
        else:
            loc1Resp = (rfWeightNeuronFits[n2][0] * rfWeightNeuronFits[n2][2 + site]) / \
                       (rfWeightNeuronFits[n2][2 + site] + rfWeightNeuronFits[n2][6])
        loc1Resp = loc1Resp + rfWeightNeuronFits[n2][7]
        n2Sel = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)

        masterPairSel.append(np.sign(n1Sel) * np.sign(n2Sel) * np.sqrt(abs(n1Sel) * abs(n2Sel)))

        if np.sign(n1Sel) * np.sign(n2Sel) == 1:
            sameSelCorr.append(pairStimCorr)
            n1NonprefWeight = (n1w / (1 + n1w)) if n1Sel >= 0 else (1 / (1 + n1w))
            n2NonprefWeight = (n2w / (1 + n2w)) if n2Sel >= 0 else (1 / (1 + n2w))
            sameSelWeight.append(np.sqrt(n1NonprefWeight * n2NonprefWeight))
        else:
            oppositeSelCorr.append(pairStimCorr)
            n1NonprefWeight = (n1w / (1 + n1w)) if n1Sel >= 0 else (1 / (1 + n1w))
            n2NonprefWeight = (n2w / (1 + n2w)) if n2Sel >= 0 else (1 / (1 + n2w))
            oppositeSelWeight.append(np.sqrt(n1NonprefWeight * n2NonprefWeight))

        # ---- Gate the four fields by R² ----
        deltaPD_keep = deltaPD if r2_ok else np.nan
        avg_pair_weight_2w_keep = avg_pair_weight_2w if r2_ok else np.nan
        bc_keep = bc_val if r2_ok else np.nan
        bd_keep = bd_val if r2_ok else np.nan

        # ---- add_record call (Block 1 example) ----
        add_record(block="centerN_offsetsPN", i=i, comb=comb,
                   pairStimCorr=pairStimCorr,
                   pairWeightGeomMean=pairWeightGeomMean,
                   n1Sel=n1Sel, n2Sel=n2Sel,
                   sim_SI=sim_val, sim_SI01=sim01_val,
                   delta_PD_abs_deg=deltaPD_keep,  # <— gated
                   avg_pair_weight_2w=avg_pair_weight_2w_keep,  # <— gated
                   monkey=monkey_name,
                   singleCorr_center=singleCorr_center,
                   singleCorr_periph=singleCorr_periph,
                   singleWeightGeom_center=singleWeightGeom_center,
                   singleWeightGeom_periph=singleWeightGeom_periph,
                   BC_dir_tuning=bc_keep,  # <— gated
                   BDir=bd_keep,  # <— gated
                   n1w=n1w,
                   n2w=n2w,
                   pairMinShareGeom=pairMinShareGeom,
                   n1MinShare=n1MinShare,
                   n2MinShare=n2MinShare,
                   n1Sig=sig1,
                   n2Sig=sig2,
                   pairCov=pairDCov[0][1],
                   singleCCov=singleCenterCov[0][1],
                   singlePCov=singlePeriCov[0][1],
                   pairMeanFR_paired_hz=pairMeanFR_paired_hz,
                   pairStimCorr_transient=pairStimCorr_T,
                   pairStimCorr_sustained=pairStimCorr_S,
                   singleCorr_center_transient=singleCorr_center_T,
                   singleCorr_periph_transient=singleCorr_periph_T,
                   singleCorr_center_sustained=singleCorr_center_S,
                   singleCorr_periph_sustained=singleCorr_periph_S,
                   pairMeanFR_single_center_hz=pairMeanFR_single_center_hz,
                   pairMeanFR_single_periph_hz=pairMeanFR_single_periph_hz
                   )

    # =========================
    # Block 2: i = 19..26 (center=P offsets P or N)
    # =========================
    for i in range(19, 27):
        n1SpkMatIdx = corr_index[comb[0]]
        n2SpkMatIdx = corr_index[comb[1]]
        n1SpikeMat = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, i]
        n2SpikeMat = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, i]
        pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat,
                                                             method='pearson')

        # --- TRANSIENT (50–100ms) paired correlation ---
        n1SpikeMat_T = masterCorrUnitsSpikeCountMatTransient[n1SpkMatIdx][:, i]
        n2SpikeMat_T = masterCorrUnitsSpikeCountMatTransient[n2SpkMatIdx][:, i]
        pairStimCorr_T, _, _ = pairCorrExclude3SD(n1SpikeMat_T, n2SpikeMat_T, method='pearson')

        # --- SUSTAINED (100ms → offLatency) paired correlation ---
        n1SpikeMat_S = masterCorrUnitsSpikeCountMatSustained[n1SpkMatIdx][:, i]
        n2SpikeMat_S = masterCorrUnitsSpikeCountMatSustained[n2SpkMatIdx][:, i]
        pairStimCorr_S, _, _ = pairCorrExclude3SD(n1SpikeMat_S, n2SpikeMat_S, method='pearson')

        # --- PAIRED avg FR across the two neurons (Hz) ---
        # Each is mean spike count over trials -> convert to Hz -> average across neurons
        n1_mean_spk_paired = float(np.mean(n1SpikeMat))
        n2_mean_spk_paired = float(np.mean(n2SpikeMat))
        n1_fr_paired_hz = n1_mean_spk_paired * HZ_PER_COUNT
        n2_fr_paired_hz = n2_mean_spk_paired * HZ_PER_COUNT
        pairMeanFR_paired_hz = 0.5 * (n1_fr_paired_hz + n2_fr_paired_hz)

        # Correct single indices and site
        center_single_idx, periph_single_idx, site = single_indices_for_pair(i)

        w0 = 1.0
        n1w = rfWeightNeuronFits[n1][2 + site]
        sig1 = rfWeightNeuronFits[n1][6]
        n2w = rfWeightNeuronFits[n2][2 + site]
        sig2 = rfWeightNeuronFits[n2][6]
        if n1w <= 0 or n2w <= 0:
            continue

        # --- NEW: min-share per neuron and geometric mean across pair ---
        # For each neuron: share = min(w0, w1) / (w0 + w1 + sigma); here w0=1, w1 = that neuron's site weight
        n1MinShare = min(1.0, float(n1w)) / (1.0 + float(n1w) + float(sig1))
        n2MinShare = min(1.0, float(n2w)) / (1.0 + float(n2w) + float(sig2))
        pairMinShareGeom = np.sqrt(n1MinShare * n2MinShare)

        pairWeightGeomMean = shared_geom_share_rectified(w1=n1w, w2=n2w, sigma1=sig1, sigma2=sig2)
        masterPairWeight.append(pairWeightGeomMean)
        masterCorr.append(pairStimCorr)

        # --- single-stim correlations ---
        n1_single_center = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, center_single_idx]
        n2_single_center = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, center_single_idx]
        singleCorr_center, singleCenterCov, _ = pairCorrExclude3SD(n1_single_center, n2_single_center,
                                                                   method='pearson')

        # --- SINGLE (center) avg FR across the two neurons (Hz) ---
        n1_mean_spk_single_c = float(np.mean(n1_single_center))
        n2_mean_spk_single_c = float(np.mean(n2_single_center))
        n1_fr_single_center_hz = n1_mean_spk_single_c * HZ_PER_COUNT
        n2_fr_single_center_hz = n2_mean_spk_single_c * HZ_PER_COUNT
        pairMeanFR_single_center_hz = 0.5 * (n1_fr_single_center_hz + n2_fr_single_center_hz)

        n1_single_periph = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, periph_single_idx]
        n2_single_periph = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, periph_single_idx]
        singleCorr_periph, singlePeriCov, _ = pairCorrExclude3SD(n1_single_periph, n2_single_periph,
                                                                 method='pearson')

        # --- TRANSIENT singles ---
        n1_single_center_T = masterCorrUnitsSpikeCountMatTransient[n1SpkMatIdx][:, center_single_idx]
        n2_single_center_T = masterCorrUnitsSpikeCountMatTransient[n2SpkMatIdx][:, center_single_idx]
        singleCorr_center_T, _, _ = pairCorrExclude3SD(n1_single_center_T, n2_single_center_T, method='pearson')

        n1_single_periph_T = masterCorrUnitsSpikeCountMatTransient[n1SpkMatIdx][:, periph_single_idx]
        n2_single_periph_T = masterCorrUnitsSpikeCountMatTransient[n2SpkMatIdx][:, periph_single_idx]
        singleCorr_periph_T, _, _ = pairCorrExclude3SD(n1_single_periph_T, n2_single_periph_T, method='pearson')

        # --- SUSTAINED singles ---
        n1_single_center_S = masterCorrUnitsSpikeCountMatSustained[n1SpkMatIdx][:, center_single_idx]
        n2_single_center_S = masterCorrUnitsSpikeCountMatSustained[n2SpkMatIdx][:, center_single_idx]
        singleCorr_center_S, _, _ = pairCorrExclude3SD(n1_single_center_S, n2_single_center_S, method='pearson')

        n1_single_periph_S = masterCorrUnitsSpikeCountMatSustained[n1SpkMatIdx][:, periph_single_idx]
        n2_single_periph_S = masterCorrUnitsSpikeCountMatSustained[n2SpkMatIdx][:, periph_single_idx]
        singleCorr_periph_S, _, _ = pairCorrExclude3SD(n1_single_periph_S, n2_single_periph_S, method='pearson')

        # --- SINGLE (periph) avg FR across the two neurons (Hz) ---
        n1_mean_spk_single_p = float(np.mean(n1_single_periph))
        n2_mean_spk_single_p = float(np.mean(n2_single_periph))
        n1_fr_single_periph_hz = n1_mean_spk_single_p * HZ_PER_COUNT
        n2_fr_single_periph_hz = n2_mean_spk_single_p * HZ_PER_COUNT
        pairMeanFR_single_periph_hz = 0.5 * (n1_fr_single_periph_hz + n2_fr_single_periph_hz)

        # --- single-stim weight geom-means ---
        s1_center = 1.0 / (1.0 + sig1)
        s2_center = 1.0 / (1.0 + sig2)
        singleWeightGeom_center = np.sqrt(s1_center * s2_center)

        s1_periph = n1w / (n1w + sig1)
        s2_periph = n2w / (n2w + sig2)
        singleWeightGeom_periph = np.sqrt(s1_periph * s2_periph)

        # Existing avg of each neuron's (nW, w0=1)
        avg_n1_2w = 0.5 * (float(n1w) + w0)
        avg_n2_2w = 0.5 * (float(n2w) + w0)
        avg_pair_weight_2w = 0.5 * (avg_n1_2w + avg_n2_2w)

        # ----- Selectivity (center=P) -----
        loc0Resp = rfWeightNeuronFits[n1][0] / (1 + rfWeightNeuronFits[n1][6])
        loc0Resp = loc0Resp + rfWeightNeuronFits[n1][7]
        if i % 2 == 0:
            loc1Resp = (rfWeightNeuronFits[n1][0] * rfWeightNeuronFits[n1][2 + site]) / \
                       (rfWeightNeuronFits[n1][2 + site] + rfWeightNeuronFits[n1][6])
        else:
            loc1Resp = (rfWeightNeuronFits[n1][1] * rfWeightNeuronFits[n1][2 + site]) / \
                       (rfWeightNeuronFits[n1][2 + site] + rfWeightNeuronFits[n1][6])
        loc1Resp = loc1Resp + rfWeightNeuronFits[n1][7]
        n1Sel = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)

        loc0Resp = rfWeightNeuronFits[n2][0] / (1 + rfWeightNeuronFits[n2][6])
        loc0Resp = loc0Resp + rfWeightNeuronFits[n2][7]
        if i % 2 == 0:
            loc1Resp = (rfWeightNeuronFits[n2][0] * rfWeightNeuronFits[n2][2 + site]) / \
                       (rfWeightNeuronFits[n2][2 + site] + rfWeightNeuronFits[n2][6])
        else:
            loc1Resp = (rfWeightNeuronFits[n2][1] * rfWeightNeuronFits[n2][2 + site]) / \
                       (rfWeightNeuronFits[n2][2 + site] + rfWeightNeuronFits[n2][6])
        loc1Resp = loc1Resp + rfWeightNeuronFits[n2][7]
        n2Sel = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)

        masterPairSel.append(np.sign(n1Sel) * np.sign(n2Sel) * np.sqrt(abs(n1Sel) * abs(n2Sel)))

        if np.sign(n1Sel) * np.sign(n2Sel) == 1:
            sameSelCorr.append(pairStimCorr)
            n1NonprefWeight = (n1w / (1 + n1w)) if n1Sel >= 0 else (1 / (1 + n1w))
            n2NonprefWeight = (n2w / (1 + n2w)) if n2Sel >= 0 else (1 / (1 + n2w))
            sameSelWeight.append(np.sqrt(n1NonprefWeight * n2NonprefWeight))
        else:
            oppositeSelCorr.append(pairStimCorr)
            n1NonprefWeight = (n1w / (1 + n1w)) if n1Sel >= 0 else (1 / (1 + n1w))
            n2NonprefWeight = (n2w / (1 + n2w)) if n2Sel >= 0 else (1 / (1 + n2w))
            oppositeSelWeight.append(np.sqrt(n1NonprefWeight * n2NonprefWeight))

        # ---- Gate the four fields by R² ----
        deltaPD_keep = deltaPD if r2_ok else np.nan
        avg_pair_weight_2w_keep = avg_pair_weight_2w if r2_ok else np.nan
        bc_keep = bc_val if r2_ok else np.nan
        bd_keep = bd_val if r2_ok else np.nan

        # ---- add_record call (Block 1 example) ----
        add_record(block="centerP_offsetsPN", i=i, comb=comb,
                   pairStimCorr=pairStimCorr,
                   pairWeightGeomMean=pairWeightGeomMean,
                   n1Sel=n1Sel, n2Sel=n2Sel,
                   sim_SI=sim_val, sim_SI01=sim01_val,
                   delta_PD_abs_deg=deltaPD_keep,  # <— gated
                   avg_pair_weight_2w=avg_pair_weight_2w_keep,  # <— gated
                   monkey=monkey_name,
                   singleCorr_center=singleCorr_center,
                   singleCorr_periph=singleCorr_periph,
                   singleWeightGeom_center=singleWeightGeom_center,
                   singleWeightGeom_periph=singleWeightGeom_periph,
                   BC_dir_tuning=bc_keep,  # <— gated
                   BDir=bd_keep,  # <— gated
                   n1w=n1w,
                   n2w=n2w,
                   pairMinShareGeom=pairMinShareGeom,
                   n1MinShare=n1MinShare,
                   n2MinShare=n2MinShare,
                   n1Sig=sig1,
                   n2Sig=sig2,
                   pairCov=pairDCov[0][1],
                   singleCCov=singleCenterCov[0][1],
                   singlePCov=singlePeriCov[0][1],
                   pairStimCorr_transient=pairStimCorr_T,
                   pairStimCorr_sustained=pairStimCorr_S,
                   singleCorr_center_transient=singleCorr_center_T,
                   singleCorr_periph_transient=singleCorr_periph_T,
                   singleCorr_center_sustained=singleCorr_center_S,
                   singleCorr_periph_sustained=singleCorr_periph_S,
                   pairMeanFR_paired_hz=pairMeanFR_paired_hz,
                   pairMeanFR_single_center_hz=pairMeanFR_single_center_hz,
                   pairMeanFR_single_periph_hz=pairMeanFR_single_periph_hz
                   )

# Wrap up: arrays if you still need them
masterCorr = np.array(masterCorr)
masterPairWeight = np.array(masterPairWeight)
masterPairSel = np.array(masterPairSel)
sameSelCorr = np.array(sameSelCorr)
sameSelWeight = np.array(sameSelWeight)
oppositeSelCorr = np.array(oppositeSelCorr)
oppositeSelWeight = np.array(oppositeSelWeight)

# DataFrame for analysis/joins/plots
corr_df = pd.DataFrame(corr_records)

savemat("corr_dataMTND.mat", {
    "nonPrefSupp": corr_df["pair_weight_geommean"].values,
    "selectivity": corr_df["pairSel_geom"].values,
    "bhattCoef": corr_df["BC_dir_tuning"].values,
    "deltaPD": corr_df["delta_PD_abs_deg"].values,
    "pairedCorr": corr_df['corr'].values,
    "pairedCorr_transient": corr_df["corr_transient"].values,
    "pairedCorr_sustained": corr_df["corr_sustained"].values,
    "singleCorrLoc0": corr_df['singleCorr_center'].values,
    "singleCorrLoc1": corr_df['singleCorr_periph'].values,
    "singleCorrLoc0_transient": corr_df["singleCorr_center_transient"].values,
    "singleCorrLoc1_transient": corr_df["singleCorr_periph_transient"].values,
    "singleCorrLoc0_sustained": corr_df["singleCorr_center_sustained"].values,
    "singleCorrLoc1_sustained": corr_df["singleCorr_periph_sustained"].values,
    "singleCenterNonPrefSupp": corr_df['singleWeightGeom_center'].values,
    "singlePeriphNonPrefSupp": corr_df['singleWeightGeom_periph'].values,
    "monkey": corr_df['monkey'].values.astype('U'),  # save as strings
    "pairMinShareGeom": corr_df["pairMinShareGeom"].values,
    "pairCov": corr_df["pairCov"].values,
    "singleCenterCov": corr_df["singleCenterCov"].values,
    "singlePeriphCov": corr_df["singlePeriphCov"].values,
    "pairMeanFR_paired_hz": corr_df["pairMeanFR_paired_hz"].values,
    "pairMeanFR_single_center_hz": corr_df["pairMeanFR_single_center_hz"].values,
    "pairMeanFR_single_periph_hz": corr_df["pairMeanFR_single_periph_hz"].values,
})

# paired
x = corr_df['pair_weight_geommean'].values
y = corr_df['BC_dir_tuning'].values
c = corr_df['corr'].values

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

# single
c = np.concatenate((corr_df['singleCorr_center'].values,
                    corr_df['singleCorr_periph'].values), axis=0)
x = np.concatenate((corr_df['singleWeightGeom_center'].values,
                    corr_df['singleWeightGeom_periph'].values), axis=0)
y = np.concatenate((corr_df['delta_PD_abs_deg'].values,
                    corr_df['delta_PD_abs_deg'].values), axis=0)

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

# plot correlations vs Delta PD
for filler in range(1):
    # Extract arrays
    deltaPD = corr_df["delta_PD_abs_deg"].to_numpy(dtype=float)
    corr_paired = corr_df["corr"].to_numpy(dtype=float)
    corr_center = corr_df["singleCorr_center"].to_numpy(dtype=float)
    corr_periph = corr_df["singleCorr_periph"].to_numpy(dtype=float)

    # Stack (3 * N)
    x_all = np.concatenate([deltaPD, deltaPD, deltaPD])
    y_all = np.concatenate([corr_paired, corr_center, corr_periph])

    # Keep only finite values
    msk = np.isfinite(x_all) & np.isfinite(y_all)
    x_all = x_all[msk]
    y_all = y_all[msk]

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # ========= Left panel: fixed-width bins =========
    axs[0].scatter(x_all, y_all, alpha=0.3, s=10, color="gray", label="data")

    nbins = 10
    bins = np.linspace(0, 180, nbins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_means = np.full(nbins, np.nan)

    for j in range(nbins):
        m = (x_all >= bins[j]) & (x_all < bins[j + 1])
        if np.any(m):
            bin_means[j] = np.nanmean(y_all[m])

    valid = np.isfinite(bin_means)
    axs[0].plot(bin_centers[valid], bin_means[valid], "-ok", lw=2, label="bin mean")

    axs[0].set_xlim(0, 180)
    axs[0].set_xlabel("ΔPD (deg)")
    axs[0].set_ylabel("Spike count correlation")
    axs[0].set_title("Fixed-width bins")
    axs[0].legend()

    # ========= Right panel: equal-count bins =========
    axs[1].scatter(x_all, y_all, alpha=0.3, s=10, color="gray", label="data")

    nbins_eq = 10
    # Quantile-based bin edges
    quantiles = np.linspace(0, 1, nbins_eq + 1)
    bin_edges = np.quantile(x_all, quantiles)
    bin_centers_eq = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_means_eq = np.full(nbins_eq, np.nan)

    for j in range(nbins_eq):
        m = (x_all >= bin_edges[j]) & (x_all < bin_edges[j + 1])
        if np.any(m):
            bin_means_eq[j] = np.nanmean(y_all[m])

    valid_eq = np.isfinite(bin_means_eq)
    axs[1].plot(bin_centers_eq[valid_eq], bin_means_eq[valid_eq], "-ok", lw=2, label="bin mean")

    axs[1].set_xlim(0, 180)
    axs[1].set_xlabel("ΔPD (deg)")
    axs[1].set_title("Equal-count bins")
    axs[1].legend()

    plt.tight_layout()
    plt.show()

# plot correlations vs BC
for filler in range(1):
    # Extract arrays
    bc = corr_df["BC_dir_tuning"].to_numpy(dtype=float)
    corr_paired = corr_df["corr"].to_numpy(dtype=float)
    corr_center = corr_df["singleCorr_center"].to_numpy(dtype=float)
    corr_periph = corr_df["singleCorr_periph"].to_numpy(dtype=float)

    # Stack (3 * N)
    x_all = np.concatenate([bc, bc, bc])
    y_all = np.concatenate([corr_paired, corr_center, corr_periph])

    # Keep only finite points and BC within [0,1]
    msk = np.isfinite(x_all) & np.isfinite(y_all) & (x_all >= 0) & (x_all <= 1)
    x_all = x_all[msk]
    y_all = y_all[msk]

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # ========= Left panel: fixed-width bins =========
    axs[0].scatter(x_all, y_all, alpha=0.3, s=10, color="gray", label="data")

    nbins = 12
    bins = np.linspace(0.0, 1.0, nbins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_means = np.full(nbins, np.nan)

    for j in range(nbins):
        m = (x_all >= bins[j]) & (x_all < bins[j + 1])
        if np.any(m):
            bin_means[j] = np.nanmean(y_all[m])

    valid = np.isfinite(bin_means)
    axs[0].plot(bin_centers[valid], bin_means[valid], "-ok", lw=2, label="bin mean")

    axs[0].set_xlim(0, 1)
    axs[0].set_xlabel("Bhattacharyya Coefficient (0–1)")
    axs[0].set_ylabel("Spike count correlation")
    axs[0].set_title("Fixed-width bins")
    axs[0].legend()

    # ========= Right panel: equal-count bins =========
    axs[1].scatter(x_all, y_all, alpha=0.3, s=10, color="gray", label="data")

    nbins_eq = 12
    # Quantile bin edges
    quantiles = np.linspace(0, 1, nbins_eq + 1)
    bin_edges = np.quantile(x_all, quantiles)
    bin_centers_eq = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_means_eq = np.full(nbins_eq, np.nan)

    for j in range(nbins_eq):
        m = (x_all >= bin_edges[j]) & (x_all < bin_edges[j + 1])
        if np.any(m):
            bin_means_eq[j] = np.nanmean(y_all[m])

    valid_eq = np.isfinite(bin_means_eq)
    axs[1].plot(bin_centers_eq[valid_eq], bin_means_eq[valid_eq], "-ok", lw=2, label="bin mean")

    axs[1].set_xlim(0, 1)
    axs[1].set_xlabel("Bhattacharyya Coefficient (0–1)")
    axs[1].set_title("Equal-count bins")
    axs[1].legend()

    plt.tight_layout()
    plt.show()

# binned plot of correlations vs pair weight
for filler in range(1):
    # 1) filter to 0–2 on x
    mask = (masterPairWeight >= 0) & (masterPairWeight <= 2)
    x = masterPairWeight[mask]
    y = masterCorr[mask]


    # 2) choose bins
    bin_width = 0.1
    bins = np.arange(0, 2 + bin_width, bin_width)  # inclusive of 2
    bin_centers = (bins[:-1] + bins[1:]) / 2
    n_bins = len(bins) - 1

    # 3) compute bin means and SEM
    bin_idx = np.digitize(x, bins) - 1
    valid_points = (bin_idx >= 0) & (bin_idx < n_bins)
    bin_idx = bin_idx[valid_points]
    xv = x[valid_points]
    yv = y[valid_points]

    counts = np.bincount(bin_idx, minlength=n_bins)
    y_sums = np.bincount(bin_idx, weights=yv, minlength=n_bins)
    y2_sums = np.bincount(bin_idx, weights=yv ** 2, minlength=n_bins)

    bin_means = np.full(n_bins, np.nan)
    bin_sems = np.full(n_bins, np.nan)

    nonempty = counts > 0
    bin_means[nonempty] = y_sums[nonempty] / counts[nonempty]

    # unbiased variance (ddof=1) for bins with n >= 2
    ge2 = counts >= 2
    bin_vars = np.full(n_bins, np.nan)
    bin_vars[ge2] = (y2_sums[ge2] - counts[ge2] * (bin_means[ge2] ** 2)) / (counts[ge2] - 1)
    bin_sems[ge2] = np.sqrt(bin_vars[ge2] / counts[ge2])

    # --- compute per-bin medians ---
    bin_medians = np.full(n_bins, np.nan)
    for k in range(n_bins):
        if counts[k] > 0:
            bin_medians[k] = np.median(yv[bin_idx == k])

    # 4) plot raw + overlay means with SEM error bars + median line
    plt.scatter(x, y, s=12, alpha=0.6, label="raw")

    plt.errorbar(
        bin_centers[nonempty],
        bin_means[nonempty],
        yerr=bin_sems[nonempty],
        fmt='o-',
        linewidth=2,
        color='black',
        capsize=3,
        label=f"bin means ± SEM (bin={bin_width})",
    )

    # median line in magenta
    plt.plot(
        bin_centers[nonempty],
        bin_medians[nonempty],
        marker='o',
        linewidth=2,
        color='red',
        label="bin median",
    )

    plt.xlabel("Pair Weight (geometric mean)")
    plt.ylabel("Spike Count Correlations")
    plt.title("Raw scatter with binned means ± SEM and median")
    plt.legend()
    plt.show()

# equally populated bins
for filler in range(1):
    # 1) filter to 0–2 on x
    mask = (masterPairWeight >= 0) & (masterPairWeight <= 2)
    x = masterPairWeight[mask]
    y = masterCorr[mask]

    # x = sameSelWeight
    # y = sameSelCorr

    # x = oppositeSelWeight
    # y = oppositeSelCorr

    # 2) choose number of equal-count bins
    n_bins_target = 10  # tweak as you like

    # quantile-based bin edges (equal-count bins)
    edges = np.quantile(x, np.linspace(0, 1, n_bins_target + 1))

    # handle duplicates (ties) in edges → drop duplicates, reduce bin count
    edges = np.unique(edges)
    n_bins = max(1, len(edges) - 1)

    # if everything collapsed to one edge (e.g., constant x), just plot scatter
    if n_bins < 1:
        plt.scatter(x, y, s=12, alpha=0.6, label="raw")
        plt.xlabel("masterPairWeight (0–2)")
        plt.ylabel("masterCorr")
        plt.title("Raw scatter (binning unavailable due to constant x)")
        plt.legend()
        plt.show()
    else:
        # bin centers (midpoints in x; widths vary with quantiles)
        bin_centers = (edges[:-1] + edges[1:]) / 2

        # 3) compute bin means and SEM
        # use searchsorted like digitize with right-closed bins
        bin_idx = np.searchsorted(edges, x, side='right') - 1
        valid_points = (bin_idx >= 0) & (bin_idx < n_bins)
        bin_idx = bin_idx[valid_points]
        xv = x[valid_points]
        yv = y[valid_points]

        counts = np.bincount(bin_idx, minlength=n_bins)
        y_sums = np.bincount(bin_idx, weights=yv, minlength=n_bins)
        y2_sums = np.bincount(bin_idx, weights=yv ** 2, minlength=n_bins)

        bin_means = np.full(n_bins, np.nan)
        bin_sems = np.full(n_bins, np.nan)

        nonempty = counts > 0
        bin_means[nonempty] = y_sums[nonempty] / counts[nonempty]

        # unbiased variance (ddof=1) for bins with n >= 2
        ge2 = counts >= 2
        bin_vars = np.full(n_bins, np.nan)
        bin_vars[ge2] = (y2_sums[ge2] - counts[ge2] * (bin_means[ge2] ** 2)) / (counts[ge2] - 1)
        bin_sems[ge2] = np.sqrt(bin_vars[ge2] / counts[ge2])

        # 4) plot raw + overlay means with SEM error bars
        plt.scatter(x, y, s=12, alpha=0.6, label="raw")
        plt.errorbar(
            bin_centers[nonempty],
            bin_means[nonempty],
            yerr=bin_sems[nonempty],
            fmt='o-',
            linewidth=2,
            capsize=3,
            color='black',
            label=f"equal-count bin means ± SEM (bins≈{n_bins_target})"
        )
        plt.xlabel("masterPairWeight (0–2)")
        plt.ylabel("masterCorr")
        plt.title("Raw scatter with equal-count binned means ± SEM")
        plt.legend()
        plt.show()

# 3d plot
fig = plt.figure()
cm = plt.colormaps.get_cmap('seismic')
x = masterPairWeight  # non-pref suppression
y = masterPairSel  # neuron pair selectivity
z = masterCorr  # noise correlations for paired Gabors
sc = plt.scatter(x, y, c=z, s=15, cmap=cm, vmin=-1, vmax=1)

# plot characteristics
plt.colorbar(sc)
plt.xlabel('non-pref suppression')
plt.ylabel('selectivity')
plt.show()

########## PLOTTING ###########
for filler in range(1):
    # collapse to average corr per pair
    pair_means = corr_df.groupby(['unit_i', 'unit_j']).agg(
        mean_corr=('corr', 'mean'),
        mean_weight=('pair_weight_geommean', 'mean'),
        SI=('similarity_SI', 'first'),   # same for all rows of that pair
        SI01=('similarity_SI01', 'first')
    ).reset_index()

    # 1) Define user bins (edit as you like)
    # Example: 8 bins uniformly from -1 to 1
    bins = np.linspace(-1.0, 1.0, 9)   # <- change this to your preferred bin edges

    # 2) Compute bin index for each pair (ignore NaNs in SI)
    valid = pair_means['SI'].notna() & pair_means['mean_corr'].notna()
    si_vals = pair_means.loc[valid, 'SI'].to_numpy()
    mc_vals = pair_means.loc[valid, 'mean_corr'].to_numpy()

    # digitize: returns 1..len(bins) - we shift to 0-based bin index
    bin_idx = np.digitize(si_vals, bins, right=False) - 1

    # 3) Aggregate mean and SEM per bin
    centers, means, sems, counts = [], [], [], []
    for b in range(len(bins) - 1):
        m = bin_idx == b
        if not np.any(m):
            continue
        y = mc_vals[m]
        centers.append(0.5 * (bins[b] + bins[b+1]))
        means.append(np.mean(y))
        # SEM with ddof=1 if count>1; fallback to 0 for singletons
        if y.size > 1:
            sems.append(np.std(y, ddof=1) / np.sqrt(y.size))
        else:
            sems.append(0.0)
        counts.append(y.size)

    centers = np.array(centers)
    means = np.array(means)
    sems = np.array(sems)
    counts = np.array(counts)

    # 4) Overlay line + error bars on top of your scatter
    plt.figure(figsize=(6,5))
    plt.scatter(pair_means['SI'], pair_means['mean_corr'], alpha=0.7)
    plt.errorbar(centers, means, yerr=sems, fmt='o-', capsize=4, color='black')
    plt.xlabel("Similarity (cos Δθ, -1 to 1)")
    plt.ylabel("Mean correlation across conditions")
    plt.title("Pairwise correlation vs. tuning similarity (with binned means)")
    plt.axhline(0, color='k', linestyle='--', alpha=0.5)
    plt.axvline(0, color='k', linestyle='--', alpha=0.5)
    plt.show()

    # # (optional) print how many pairs landed in each bin
    # for c, left, right in zip(counts, centers - (bins[1]-bins[0])/2, centers + (bins[1]-bins[0])/2):
    #     print(f"Bin [{left:.2f}, {right:.2f}): n={c}")


# ===============================================
# Scatter + binned means, concatenating single & paired
# - Concatenates: corr, singleCorr_center, singleCorr_periph
# - 2×2 subplots with requested splits and bin-mean overlays
# - Colors/markers show condition_type (paired / single_center / single_periph)
# Paste & run (expects `corr_df` in scope)
# ===============================================

# ==== Helpers ====
def _finite_xy(df, xcol, ycol):
    """Return finite (x, y) arrays with NaNs/Infs removed."""
    x = df[xcol].to_numpy()
    y = df[ycol].to_numpy()
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def binned_means(x, y, n_bins=N_BINS, use_quantiles=False):
    """
    Compute bin centers and y-means per x-bin.
    - If use_quantiles=True, bins are equal-count (robust to outliers).
    - Otherwise, uniform bins from 1st to 99th percentile (robust range).
    Returns: (bin_centers, y_means, counts)
    """
    x = np.asarray(x)
    y = np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) == 0:
        return np.array([]), np.array([]), np.array([])

    if use_quantiles:
        qs = np.linspace(0, 1, n_bins + 1)
        edges = np.unique(np.quantile(x, qs))
        if edges.size < 3:
            lo, hi = np.min(x), np.max(x)
            if lo == hi:
                return np.array([lo]), np.array([np.mean(y)]), np.array([len(y)])
            edges = np.linspace(lo, hi, n_bins + 1)
    else:
        lo, hi = np.nanpercentile(x, 1), np.nanpercentile(x, 99)
        if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
            lo, hi = np.min(x), np.max(x)
        edges = np.linspace(lo, hi, n_bins + 1)

    inds = np.digitize(x, edges, right=False)
    centers, means, counts = [], [], []
    for b in range(1, len(edges)):
        mask = inds == b
        if np.any(mask):
            xbin = x[mask]
            ybin = y[mask]
            centers.append((edges[b-1] + edges[b]) / 2.0)
            means.append(np.mean(ybin))
            counts.append(np.sum(mask))
    return np.asarray(centers), np.asarray(means), np.asarray(counts)


def scatter_with_binned_means(ax, df_sub, xcol, ycol, xlabel, ylabel, title):
    """
    Scatter all points in df_sub (colored by condition_type) and overlay binned mean line across ALL points.
    """
    # Color/marker per condition (optional, helps visual separation)
    style_map = {
        'paired':        dict(color=None, marker='o'),
        'single_center': dict(color=None, marker='s'),
        'single_periph': dict(color=None, marker='^'),
    }
    # Plot each condition subset
    for cond, style in style_map.items():
        d = df_sub[df_sub['condition_type'] == cond]
        if len(d):
            x, y = _finite_xy(d, xcol, ycol)
            if len(x):
                ax.scatter(x, y, s=MARKER_SIZE, alpha=ALPHA_SCATTER,
                           marker=style['marker'], label=cond)

    # Binned means over ALL points combined
    x_all, y_all = _finite_xy(df_sub, xcol, ycol)
    cx, my, cnt = binned_means(x_all, y_all, n_bins=N_BINS, use_quantiles=False)
    if len(cx) > 0:
        ax.plot(cx, my, linewidth=LINE_WIDTH, color='black', zorder=5)
        ax.scatter(cx, my, color='black', s=MARKER_SIZE+4, zorder=6)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    ax.legend(frameon=False, fontsize=9, loc='best')


# ==== CONFIG ====
N_BINS = 10              # number of bins for the line-of-means
ALPHA_SCATTER = 0.10     # transparency for scatter points
MARKER_SIZE = 12         # scatter marker size
LINE_WIDTH = 2.0         # thickness for binned-mean line
bc_thresh = 0.50         # threshold for BC_dir_tuning split
pw_thresh = 0.50         # threshold for pair_weight_geommean split

# ==== Prepare concatenated data (paired + singles) ====
df = corr_df.copy()
df = df.replace([np.inf, -np.inf], np.nan)

for filler in range(1):
    # Make sure needed columns exist
    needed = [
        'corr', 'singleCorr_center', 'singleCorr_periph',
        'pair_weight_geommean', 'BC_dir_tuning'
    ]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Long reshape to concatenate all three correlations
    df_long = pd.melt(
        df,
        id_vars=[c for c in df.columns if c not in ['corr', 'singleCorr_center', 'singleCorr_periph']],
        value_vars=['corr', 'singleCorr_center', 'singleCorr_periph'],
        var_name='condition_type',
        value_name='corr_value'
    )
    df_long['condition_type'] = df_long['condition_type'].map({
        'corr': 'paired',
        'singleCorr_center': 'single_center',
        'singleCorr_periph': 'single_periph'
    })
    df_long = df_long[np.isfinite(df_long['corr_value'])].copy()

    # Subsets for the requested splits (based on original x variables)
    mask_bc_lo = df_long['BC_dir_tuning'] < bc_thresh
    mask_bc_hi = df_long['BC_dir_tuning'] >= bc_thresh

    mask_pw_lo = df_long['pair_weight_geommean'] < pw_thresh
    mask_pw_hi = df_long['pair_weight_geommean'] >= pw_thresh

    # ==== Plot ====
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    # Row 1, Col 1: corr_value vs pair_weight_geommean for BC_dir_tuning < bc_thresh
    scatter_with_binned_means(
        axes[0, 0],
        df_long[mask_bc_lo],
        xcol='pair_weight_geommean',
        ycol='corr_value',
        xlabel='pair_weight_geommean',
        ylabel='correlation',
        title=f'All conditions vs pair_weight_geommean | BC_dir_tuning < {bc_thresh}'
    )

    # Row 1, Col 2: corr_value vs pair_weight_geommean for BC_dir_tuning >= bc_thresh
    scatter_with_binned_means(
        axes[0, 1],
        df_long[mask_bc_hi],
        xcol='pair_weight_geommean',
        ycol='corr_value',
        xlabel='pair_weight_geommean',
        ylabel='correlation',
        title=f'All conditions vs pair_weight_geommean | BC_dir_tuning ≥ {bc_thresh}'
    )

    # Row 2, Col 1: corr_value vs BC_dir_tuning for pair_weight_geommean < pw_thresh
    scatter_with_binned_means(
        axes[1, 0],
        df_long[mask_pw_lo],
        xcol='BC_dir_tuning',
        ycol='corr_value',
        xlabel='BC_dir_tuning',
        ylabel='correlation',
        title=f'All conditions vs BC_dir_tuning | pair_weight_geommean < {pw_thresh}'
    )

    # Row 2, Col 2: corr_value vs BC_dir_tuning for pair_weight_geommean ≥ pw_thresh
    scatter_with_binned_means(
        axes[1, 1],
        df_long[mask_pw_hi],
        xcol='BC_dir_tuning',
        ycol='corr_value',
        xlabel='BC_dir_tuning',
        ylabel='correlation',
        title=f'All conditions vs BC_dir_tuning | pair_weight_geommean ≥ {pw_thresh}'
    )

    plt.show()

########################################################################################################################
################################################# Evaluate BC vs deltaPD ###############################################
########################################################################################################################

# ============================== IMPORTS ==============================
import os
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import i0, i1


# ========================= USER-ADJUSTABLE PARAMS ====================
LOW_BC_THRESH = 0.60
HIGH_BC_THRESH = 0.85
R2_MIN = 0.0          # raise to 0.7/0.9 for stricter quality gating
N_DIRS_MIN = 12       # sanity check: minimum # of directions in a tuning curve


# ===================== Von Mises model & helpers =====================
def von_mises_func(theta, baseline, amplitude, pref_dir, kappa):
    """
    theta, pref_dir in radians
    normalized von Mises bump + baseline
    """
    return baseline + amplitude * np.exp(kappa * np.cos(theta - pref_dir)) / (2 * np.pi * i0(kappa))


def r2_score_safe(y, yhat):
    """Robust R^2 that handles NaNs and near-constant vectors."""
    y = np.asarray(y, float); yhat = np.asarray(yhat, float)
    m = np.isfinite(y) & np.isfinite(yhat)
    if m.sum() < 2:
        return np.nan
    y_ = y[m]; yh = yhat[m]
    ss_res = np.sum((y_ - yh)**2)
    ss_tot = np.sum((y_ - np.mean(y_))**2)
    if ss_tot <= 1e-12:
        return 1.0 if ss_res <= 1e-12 else 0.0
    return 1.0 - ss_res/ss_tot


def fit_von_mises_single(responses, theta_rad):
    """
    Fit baseline + amplitude * vonMises.
    Returns: baseline, amplitude, pref_dir(rad in [0,2π)), kappa
    """
    responses = np.asarray(responses, dtype=float)
    theta_rad = np.asarray(theta_rad, dtype=float)
    if np.sum(np.isfinite(responses)) < 4:
        return np.nan, np.nan, np.nan, np.nan

    baseline_guess = np.nanmin(responses)
    amp_guess = max(np.nanmax(responses) - baseline_guess, 0.0)
    pref_guess = theta_rad[np.nanargmax(responses)]
    kappa_guess = 2.0
    p0 = [baseline_guess, amp_guess, pref_guess, kappa_guess]
    bounds = ([-np.inf, 0.0, 0.0, 1e-6],
              [ np.inf,  np.inf, 2*np.pi, 1000.0])
    try:
        params, _ = curve_fit(
            von_mises_func, theta_rad, responses, p0=p0,
            bounds=bounds, maxfev=20000
        )
        b, A, mu, k = params
        mu = np.mod(mu, 2*np.pi)
        return float(b), float(A), float(mu), float(k)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


# ========================== Feature helpers ==========================
def vm_fwhm_from_kappa(kappa):
    """Exact FWHM (radians) of a von Mises bump."""
    kappa = float(kappa)
    x = np.clip(1.0 - np.log(2.0)/max(kappa, 1e-9), -1.0, 1.0)
    return 2.0 * np.arccos(x)


def pmf_from_curve(y, eps=1e-12):
    y = np.clip(np.asarray(y, float), 0, None)
    s = np.sum(y) + eps
    return (y + eps) / s


def bhattacharyya(p, q):
    p = np.clip(np.asarray(p, float), 0, None)
    q = np.clip(np.asarray(q, float), 0, None)
    return float(np.sum(np.sqrt(p * q)))


def circular_delta_abs_deg(a_deg, b_deg):
    d = np.abs((a_deg - b_deg) % 360.0)
    return float(min(d, 360.0 - d))


def entropy_bits(p, eps=1e-12):
    p = np.clip(np.asarray(p,float), eps, 1.0)
    return float(-(p * (np.log(p)/np.log(2))).sum())


def equal_width_surrogate_bc(delta_pd_deg, kappa_bar, n_dirs=72):
    """
    Predicted BC if both neurons had SAME width kappa_bar, zero baseline, same amplitude.
    Captures ΔPD + width only (no baseline/tails).
    """
    th = np.linspace(0, 2*np.pi, n_dirs, endpoint=False)
    mu1 = 0.0
    mu2 = np.deg2rad(delta_pd_deg)
    y1 = np.exp(kappa_bar * np.cos(th - mu1)) / (2*np.pi*i0(kappa_bar))
    y2 = np.exp(kappa_bar * np.cos(th - mu2)) / (2*np.pi*i0(kappa_bar))
    p1, p2 = pmf_from_curve(y1), pmf_from_curve(y2)
    return bhattacharyya(p1, p2)


# ===================== Session indexing & curve fits ==================
def build_session_labels(masterGoodUnits):
    """Return date labels like 'Meetz_230607' / 'Akshan_24xxxx' for dates with >1 good unit."""
    dates = [u.split('_')[0] for u in masterGoodUnits]
    date_counts = Counter(dates)
    multi_unit_dates = [d for d, c in date_counts.items() if c > 1]
    labels = []
    for d in multi_unit_dates:
        if d.startswith("23"):
            labels.append(f"Meetz_{d}")
        elif d.startswith("24"):
            labels.append(f"Akshan_{d}")
        else:
            labels.append(d)
    return labels


def fit_dir_curves_for_units(base_dir, date_label, masterAllUnits, masterGoodUnits):
    """
    Returns: dict unit -> dict with keys:
      idx, angles_rad, yi, fit=(b,A,mu,k), yhat, R2
    Only for the good units present on that date.
    """
    monkeyName, sessionDate = date_label.split('_')
    sesh_dir = os.path.join(base_dir, monkeyName, sessionDate, "Direction Tuning")
    if not os.path.isdir(sesh_dir):
        return {}

    cwd0 = os.getcwd()
    try:
        os.chdir(sesh_dir)
        dirData = np.load('unitsDirTuningMat.npy')
    finally:
        os.chdir(cwd0)

    if dirData.ndim != 2 or dirData.shape[1] < N_DIRS_MIN:
        return {}

    units_today = [u for u in masterAllUnits  if u.startswith(sessionDate + "_")]
    good_today  = [u for u in masterGoodUnits if u.startswith(sessionDate + "_")]
    unit_to_row = {u: i for i, u in enumerate(units_today)}
    good_today  = [u for u in good_today if u in unit_to_row]
    if len(good_today) == 0:
        return {}

    n_neurons_today, n_dirs = dirData.shape
    angles_deg = np.linspace(0, 360, n_dirs, endpoint=False)
    angles_rad = np.deg2rad(angles_deg)

    out = {}
    for u in good_today:
        idx = unit_to_row[u]
        yi = dirData[idx, :].astype(float)
        b, A, mu, k = fit_von_mises_single(yi, angles_rad)
        if np.all(np.isfinite([b, A, mu, k])):
            yhat = von_mises_func(angles_rad, b, A, mu, k)
            R2 = r2_score_safe(yi, yhat)
        else:
            yhat, R2 = np.full_like(angles_rad, np.nan), np.nan
        out[u] = dict(idx=idx, angles_rad=angles_rad, yi=yi, fit=(b,A,mu,k), yhat=yhat, R2=R2)
    return out


# =================== Utilities to select pairs by BC ==================
def summarize_pairs_by_bc(corr_df, low_thresh=LOW_BC_THRESH, high_thresh=HIGH_BC_THRESH):
    df = corr_df.copy()
    df["pair_key"] = df.apply(lambda r: tuple(sorted((str(r["unit_i"]), str(r["unit_j"])))), axis=1)
    bc_med = df.groupby("pair_key")["BC_dir_tuning"].median()
    groups = {}
    for pk, bc in bc_med.items():
        if pd.notna(bc) and bc < low_thresh:
            groups[pk] = "low"
        elif pd.notna(bc) and bc > high_thresh:
            groups[pk] = "high"
    return groups  # dict: pair_key -> "low"/"high"


# ================= Build pair-level features (means & diffs) ==========
def build_features_for_pairs(base_dir, masterAllUnits, masterGoodUnits, pair_groups,
                             r2_min=R2_MIN):
    """
    For each selected pair, (re)fit both units from the session file and compute PAIR features:
      - deltaPD_deg
      - pair_mean_FWHM_deg, pair_absdiff_FWHM_deg
      - pair_mean_kappa, pair_absdiff_kappa
      - pair_mean_baselineFrac, pair_absdiff_baselineFrac
      - pair_mean_antiFrac, pair_absdiff_antiFrac
      - pair_mean_entropy_bits, pair_absdiff_entropy_bits
      - BC_surrogate_equalWidth (ΔPD + mean κ only)
      - BC_fit_fromVM (from PMFs of fitted curves)
      - BC_pair_median (observed, merged from corr_df)
    """
    # Cache per session to avoid refitting
    date_labels = build_session_labels(masterGoodUnits)
    fits_cache = {}

    rows = []
    for pk, grp in pair_groups.items():
        u1, u2 = pk
        date = u1.split('_')[0]
        if date != u2.split('_')[0]:
            continue  # safety

        monkeyName = "Meetz" if date.startswith("23") else ("Akshan" if date.startswith("24") else "Unknown")
        date_label = f"{monkeyName}_{date}"

        if date_label not in fits_cache:
            fits_cache[date_label] = fit_dir_curves_for_units(base_dir, date_label, masterAllUnits, masterGoodUnits)
        cur = fits_cache[date_label]
        if (u1 not in cur) or (u2 not in cur):
            continue

        d1, d2 = cur[u1], cur[u2]
        if not (np.isfinite(d1["R2"]) and np.isfinite(d2["R2"]) and d1["R2"] >= r2_min and d2["R2"] >= r2_min):
            continue

        angles = d1["angles_rad"]
        b1, A1, mu1, k1 = d1["fit"]; yhat1 = d1["yhat"]
        b2, A2, mu2, k2 = d2["fit"]; yhat2 = d2["yhat"]

        # --- Preferred-direction responses ---
        prefResp1 = float(np.max(yhat1))  # absolute peak
        prefResp2 = float(np.max(yhat2))
        prefResp1_above = prefResp1 - float(b1)  # peak above baseline
        prefResp2_above = prefResp2 - float(b2)

        pair_mean_prefResp = np.nanmean([prefResp1, prefResp2])
        pair_absdiff_prefResp = np.abs(prefResp1 - prefResp2)
        pair_mean_prefResp_minus_base = np.nanmean([prefResp1_above, prefResp2_above])
        pair_absdiff_prefResp_minus_base = np.abs(prefResp1_above - prefResp2_above)

        # PD / ΔPD
        pd1 = np.rad2deg(mu1) % 360.0
        pd2 = np.rad2deg(mu2) % 360.0
        dPD = circular_delta_abs_deg(pd1, pd2)

        # Widths / kappas
        fwhm1 = np.rad2deg(vm_fwhm_from_kappa(k1))
        fwhm2 = np.rad2deg(vm_fwhm_from_kappa(k2))
        pair_mean_FWHM = np.nanmean([fwhm1, fwhm2])
        pair_absdiff_FWHM = np.abs(fwhm1 - fwhm2)

        pair_mean_kappa = np.nanmean([k1, k2])
        pair_absdiff_kappa = np.abs(k1 - k2)

        # Baseline fraction & anti-pref fraction
        peak1 = float(np.max(yhat1)); peak2 = float(np.max(yhat2))
        anti1 = float(von_mises_func(mu1 + np.pi, b1, A1, mu1, k1))
        anti2 = float(von_mises_func(mu2 + np.pi, b2, A2, mu2, k2))
        antiFrac1 = anti1/peak1 if peak1 > 1e-12 else np.nan
        antiFrac2 = anti2/peak2 if peak2 > 1e-12 else np.nan
        baseFrac1 = b1/(b1 + A1 + 1e-12) if np.isfinite(b1+A1) else np.nan
        baseFrac2 = b2/(b2 + A2 + 1e-12) if np.isfinite(b2+A2) else np.nan

        pair_mean_baseFrac = np.nanmean([baseFrac1, baseFrac2])
        pair_absdiff_baseFrac = np.abs(baseFrac1 - baseFrac2)

        pair_mean_antiFrac = np.nanmean([antiFrac1, antiFrac2])
        pair_absdiff_antiFrac = np.abs(antiFrac1 - antiFrac2)

        # Entropy on fitted PMFs
        p1 = pmf_from_curve(yhat1)
        p2 = pmf_from_curve(yhat2)
        H1 = entropy_bits(p1); H2 = entropy_bits(p2)
        pair_mean_entropy = np.nanmean([H1, H2])
        pair_absdiff_entropy = np.abs(H1 - H2)

        # BCs
        bc_fit = bhattacharyya(p1, p2)
        bc_surr = equal_width_surrogate_bc(dPD, pair_mean_kappa, n_dirs=len(angles))

        rows.append({
            "pair_key": pk, "date": date, "monkey": monkeyName, "group": grp,
            "deltaPD_deg": dPD,
            "pair_mean_FWHM_deg": pair_mean_FWHM,
            "pair_absdiff_FWHM_deg": pair_absdiff_FWHM,
            "pair_mean_kappa": pair_mean_kappa,
            "pair_absdiff_kappa": pair_absdiff_kappa,
            "pair_mean_baselineFrac": pair_mean_baseFrac,
            "pair_absdiff_baselineFrac": pair_absdiff_baseFrac,
            "pair_mean_antiFrac": pair_mean_antiFrac,
            "pair_absdiff_antiFrac": pair_absdiff_antiFrac,
            "pair_mean_entropy_bits": pair_mean_entropy,
            "pair_absdiff_entropy_bits": pair_absdiff_entropy,
            "pair_mean_prefResp": pair_mean_prefResp,
            "pair_absdiff_prefResp": pair_absdiff_prefResp,
            "pair_mean_prefResp_minus_base": pair_mean_prefResp_minus_base,
            "pair_absdiff_prefResp_minus_base": pair_absdiff_prefResp_minus_base,
            "BC_fit_fromVM": bc_fit,
            "BC_surrogate_equalWidth": bc_surr,
            "R2_i": d1["R2"], "R2_j": d2["R2"]
        })

    feat_df = pd.DataFrame(rows)

    # Merge observed median BC from corr_df
    cdf = corr_df.copy()
    cdf["pair_key"] = cdf.apply(lambda r: tuple(sorted((str(r["unit_i"]), str(r["unit_j"])))), axis=1)
    bc_med = cdf.groupby("pair_key")["BC_dir_tuning"].median().rename("BC_pair_median")
    if len(feat_df):
        feat_df = feat_df.merge(bc_med.reset_index(), on="pair_key", how="left")

    return feat_df


# ============================== GUARDS ===============================
required_globals = ["base_dir", "masterAllUnits", "masterGoodUnits", "corr_df",
                    "fit_dir_curves_for_units", "build_session_labels",
                    "vm_fwhm_from_kappa", "pmf_from_curve", "bhattacharyya",
                    "equal_width_surrogate_bc", "circular_delta_abs_deg"]
for _name in required_globals:
    if _name not in globals():
        raise RuntimeError(f"Missing required variable `{_name}` or helper in the workspace.")

# =================== RUN: build features for low vs high ===================
pair_groups = summarize_pairs_by_bc(corr_df,
                                    low_thresh=LOW_BC_THRESH,
                                    high_thresh=HIGH_BC_THRESH)
print(f"Pairs selected -> low: {sum(np.array(list(pair_groups.values()))=='low')}, "
      f"high: {sum(np.array(list(pair_groups.values()))=='high')}")

features_df = build_features_for_pairs(base_dir,
                                       masterAllUnits,
                                       masterGoodUnits,
                                       pair_groups,
                                       r2_min=R2_MIN)

print(f"Feature table built with N={len(features_df)} pairs "
      f"({features_df['group'].value_counts().to_dict() if len(features_df) else {}})")

# ============================== PLOTTING ==============================
if len(features_df):
    # 1) BC vs ΔPD colored by pair_mean_FWHM
    plt.figure(figsize=(6, 5))
    sc = plt.scatter(features_df["deltaPD_deg"], features_df["BC_pair_median"],
                     c=features_df["pair_mean_FWHM_deg"], s=28, alpha=0.9, cmap="viridis")
    plt.colorbar(sc, label="Pair mean FWHM (deg)")
    plt.axvline(120, color='k', ls='--', alpha=0.4)
    plt.axhline(LOW_BC_THRESH, color='r', ls='--', alpha=0.5, label=f"BC={LOW_BC_THRESH}")
    plt.axhline(HIGH_BC_THRESH, color='g', ls='--', alpha=0.5, label=f"BC={HIGH_BC_THRESH}")
    plt.xlabel("ΔPD (deg)"); plt.ylabel("Observed BC (median)"); plt.legend()
    plt.title("Observed BC vs ΔPD (colored by pair mean width)")
    plt.tight_layout()

    # 2) Observed BC vs surrogate (ΔPD + equal width only)
    plt.figure(figsize=(6,5))
    x = features_df["BC_pair_median"]; y = features_df["BC_surrogate_equalWidth"]
    colors = np.where(features_df["group"]=="low", "tab:red", "tab:blue")
    plt.scatter(x, y, c=colors, s=30, alpha=0.8)
    lim = [0,1]; plt.plot(lim, lim, 'k--', alpha=0.5)
    plt.xlabel("Observed BC (median)"); plt.ylabel("Surrogate BC (ΔPD + mean κ, no baseline)")
    plt.title("How much of BC is explained by ΔPD + width?")
    plt.tight_layout()

    # 3) Boxplots by group for PAIR-LEVEL features
    def box_by_group(col, label):
        plt.figure(figsize=(5.2,4.2))
        data_low  = features_df.loc[features_df["group"]=="low",  col].dropna()
        data_high = features_df.loc[features_df["group"]=="high", col].dropna()
        bp = plt.boxplot([data_low, data_high], labels=["low BC", "high BC"],
                         patch_artist=True, widths=0.6)
        for patch, color in zip(bp['boxes'], ["#f28e8e", "#8eb3f2"]):
            patch.set_facecolor(color)
        plt.ylabel(label); plt.title(label + " (low vs high BC)"); plt.tight_layout()

    pair_cols = [
        ("deltaPD_deg", "ΔPD (deg)"),
        ("pair_mean_FWHM_deg", "Pair mean FWHM (deg)"),
        ("pair_absdiff_FWHM_deg", "Pair |Δ| FWHM (deg)"),
        ("pair_mean_kappa", "Pair mean κ"),
        ("pair_absdiff_kappa", "Pair |Δ| κ"),
        ("pair_mean_baselineFrac", "Pair mean baseline fraction"),
        ("pair_absdiff_baselineFrac", "Pair |Δ| baseline fraction"),
        ("pair_mean_antiFrac", "Pair mean anti-pref/peak"),
        ("pair_absdiff_antiFrac", "Pair |Δ| anti-pref/peak"),
        ("pair_mean_entropy_bits", "Pair mean entropy (bits)"),
        ("pair_absdiff_entropy_bits", "Pair |Δ| entropy (bits)"),
        ("BC_fit_fromVM", "BC from fitted VM"),
        ("BC_surrogate_equalWidth", "BC surrogate (ΔPD + mean κ)"),
        ("pair_mean_prefResp", "Pair mean preferred response"),
        ("pair_mean_prefResp_minus_base", "Pair mean preferred response (above baseline)"),
    ]
    for col, label in pair_cols:
        if col in features_df.columns:
            box_by_group(col, label)

    plt.show()

    # 4) Quick numeric summary table (means ± SEM per group)
    def summarize(cols):
        out = {}
        for g, gdf in features_df.groupby("group"):
            for c in cols:
                v = gdf[c].dropna().values
                if v.size:
                    out[(g,c)] = (np.mean(v), np.std(v, ddof=1)/np.sqrt(v.size))
        return out

    cols_to_report = ["deltaPD_deg", "pair_mean_FWHM_deg", "pair_mean_kappa",
                      "pair_mean_baselineFrac", "pair_mean_antiFrac",
                      "pair_mean_entropy_bits", "pair_mean_prefResp",
                      "pair_mean_prefResp_minus_base",
                      "BC_surrogate_equalWidth"]
    summary = summarize(cols_to_report)
    print("\n=== Pair-level feature means ± SEM ===")
    for (g, c), (mu, se) in summary.items():
        print(f"{g:>4} | {c:>26}: {mu:8.3f} ± {se:6.3f}")

# The main output DataFrame you’ll use downstream:
# features_df  (one row per PAIR with columns shown above)
# If useful:
# features_df.to_csv("pair_features_low_vs_high_pairLevel.csv", index=False)

# =============================== STAT TESTS =======================================


########################################################################################################################
########################################### Population Average Von Mises Fit ###########################################
########################################################################################################################


# -------------- helper: normalize a tuning curve against its baseline --------------
def baseline_peak_normalize(yhat, baseline):
    """
    Return curve in [0,1] where 0 = baseline, 1 = peak above baseline.
    If peak==baseline (flat), returns zeros.
    """
    peak = np.nanmax(yhat)
    denom = (peak - baseline)
    if denom <= 1e-12:
        return np.zeros_like(yhat)
    return (yhat - baseline) / denom


# Hyperparams
R2_MIN_FOR_SUPER = 0.00   # e.g., 0.7 or 0.9 if you want to quality-gate the population plot
N_GRID = 361             # resolution for the relative-angle grid

# Relative angle grid centered at PD=0 (in radians)
rel_grid = np.linspace(-np.pi, np.pi, N_GRID, endpoint=True)

# Accumulators across ALL sessions
super_curves = []        # list of arrays, each length = N_GRID
kept_meta = []           # (date, unit_id, r2) for reference/debug
kept_PD_rad = []   # absolute PDs (radians) for units that pass the R^2 gate
kept_kappa = []

# ========== Pass 2 over sessions to build super curves (reuse the fits you already compute) ==========
# If you still have the loop variables in scope, you can replicate the fitting part minimally.
# Otherwise, we can refit here quickly using dirData and the same helpers.

# Recreate the date/session iteration (same as above)
super_date_labels = date_labels  # reuse from your earlier block

for date in super_date_labels:
    monkeyName, sessionDate = date.split('_')
    sesh_dir = os.path.join(base_dir, monkeyName, sessionDate, "Direction Tuning")
    os.chdir(sesh_dir)

    dirData = np.load('unitsDirTuningMat.npy')
    units_today = [u for u in masterAllUnits if u.startswith(sessionDate + "_")]
    goodUnits_today = [u for u in masterGoodUnits if u.startswith(sessionDate + "_")]

    unit_to_row = {u: i for i, u in enumerate(units_today)}
    goodUnits_today = [u for u in goodUnits_today if u in unit_to_row]
    if len(goodUnits_today) == 0:
        os.chdir(base_dir)
        continue

    # Angle grid for the raw data
    n_neurons_today, n_dirs = dirData.shape
    angles_deg = np.linspace(0, 360, n_dirs, endpoint=False)
    angles_rad = np.deg2rad(angles_deg)

    for u in goodUnits_today:
        idx = unit_to_row[u]
        yi = dirData[idx, :].astype(float)

        # Fit the von Mises (reuse your helper)
        b, A, prefD, k = fit_von_mises_single(yi, angles_rad)
        if not (np.isfinite(b) and np.isfinite(A) and np.isfinite(prefD) and np.isfinite(k)):
            continue

        # Compute R^2 for gating
        yhat = von_mises_func(angles_rad, b, A, prefD, k)
        r2 = r2_score_safe(yi, yhat)
        if not np.isfinite(r2) or r2 < R2_MIN_FOR_SUPER:
            continue

        # Evaluate the fitted curve on a fine grid *relative to PD* (PD aligned to 0)
        # i.e., model at absolute angles (prefD + rel_grid)
        yfit_fine = von_mises_func(prefD + rel_grid, b, A, prefD, k)

        # Normalize: baseline->0, peak->1 (unit area not enforced; this is peak normalization)
        yfit_fine_norm = baseline_peak_normalize(yfit_fine, baseline=b)

        # Save absolute PD for the rose plot (only if it passes the R² gate)
        kept_PD_rad.append(float(prefD))
        kept_kappa.append(float(k))

        super_curves.append(yfit_fine_norm)
        kept_meta.append((sessionDate, u, r2))

    os.chdir(base_dir)

# Stack and compute population stats
if len(super_curves) == 0:
    print("No curves passed the filter for the super plot.")
else:
    super_mat = np.vstack(super_curves)                 # (n_units_kept, N_GRID)
    mean_curve = np.nanmean(super_mat, axis=0)
    sem_curve  = np.nanstd(super_mat, axis=0, ddof=1) / np.sqrt(super_mat.shape[0])

    # ========================= Plot 1: Spaghetti + mean ± SEM in Cartesian =========================
    fig, ax = plt.subplots(figsize=(7, 4))
    # plot all normalized, aligned curves (light)
    for row in super_mat:
        ax.plot(np.rad2deg(rel_grid), row, alpha=0.08, linewidth=1)

    # mean ± sem
    ax.plot(np.rad2deg(rel_grid), mean_curve, linewidth=2.5, label='Population mean')
    ax.fill_between(np.rad2deg(rel_grid), mean_curve - sem_curve, mean_curve + sem_curve, alpha=0.25, label='± SEM')

    ax.set_xlabel("Angle relative to PD (deg)")
    ax.set_ylabel("Normalized response (baseline→0, peak→1)")
    ax.set_title(f"Population-aligned von Mises (n={super_mat.shape[0]} units; R² ≥ {R2_MIN_FOR_SUPER})")
    ax.set_xlim([-180, 180])
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    plt.tight_layout()

    # ========================= Plot 2: Polar plot of the population mean =========================
    # Build a closed curve for polar (wrap ends to the same point)
    theta_polar = np.mod(rel_grid, 2*np.pi)
    order = np.argsort(theta_polar)
    theta_sorted = theta_polar[order]
    mean_sorted  = mean_curve[order]
    sem_sorted   = sem_curve[order]

    fig2 = plt.figure(figsize=(5.5, 5.5))
    axp = fig2.add_subplot(111, projection='polar')
    axp.plot(theta_sorted, mean_sorted, linewidth=2.5, label='Population mean')

    # Optional: add SEM guides (two dashed lines). Polar fill_between is tricky; lines are cleaner.
    axp.plot(theta_sorted, mean_sorted + sem_sorted, linestyle='--', linewidth=1.2, alpha=0.6, label='± SEM')
    axp.plot(theta_sorted, np.maximum(0, mean_sorted - sem_sorted), linestyle='--', linewidth=1.2, alpha=0.6)

    axp.set_title("Population-aligned von Mises (polar)", va='bottom')
    axp.set_theta_zero_location('E')  # PD at 0° on the right
    axp.set_theta_direction(-1)       # increasing clockwise (optional; set to +1 for CCW)
    axp.set_rlabel_position(90)       # move radial labels out of the way
    axp.set_rticks([0.25, 0.5, 0.75, 1.0])
    axp.grid(True, alpha=0.3)
    axp.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

    plt.tight_layout()

    # If you want to keep the population arrays for later analysis:
    population_super = {
        "rel_grid_rad": rel_grid,
        "mean_curve": mean_curve,
        "sem_curve": sem_curve,
        "n_units": super_mat.shape[0]
    }

# ========================= Plot 3: Polar "rose" of preferred directions =========================
if len(kept_PD_rad) > 0:
    PD = np.asarray(kept_PD_rad, dtype=float)

    # Number of bins for the rose (e.g., 24 bins = 15° each)
    NBINS = 24
    bins = np.linspace(0.0, 2*np.pi, NBINS + 1, endpoint=True)

    counts, edges = np.histogram(PD % (2*np.pi), bins=bins)
    theta_centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    # Normalize if you prefer density; set to False to plot raw counts
    as_density = False
    heights = counts.astype(float)
    if as_density:
        area = heights.sum()
        if area > 0:
            heights = heights / area

    fig3 = plt.figure(figsize=(5.8, 5.8))
    axr = fig3.add_subplot(111, projection='polar')

    # Polar bar plot (rose)
    axr.bar(theta_centers, heights, width=widths, align='center', alpha=0.8, edgecolor='none')

    # Circular mean vector
    C = np.cos(PD).mean()
    S = np.sin(PD).mean()
    mean_dir = np.arctan2(S, C) % (2*np.pi)
    R = np.hypot(C, S)  # resultant length in [0,1]

    # Scale the mean vector to the plot’s max radius for visibility
    r_max = (heights.max() if heights.size > 0 else 1.0)
    axr.plot([mean_dir, mean_dir], [0, r_max], linewidth=3, label=f"Mean dir (R={R:.2f})")

    axr.set_title(f"Preferred-direction distribution (n={len(PD)})", va='bottom')
    axr.set_theta_zero_location('E')  # 0° on the right
    axr.set_theta_direction(-1)       # clockwise increase (optional)
    # axr.set_rticks([])                # hide radial ticks for a cleaner rose
    axr.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
else:
    print("No PDs passed the R² gate; skipping PD rose plot.")

plt.show()

if len(kept_kappa) > 0:
    mean_kappa = np.mean(kept_kappa)
    median_kappa = np.median(kept_kappa)
    print(f"Average κ across population: {mean_kappa:.2f}")
    print(f"Median κ across population: {median_kappa:.2f}")
else:
    print("No κ values collected.")

########################################################################################################################
################################################# extra code ###########################################################
########################################################################################################################

import numpy as np
import matplotlib.pyplot as plt

# === Filter for your condition ===
filterVal = 0.1
mask = (
    (corr_df["singleWeightGeom_center"] > filterVal) &
    (corr_df["singleWeightGeom_periph"] > filterVal) &
    (corr_df["pair_weight_geommean"] > filterVal)
)

filtered = corr_df[mask]

# === Concatenate n1w and n2w ===
weights = np.concatenate([filtered["n1w"].to_numpy(),
                          filtered["n2w"].to_numpy()])

# === Restrict to values < 5 ===
weights = weights[weights < 1]

print(f"Total weights used: {len(weights)}")

# === Plot distribution ===
plt.figure(figsize=(7, 5))
plt.hist(weights, bins=30, alpha=0.7, color="steelblue", edgecolor="black")
plt.xlabel("Weight value (< 5)")
plt.ylabel("Count")
plt.title("Distribution of n1w and n2w (concatenated)")
plt.show()

