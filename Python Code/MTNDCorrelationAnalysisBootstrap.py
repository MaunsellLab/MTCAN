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
# NEW: sliding-window spike counts (per unit: [n_windows, blocksDone, nConds])
masterCorrUnitsSpikeCountMatSliding = []
slidingWindowsMS = None  # will be set per session (and should be consistent)

# ================= GLOBAL PSTH AGGREGATION (across sessions) =================
# Will accumulate normalized and raw PSTHs across (good unit × condition)
masterPSTH_norm_list = []     # each entry: [T] normalized PSTH for one (unit, stimIndex)
masterPSTH_raw_list  = []     # each entry: [T] raw spikes/s PSTH for one (unit, stimIndex)

# Max firing rates (spikes/s) for each (unit, stimIndex)
masterPSTH_maxFR = []         # float
masterPSTH_unitKeys = []      # string like '230607_149'
masterPSTH_stimIdx = []       # int stimIndex

# (optional) store metadata
masterPSTH_histPrePostMS = None
masterPSTH_trueStimDurMS = None
# ============================================================================

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

    # ================= NEW: SLIDING WINDOWS =================
    # Settings (ms)
    SLIDE_ONSET_MS = 50  # onset latency
    SLIDE_WIN_MS = 80  # fixed window length
    SLIDE_STEP_MS = 20  # step size
    # stimulus + offLatency locked end time (ms relative to stim onset)
    offLatencyMS = int(round(offLatency * 1000))
    maxEndMS = int(trueStimDurMS) + offLatencyMS  # windows can extend into offLatency
    # window starts (ms) from 50ms to maxEndMS - 50ms inclusive, step 10ms
    startsMS = np.arange(SLIDE_ONSET_MS, maxEndMS - SLIDE_WIN_MS + 1, SLIDE_STEP_MS, dtype=int)
    endsMS = startsMS + SLIDE_WIN_MS
    nWin = len(startsMS)
    # Save window definitions (ms) for later use/metadata
    slidingWindowsMS = np.column_stack([startsMS, endsMS])  # shape (nWin, 2)
    # Allocate sliding-window spike count matrix:
    # shape = (nWin, nUnits, blocksDone+1, nConds)
    spikeCountMatSliding = np.zeros((nWin, len(units), blocksDone + 1, numSteps * 7 + 3), dtype=np.float32)
    # =======================================================

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

                            # ================= NEW: SLIDING WINDOW COUNTS =================
                            # Compute once per stim presentation (speed): absolute edges in seconds
                            # Use stimulus+offLatency locked windows (relative to stimOnTimeS)
                            # Count spikes in each window [start, end] (same inclusive convention as your existing code)
                            for w in range(nWin):
                                w_start_s = stimOnTimeS + (startsMS[w] / 1000.0)
                                w_end_s = stimOnTimeS + (endsMS[w] / 1000.0)

                                # Count spikes in this sliding window
                                stimSpikes_w = np.where((unitTimeStamps >= w_start_s) &
                                                        (unitTimeStamps <= w_end_s))[0]
                                spikeCountMatSliding[w, unitCount, stCount, stimIndex] = len(stimSpikes_w)
                            # ==============================================================

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
            # NEW: sliding windows -> per unit array shape (nWin, blocksDone, nConds)
            masterCorrUnitsSpikeCountMatSliding.append(spikeCountMatSliding[:, idx, :blocksDone, :])

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

        # ================= GRAND PSTH (blocksDone division; 1 presentation/cond/block) =================
        nConds = numSteps * 7 + 3

        # Save metadata once (assumes consistent across sessions; if not, you’ll want to resample)
        if masterPSTH_histPrePostMS is None:
            masterPSTH_histPrePostMS = histPrePostMS
            masterPSTH_trueStimDurMS = int(trueStimDurMS)

        stim_start = histPrePostMS
        stim_end = histPrePostMS + int(trueStimDurMS)

        for uid in corrUnits:
                uidx = np.where(units == uid)[0][0]
                unitKey = f"{seshDate}_{uid}"

                for stimIdx in range(nConds):
                    # counts per 1 ms bin across all blocks (since 1 presentation per block)
                    hist_counts = spikeHists[uidx, stimIdx, :].astype(np.float32)

                    # Convert -> spikes/s (per-bin)
                    psth_rate = hist_counts * 1000.0 / float(blocksDone)

                    # max FR within stimulus window (exclude pre-period)
                    if stim_end <= psth_rate.size:
                        maxFR = float(np.max(psth_rate[stim_start:stim_end]))
                    else:
                        maxFR = float(np.max(psth_rate))

                    masterPSTH_maxFR.append(maxFR)
                    masterPSTH_unitKeys.append(unitKey)
                    masterPSTH_stimIdx.append(int(stimIdx))

                    # Normalize by per-(unit,condition) max
                    if maxFR > 0:
                        psth_norm = psth_rate / maxFR
                    else:
                        psth_norm = np.full_like(psth_rate, np.nan, dtype=np.float32)

                    masterPSTH_raw_list.append(psth_rate.astype(np.float32))
                    masterPSTH_norm_list.append(psth_norm.astype(np.float32))
        # =============================================================================

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

# ================= FINAL GRAND AVERAGES =================
masterPSTH_norm_arr = np.array(masterPSTH_norm_list)  # shape: [NunitCond, T]
masterPSTH_raw_arr  = np.array(masterPSTH_raw_list)   # shape: [NunitCond, T]

grandAvgPSTH_norm = np.nanmean(masterPSTH_norm_arr, axis=0)  # [T]
grandAvgPSTH_raw  = np.nanmean(masterPSTH_raw_arr,  axis=0)  # [T]

print("Grand PSTH computed from N (unit×condition) =", masterPSTH_norm_arr.shape[0])
print("PSTH length T =", grandAvgPSTH_norm.shape[0])

# Also useful summaries of max firing rates
masterPSTH_maxFR = np.array(masterPSTH_maxFR, dtype=float)
print("Mean max FR (sp/s) across unit×condition =", np.nanmean(masterPSTH_maxFR))
print("Median max FR (sp/s) across unit×condition =", np.nanmedian(masterPSTH_maxFR))
# =========================================================

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
############################################## Standard Normalization Model  ###########################################
########################################################################################################################

for filler in range(1):
    """
    standard_norm_model.py
    ======================
    Fits the standard (Carandini & Heeger-style) normalization model to each neuron
    alongside the existing RF-weighted heuristic model, then compares R².
    
    Standard model:
        R = (c_center * L_loc_center + c_periph * L_loc_periph) / 
            (a_center * c_center + a_periph * c_periph + sigma) + b
    
    Where:
      - Each of the 5 locations (0=center, 1-4=periphery) has its own Lp and Lnp
      - a_center (loc 0) is FIXED to 1.0; a_1..a_4 are free (4 free alphas)
      - One shared sigma across all locations
      - One shared baseline b
    
    Free parameters (15 total):
        Lp_0, Lnp_0                    — center location drives (2)
        Lp_1, Lnp_1, ..., Lp_4, Lnp_4 — peripheral location drives (8)
        a_1, a_2, a_3, a_4             — peripheral suppressive weights (4)
        sigma                           — semi-saturation constant (1)
        b                               — baseline (1)
    
    Stim encoding (same as heuristic model):
        contrast_centerSpon, contrast_peripherySpon  — 0 or 1
        locationsSpon        — -1 (no periph), 0..3 (periph site index)
        stim_type_centerSpon — 1=pref, 0=nonpref, -1=no stim
        stim_type_peripherySpon — 1=pref, 0=nonpref, -1=no stim
    
    Requires in scope (already built from main pipeline):
        prefNormalizedWithBase, nonprefNormalizedWithBase,
        pnNormalizedWithBase, npNormalizedWithBase,
        ppNormalizedWithBase, nnNormalizedWithBase,
        sponNormalizedToPrefCenter,
        contrast_centerSpon, contrast_peripherySpon,
        locationsSpon, stim_type_centerSpon, stim_type_peripherySpon,
        rfWeightNeuronFits, r2ScoresRFWeightHeuristic  (for comparison)
    """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from sklearn.metrics import r2_score
    import warnings
    warnings.filterwarnings("ignore")


    # ══════════════════════════════════════════════════════════════════════════════
    #  1. MODEL DEFINITION
    # ══════════════════════════════════════════════════════════════════════════════

    def fullStandardNorm(contrast_center, contrast_periphery, location,
                         stim_type_center, stim_type_periphery,
                         Lp_0, Lnp_0,
                         Lp_1, Lnp_1,
                         Lp_2, Lnp_2,
                         Lp_3, Lnp_3,
                         Lp_4, Lnp_4,
                         a_1, a_2, a_3, a_4,
                         sigma, b):
        """
        Standard normalization model with location-specific L terms.

        Location encoding:
            location == -1  : no peripheral stimulus (center only or blank)
            location == 0   : peripheral site 1
            location == 1   : peripheral site 2
            location == 2   : peripheral site 3
            location == 3   : peripheral site 4

        Center is always 'location 0' in the L/alpha sense (alpha fixed = 1).

        Parameters
        ----------
        Lp_0, Lnp_0        : drive terms for center (loc index 0 in stimulus space)
        Lp_1..Lp_4         : drive terms for peripheral sites 1..4
        a_1..a_4            : suppressive weights for peripheral sites (a_center = 1)
        sigma               : semi-saturation constant (shared)
        b                   : additive baseline
        """
        loc = int(location)

        # ── Center drive ──────────────────────────────────────────────────────────
        stc = int(stim_type_center)
        if stc == 1:
            L_center = Lp_0
        elif stc == 0:
            L_center = Lnp_0
        else:
            L_center = 0.0

        # ── Peripheral drive and suppression ─────────────────────────────────────
        if loc == -1:
            # No peripheral stimulus
            L_periph = 0.0
            a_periph = 0.0
        else:
            # Location-specific Lp/Lnp
            Lp_locs  = [Lp_1,  Lp_2,  Lp_3,  Lp_4]
            Lnp_locs = [Lnp_1, Lnp_2, Lnp_3, Lnp_4]
            a_locs   = [a_1,   a_2,   a_3,   a_4  ]

            stp = int(stim_type_periphery)
            if stp == 1:
                L_periph = Lp_locs[loc]
            elif stp == 0:
                L_periph = Lnp_locs[loc]
            else:
                L_periph = 0.0

            a_periph = a_locs[loc]

        # ── Response ──────────────────────────────────────────────────────────────
        # a_center fixed = 1.0
        numerator   = (contrast_center  * L_center +
                       contrast_periphery * L_periph) + (sigma * b)
        denominator = (contrast_center  * 1.0 +
                       contrast_periphery * a_periph + sigma)

        # Guard against division by zero
        if abs(denominator) < 1e-12:
            return b

        return numerator / denominator


    def model_wrapper_standard(data, Lp_0, Lnp_0,
                                Lp_1, Lnp_1,
                                Lp_2, Lnp_2,
                                Lp_3, Lnp_3,
                                Lp_4, Lnp_4,
                                a_1, a_2, a_3, a_4,
                                sigma, b):
        """Wrapper for scipy curve_fit — applies model element-wise."""
        contrast_c, contrast_p, loc, stc, stp = data
        return [fullStandardNorm(
                    float(contrast_c[k]), float(contrast_p[k]),
                    float(loc[k]), float(stc[k]), float(stp[k]),
                    Lp_0, Lnp_0,
                    Lp_1, Lnp_1,
                    Lp_2, Lnp_2,
                    Lp_3, Lnp_3,
                    Lp_4, Lnp_4,
                    a_1, a_2, a_3, a_4,
                    sigma, b)
                for k in range(len(contrast_c))]


    # ══════════════════════════════════════════════════════════════════════════════
    #  2. INITIAL GUESS BUILDER
    # ══════════════════════════════════════════════════════════════════════════════

    def build_initial_guess_standard(neuron_idx, prefNormWithBase, nonprefNormWithBase,
                                      heuristic_fit=None):
        """
        Build sensible initial guesses for the standard model.

        Strategy:
          Lp_0  : pref response at center (prefNormWithBase[:, 0])
          Lnp_0 : nonpref response at center
          Lp_k  : pref response at periph site k (prefNormWithBase[:, k])
          Lnp_k : nonpref response at periph site k
          a_k   : 1.0 (neutral start)
          sigma : 0.1 (small semi-saturation)
          b     : spontaneous rate
        """
        pref   = prefNormWithBase[neuron_idx]     # shape (5,): center + 4 periph
        nonp   = nonprefNormWithBase[neuron_idx]  # shape (5,)

        Lp_0   = max(float(pref[0]),  0.01)
        Lnp_0  = max(float(nonp[0]), 0.01)
        Lp_1   = max(float(pref[1]),  0.01)
        Lnp_1  = max(float(nonp[1]), 0.01)
        Lp_2   = max(float(pref[2]),  0.01)
        Lnp_2  = max(float(nonp[2]), 0.01)
        Lp_3   = max(float(pref[3]),  0.01)
        Lnp_3  = max(float(nonp[3]), 0.01)
        Lp_4   = max(float(pref[4]),  0.01)
        Lnp_4  = max(float(nonp[4]), 0.01)

        # Seed alphas from heuristic W if available, else 1.0
        if heuristic_fit is not None:
            a_1 = max(float(heuristic_fit[2]), 0.01)
            a_2 = max(float(heuristic_fit[3]), 0.01)
            a_3 = max(float(heuristic_fit[4]), 0.01)
            a_4 = max(float(heuristic_fit[5]), 0.01)
        else:
            a_1 = a_2 = a_3 = a_4 = 1.0

        sigma = max(float(heuristic_fit[6]), 0.01) if heuristic_fit is not None else 0.10
        b     = max(float(heuristic_fit[7]), 0.0)  if heuristic_fit is not None else 0.05

        return [Lp_0, Lnp_0,
                Lp_1, Lnp_1,
                Lp_2, Lnp_2,
                Lp_3, Lnp_3,
                Lp_4, Lnp_4,
                a_1, a_2, a_3, a_4,
                sigma, b]


    # ══════════════════════════════════════════════════════════════════════════════
    #  3. FIT LOOP
    # ══════════════════════════════════════════════════════════════════════════════

    def fit_standard_norm_population(prefNormalizedWithBase,
                                      nonprefNormalizedWithBase,
                                      pnNormalizedWithBase,
                                      npNormalizedWithBase,
                                      ppNormalizedWithBase,
                                      nnNormalizedWithBase,
                                      sponNormalizedToPrefCenter,
                                      contrast_centerSpon,
                                      contrast_peripherySpon,
                                      locationsSpon,
                                      stim_type_centerSpon,
                                      stim_type_peripherySpon,
                                      rfWeightNeuronFits=None):
        """
        Fit the standard normalization model to every neuron.

        Returns
        -------
        stdNormFits    : list of fitted parameter arrays (length = n_neurons)
                         Each entry is a 16-element array:
                         [Lp_0, Lnp_0, Lp_1, Lnp_1, ..., Lp_4, Lnp_4,
                          a_1, a_2, a_3, a_4, sigma, b]
                         or np.full(16, np.nan) if fit fails.
        r2ScoresStdNorm : np.ndarray of R² values per neuron
        """
        n_neurons = len(prefNormalizedWithBase)

        # Convert stim vectors to numpy arrays once
        cc  = np.asarray(contrast_centerSpon,    dtype=float)
        cp  = np.asarray(contrast_peripherySpon, dtype=float)
        loc = np.asarray(locationsSpon,          dtype=float)
        stc = np.asarray(stim_type_centerSpon,   dtype=float)
        stp = np.asarray(stim_type_peripherySpon, dtype=float)
        data = (cc, cp, loc, stc, stp)

        # Parameter bounds: all L terms >= 0, all a >= 0, sigma >= 1e-4, b free
        n_params = 16
        lower = [0.0] * 10 + [0.0] * 4 + [1e-4, -np.inf]   # 10 L + 4 a + sigma + b
        upper = [np.inf] * 16
        bounds = (lower, upper)

        stdNormFits     = []
        r2ScoresStdNorm = []

        for i in range(n_neurons):
            # Build response vector (same concatenation order as heuristic fit)
            resp = np.concatenate([
                prefNormalizedWithBase[i],
                nonprefNormalizedWithBase[i],
                pnNormalizedWithBase[i],
                npNormalizedWithBase[i],
                ppNormalizedWithBase[i],
                nnNormalizedWithBase[i],
                sponNormalizedToPrefCenter[i],
            ], axis=0)

            # Initial guess
            hfit = rfWeightNeuronFits[i] if rfWeightNeuronFits is not None else None
            p0 = build_initial_guess_standard(
                i, prefNormalizedWithBase, nonprefNormalizedWithBase, hfit
            )

            try:
                popt, _ = curve_fit(
                    model_wrapper_standard, data, resp,
                    p0=p0, bounds=bounds,
                    maxfev=20000, method="trf"
                )
                y_pred = model_wrapper_standard(data, *popt)
                r2 = r2_score(resp, y_pred)
                stdNormFits.append(np.array(popt))
                r2ScoresStdNorm.append(float(r2))

            except Exception as e:
                stdNormFits.append(np.full(n_params, np.nan))
                r2ScoresStdNorm.append(np.nan)
                if i < 5:   # only print first few failures to avoid spam
                    print(f"  Neuron {i}: fit failed — {e}")

            if (i + 1) % 25 == 0:
                print(f"  Fitted {i+1}/{n_neurons}  "
                      f"median R²={np.nanmedian(r2ScoresStdNorm):.3f}")

        r2ScoresStdNorm = np.array(r2ScoresStdNorm)
        print(f"\n  Done. Median R² (standard norm) = {np.nanmedian(r2ScoresStdNorm):.4f}")
        return stdNormFits, r2ScoresStdNorm


    # ══════════════════════════════════════════════════════════════════════════════
    #  4. COMPARISON & PLOTS
    # ══════════════════════════════════════════════════════════════════════════════

    def compare_model_fits(r2_heuristic, r2_standard, masterGoodUnits=None):
        """
        Compare R² distributions between heuristic and standard models.

        Returns a summary DataFrame with one row per neuron.
        """
        r2_h = np.asarray(r2_heuristic, dtype=float)
        r2_s = np.asarray(r2_standard,  dtype=float)

        delta = r2_s - r2_h   # positive = standard is better

        summary = pd.DataFrame({
            "unit":        masterGoodUnits if masterGoodUnits is not None
                           else np.arange(len(r2_h)),
            "r2_heuristic": r2_h,
            "r2_standard":  r2_s,
            "delta_r2":     delta,
        })

        print("\n" + "═"*55)
        print("  MODEL COMPARISON: R² (standard − heuristic)")
        print("═"*55)
        for label, vals in [("Heuristic (RF-weighted)", r2_h),
                             ("Standard  (loc-specific L)", r2_s)]:
            fin = vals[np.isfinite(vals)]
            print(f"  {label}:")
            print(f"    median={np.median(fin):.4f}  "
                  f"mean={np.mean(fin):.4f}  "
                  f"std={np.std(fin):.4f}  "
                  f"n={len(fin)}")

        fin_d = delta[np.isfinite(delta)]
        pct_std_better = 100 * (fin_d > 0).mean()
        from scipy.stats import wilcoxon
        stat, p = wilcoxon(fin_d)
        print(f"\n  ΔR² (standard − heuristic):")
        print(f"    median={np.median(fin_d):+.4f}  "
              f"mean={np.mean(fin_d):+.4f}")
        print(f"    Standard better in {pct_std_better:.1f}% of neurons")
        print(f"    Wilcoxon signed-rank: W={stat:.1f}, p={p:.4e}")
        print("═"*55)

        return summary


    def plot_r2_comparison(r2_heuristic, r2_standard):
        """Three-panel comparison plot."""
        r2_h = np.asarray(r2_heuristic, dtype=float)
        r2_s = np.asarray(r2_standard,  dtype=float)
        ok   = np.isfinite(r2_h) & np.isfinite(r2_s)
        delta = r2_s[ok] - r2_h[ok]

        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

        # Panel A: scatter R²_heuristic vs R²_standard
        ax = axes[0]
        ax.scatter(r2_h[ok], r2_s[ok], s=20, alpha=0.6,
                   color="#4F8FC0", edgecolors="white", linewidths=0.3)
        lim = [min(r2_h[ok].min(), r2_s[ok].min()) - 0.02,
               max(r2_h[ok].max(), r2_s[ok].max()) + 0.02]
        ax.plot(lim, lim, "k--", lw=1, alpha=0.5, label="y = x")
        ax.set_xlabel("R²  (RF-weighted heuristic)", fontsize=10)
        ax.set_ylabel("R²  (standard, loc-specific L)", fontsize=10)
        ax.set_title("R² scatter\n(one dot = one neuron)", fontsize=10)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.legend(fontsize=8, frameon=False)
        ax.spines[["top","right"]].set_visible(False)

        # Panel B: overlapping histograms
        ax = axes[1]
        bins = np.linspace(min(r2_h[ok].min(), r2_s[ok].min()),
                           max(r2_h[ok].max(), r2_s[ok].max()), 30)
        ax.hist(r2_h[ok], bins=bins, alpha=0.55, color="#E07B3F",
                label=f"Heuristic  (med={np.median(r2_h[ok]):.3f})", edgecolor="white")
        ax.hist(r2_s[ok], bins=bins, alpha=0.55, color="#4F8FC0",
                label=f"Standard   (med={np.median(r2_s[ok]):.3f})", edgecolor="white")
        ax.axvline(np.median(r2_h[ok]), color="#E07B3F", lw=2, ls="--")
        ax.axvline(np.median(r2_s[ok]), color="#4F8FC0", lw=2, ls="--")
        ax.set_xlabel("R²", fontsize=10)
        ax.set_ylabel("Neuron count", fontsize=10)
        ax.set_title("R² distributions", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        ax.spines[["top","right"]].set_visible(False)

        # Panel C: ΔR² histogram
        ax = axes[2]
        ax.hist(delta, bins=25, color="#7D5BA6", alpha=0.75, edgecolor="white")
        ax.axvline(0, color="gray", lw=1, ls="--")
        ax.axvline(np.median(delta), color="#7D5BA6", lw=2, ls="--",
                   label=f"Median ΔR² = {np.median(delta):+.4f}")
        pct = 100*(delta > 0).mean()
        ax.set_xlabel("ΔR²  (standard − heuristic)", fontsize=10)
        ax.set_ylabel("Neuron count", fontsize=10)
        ax.set_title(f"Standard better in {pct:.1f}% of neurons", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        ax.spines[["top","right"]].set_visible(False)

        fig.suptitle("Standard vs RF-weighted heuristic normalization model fit",
                     fontsize=12, fontweight="bold", y=1.02)
        plt.tight_layout()
        plt.savefig("standard_vs_heuristic_r2.png", dpi=150, bbox_inches="tight")
        plt.show()
        return fig


# ══════════════════════════════════════════════════════════════════════════════
#  5. ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════
#
# Run after rfWeightNeuronFits and r2ScoresRFWeightHeuristic are already built:
#
stdNormFits, r2ScoresStdNorm = fit_standard_norm_population(
      prefNormalizedWithBase, nonprefNormalizedWithBase,
      pnNormalizedWithBase, npNormalizedWithBase,
      ppNormalizedWithBase, nnNormalizedWithBase,
      sponNormalizedToPrefCenter,
      contrast_centerSpon, contrast_peripherySpon,
      locationsSpon, stim_type_centerSpon, stim_type_peripherySpon,
      rfWeightNeuronFits=rfWeightNeuronFits,   # seeds alphas from heuristic W
  )

summary_df = compare_model_fits(
      r2ScoresRFWeightHeuristic, r2ScoresStdNorm,
      masterGoodUnits=masterGoodUnits
  )

fig = plot_r2_comparison(r2ScoresRFWeightHeuristic, r2ScoresStdNorm)


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

# # ========================= Iterate sessions & collect similarity =========================
# base_dir = "/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN"
#
# similarity_records = []  # will accumulate dicts across all sessions
# for date in date_labels:
#     monkeyName, sessionDate = date.split('_')
#     sesh_dir = os.path.join(base_dir, monkeyName, sessionDate, "Direction Tuning")
#     os.chdir(sesh_dir)
#
#     # shape: (n_units_today, n_dirs)
#     dirData = np.load('unitsDirTuningMat.npy')
#
#     # Units present today (strings like '230607_77'), and the good ones
#     # Assumes: masterAllUnits and masterGoodUnits exist in your environment
#     units_today = [u for u in masterAllUnits if u.startswith(sessionDate + "_")]
#     goodUnits_today = [u for u in masterGoodUnits if u.startswith(sessionDate + "_")]
#
#     # Map unit string -> row index in dirData
#     # ASSUMPTION: rows of dirData correspond to units_today in order.
#     unit_to_row = {u: i for i, u in enumerate(units_today)}
#
#     # Keep only good units that actually exist in today's data
#     goodUnits_today = [u for u in goodUnits_today if u in unit_to_row]
#
#     # Need at least a pair
#     if len(goodUnits_today) < 2:
#         os.chdir(base_dir)
#         continue
#
#     # Angle grid
#     n_neurons_today, n_dirs = dirData.shape
#     angles_deg = np.linspace(0, 360, n_dirs, endpoint=False)
#     angles_rad = np.deg2rad(angles_deg)
#
#     # ----- Fit von Mises for good units; store PDs and R^2 aligned to goodUnits_today -----
#     pref_dirs_rad = np.zeros(len(goodUnits_today), dtype=float)
#     r2_list = np.full(len(goodUnits_today), np.nan, dtype=float)
#
#     for gi, u in enumerate(goodUnits_today):
#         idx = unit_to_row[u]
#         yi = dirData[idx, :].astype(float)
#         b, A, prefD, k = fit_von_mises_single(yi, angles_rad)
#         pref_dirs_rad[gi] = prefD
#
#         # Predicted tuning & R^2 (NaN-safe)
#         if np.isfinite(b) and np.isfinite(A) and np.isfinite(prefD) and np.isfinite(k):
#             yhat = von_mises_func(angles_rad, b, A, prefD, k)
#             r2_list[gi] = r2_score_safe(yi, yhat)
#         else:
#             r2_list[gi] = np.nan
#
#     pref_dirs_deg = np.rad2deg(pref_dirs_rad)
#
#     # ----- Tuning similarity (cos Δθ): SI in [-1,1], SI01 in [0,1] -----
#     c = np.cos(pref_dirs_rad).reshape(-1, 1)
#     s = np.sin(pref_dirs_rad).reshape(-1, 1)
#     SI = (c @ c.T) + (s @ s.T)     # (ngood x ngood), includes self-diagonal=1
#     SI01 = 0.5 * (SI + 1.0)
#
#     # ----- Bhattacharyya on raw tuning curves -----
#     eps = 1e-12
#     pmfs = []
#     for u in goodUnits_today:
#         idx = unit_to_row[u]
#         y = dirData[idx, :].astype(float)
#         y = np.clip(y, 0.0, None)                     # ensure nonnegative
#         p = (y + eps) / (np.sum(y + eps))             # normalize to sum=1
#         pmfs.append(p)
#     pmfs = np.vstack(pmfs)
#
#     # ----- Record only unique pairs (i < j) -----
#     for i, j in combinations(range(len(goodUnits_today)), 2):
#         # Bhattacharyya coefficient & distance
#         bc = float(np.sum(np.sqrt(pmfs[i] * pmfs[j])))    # in [0,1]
#         bd = float(-np.log(max(bc, eps)))                 # >= 0
#
#         similarity_records.append({
#             "date": sessionDate,
#             "monkey": monkeyName,
#             "unit_i": goodUnits_today[i],
#             "unit_j": goodUnits_today[j],
#             "PD_i_deg": float(pref_dirs_deg[i]),
#             "PD_j_deg": float(pref_dirs_deg[j]),
#             "SI": float(SI[i, j]),
#             "SI01": float(SI01[i, j]),
#             "BC_dir_tuning": bc,      # Bhattacharyya coefficient (0..1)
#             "BDir": bd,               # Bhattacharyya distance (-log(BC))
#             "R2_i": float(r2_list[i]),  # NEW: von Mises fit R^2 for unit_i
#             "R2_j": float(r2_list[j])   # NEW: von Mises fit R^2 for unit_j
#         })
#
#     # Return to base directory before next session
#     os.chdir(base_dir)

'''
'''
# ========================= Iterate sessions & collect similarity =========================
# ================= Signal Correlation on Fitted Responses of Von Mises ===================
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
    fitted_curves = np.full((len(goodUnits_today), n_dirs), np.nan, dtype=float)

    for gi, u in enumerate(goodUnits_today):
        idx = unit_to_row[u]
        yi = dirData[idx, :].astype(float)
        b, A, prefD, k = fit_von_mises_single(yi, angles_rad)
        pref_dirs_rad[gi] = prefD

        # Predicted tuning & R^2 (NaN-safe)
        if np.isfinite(b) and np.isfinite(A) and np.isfinite(prefD) and np.isfinite(k):
            yhat = von_mises_func(angles_rad, b, A, prefD, k)
            r2_list[gi] = r2_score_safe(yi, yhat)
            fitted_curves[gi, :] = yhat
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

    # ----- Signal correlation: Pearson r on mean tuning curves -----
    # Each row of tuning_mat is the mean response vector for a good unit
    tuning_mat = np.array([dirData[unit_to_row[u], :].astype(float)
                           for u in goodUnits_today])  # shape: (ngood, n_dirs)

    # Pearson r between all pairs; np.corrcoef returns (ngood x ngood) matrix
    sig_corr_mat = np.corrcoef(tuning_mat)  # SIG CORR On RAW VALUES NaN-safe as long as no zero-variance rows
    sig_corr_mat = np.corrcoef(fitted_curves)  # SIG CORR ON FITTED VALUES


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
            "sig_corr": float(sig_corr_mat[i, j]),  # Pearson r on mean tuning curves
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


def build_sig_corr_lookup(similarity_records):
    sig_corr_lookup = {}
    for r in similarity_records:
        a, b = r["unit_i"], r["unit_j"]
        key = tuple(sorted((a, b)))
        sig_corr_lookup[key] = r.get("sig_corr")
    return sig_corr_lookup


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
import numpy as np

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
               n1Sig=None, n2Sig=None,
               pairCov=None, singleCCov=None, singlePCov=None,
               pairVar1=None, pairVar2=None,
               singleCenterVar1=None, singleCenterVar2=None,
               singlePeriphVar1=None, singlePeriphVar2=None,
               pairMeanFR_paired_hz=None,
               pairMeanFR_single_center_hz=None,
               pairMeanFR_single_periph_hz=None,
               sig_corr=None,
               site=None,
               center_single_idx=None,
               periph_single_idx=None,
               single_center_isPref=None,
               single_periph_isPref=None,
               paired_center_isPref=None,
               n1StdSupp=None,
               n2StdSupp=None,
               pairStdSuppGeom=None):

    u1, u2 = comb
    date = u1.split("_")[0]

    def _f(x):
        return None if x is None else float(x)

    def _i(x):
        return None if x is None else int(x)

    def _b(x):
        return None if x is None else bool(x)

    def _den(v1, v2):
        # denom = sqrt(var1*var2) if both valid and positive
        if (v1 is None) or (v2 is None):
            return None
        try:
            v1 = float(v1); v2 = float(v2)
        except Exception:
            return None
        if not (np.isfinite(v1) and np.isfinite(v2)) or (v1 <= 0) or (v2 <= 0):
            return None
        return float(np.sqrt(v1 * v2))

    corr_records.append({
        "date": date,
        "unit_i": u1,
        "unit_j": u2,
        "block": block,
        "cond_index": i,

        # --- correlations ---
        "corr": float(pairStimCorr),
        "corr_transient": (None if pairStimCorr_transient is None else float(pairStimCorr_transient)),
        "corr_sustained": (None if pairStimCorr_sustained is None else float(pairStimCorr_sustained)),

        "singleCorr_center": _f(singleCorr_center),
        "singleCorr_periph": _f(singleCorr_periph),
        "singleCorr_center_transient": _f(singleCorr_center_transient),
        "singleCorr_periph_transient": _f(singleCorr_periph_transient),
        "singleCorr_center_sustained": _f(singleCorr_center_sustained),
        "singleCorr_periph_sustained": _f(singleCorr_periph_sustained),

        # --- weights / similarity ---
        "pair_weight_geommean": float(pairWeightGeomMean),
        "singleWeightGeom_center": _f(singleWeightGeom_center),
        "singleWeightGeom_periph": _f(singleWeightGeom_periph),

        "n1Sel": float(n1Sel),
        "n2Sel": float(n2Sel),
        "pairSel_geom": float(np.sign(n1Sel) * np.sign(n2Sel) * np.sqrt(abs(n1Sel) * abs(n2Sel))),

        "similarity_SI": _f(sim_SI),
        "similarity_SI01": _f(sim_SI01),
        "delta_PD_abs_deg": _f(delta_PD_abs_deg),
        "avg_pair_weight_2w": _f(avg_pair_weight_2w),
        "monkey": (None if monkey is None else str(monkey)),

        "BC_dir_tuning": _f(BC_dir_tuning),
        "BDir": _f(BDir),
        "sig_corr": _f(sig_corr),

        "n1StdSupp": _f(n1StdSupp),
        "n2StdSupp": _f(n2StdSupp),
        "pairStdSuppGeom": _f(pairStdSuppGeom),

        # --- metadata ---
        "site": _i(site),
        "center_single_idx": _i(center_single_idx),
        "periph_single_idx": _i(periph_single_idx),
        "single_center_isPref": _b(single_center_isPref),
        "single_periph_isPref": _b(single_periph_isPref),
        "paired_center_isPref": _b(paired_center_isPref),

        # --- per-neuron weights/sigmas ---
        "n1w": _f(n1w),
        "n2w": _f(n2w),
        "pairMinShareGeom": _f(pairMinShareGeom),
        "n1MinShare": _f(n1MinShare),
        "n2MinShare": _f(n2MinShare),
        "n1Sigma": _f(n1Sig),
        "n2Sigma": _f(n2Sig),

        # --- covariance + variances (counts^2) ---
        "pairCov": _f(pairCov),
        "pairVar1": _f(pairVar1),
        "pairVar2": _f(pairVar2),
        "pairDenom": _den(pairVar1, pairVar2),

        "singleCenterCov": _f(singleCCov),
        "singleCenterVar1": _f(singleCenterVar1),
        "singleCenterVar2": _f(singleCenterVar2),
        "singleCenterDenom": _den(singleCenterVar1, singleCenterVar2),

        "singlePeriphCov": _f(singlePCov),
        "singlePeriphVar1": _f(singlePeriphVar1),
        "singlePeriphVar2": _f(singlePeriphVar2),
        "singlePeriphDenom": _den(singlePeriphVar1, singlePeriphVar2),

        # --- rates (Hz) ---
        "pairMeanFR_paired_hz": _f(pairMeanFR_paired_hz),
        "pairMeanFR_single_center_hz": _f(pairMeanFR_single_center_hz),
        "pairMeanFR_single_periph_hz": _f(pairMeanFR_single_periph_hz),
    })


def cov_var_mean_from_vectors(x, y):
    """
    x,y: 1D arrays (trials).
    Returns: cov12, var1, var2, mean1, mean2 using sample cov (ddof=1).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mean1 = np.nanmean(x)
    mean2 = np.nanmean(y)

    # np.cov uses ddof=1 by default (bias=False)
    C = np.cov(x, y, bias=False)
    cov12 = float(C[0, 1])
    var1  = float(C[0, 0])
    var2  = float(C[1, 1])

    return cov12, var1, var2, float(mean1), float(mean2)


def standard_model_suppression(sel, alpha_site, a0=1.0):
    """
    Compute per-neuron standard model suppression for a paired condition.

    sel        : selectivity value (n1Sel or n2Sel)
                 positive → neuron prefers center (loc 0) → nonpref is periph
                 negative → neuron prefers periph (loc 1) → nonpref is center
    alpha_site : alpha value at the peripheral site (a_1..a_4 from stdNormFits)
    a0         : center alpha, always fixed = 1.0

    Returns the weight fraction of the nonpreferred location:
        if sel >= 0: nonpref is periph → alpha_site / (a0 + alpha_site)
        if sel <  0: nonpref is center → a0 / (a0 + alpha_site)
    """
    if (not np.isfinite(sel)) or (not np.isfinite(alpha_site)) or alpha_site <= 0:
        return np.nan
    if sel >= 0:
        return alpha_site / (a0 + alpha_site)
    else:
        return a0 / (a0 + alpha_site)


def std_model_selectivity(fit, site, center_isPref, periph_isPref):
    """
    Selectivity from raw L drives of the standard-model fit.
    fit layout: [Lp_0, Lnp_0, Lp_1, Lnp_1, ..., Lp_4, Lnp_4, a_1..a_4, sigma, b]
    site: 0..3 (peripheral site 1..4)
    Returns (L0 - L1) / (L0 + L1); positive = center drives more.
    """
    if fit is None:
        return np.nan
    L0 = fit[0] if center_isPref else fit[1]
    L1 = fit[2 + 2*site] if periph_isPref else fit[3 + 2*site]
    denom = L0 + L1
    if (not np.isfinite(denom)) or denom < 1e-12:
        return np.nan
    return float((L0 - L1) / denom)


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
sig_corr_lookup = build_sig_corr_lookup(similarity_records)

# Faster index maps for master arrays to avoid repeated np.where calls
good_index = {u: i for i, u in enumerate(masterGoodUnits)}
corr_index = {u: i for i, u in enumerate(masterCorrUnits)}


# Map unit string -> index in masterGoodUnits (already exists as good_index)
# stdNormFits[good_index[unit]] gives the 16-param fit for that unit
# alpha indices in stdNormFits: [10, 11, 12, 13] = a_1, a_2, a_3, a_4
STD_ALPHA_OFFSET = 10   # stdNormFits[i][10 + site] = alpha at that site

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
    sig_corr_val = sig_corr_lookup.get(key)

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
        block = "centerN_offsetsPN"
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

        # For paired blocks: whether CENTER stimulus is preferred
        paired_center_isPref = (block == "centerP_offsetsPN")
        # For singles: infer preference from i parity if that’s how your indices are defined
        # (If you’re not 100% sure, set these to None for now and we’ll derive later from index maps)
        single_center_isPref = (i % 2 == 1)  # <-- ONLY if your single_indices_for_pair follows this convention
        single_periph_isPref = (i % 2 == 1)  # adjust if needed

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

        if (any(v < 0 for v in (pairStimCorr, singleCorr_center, singleCorr_periph))
                and pairWeightGeomMean > 0.5
                and bc_val < 0.5):
            print(comb, pairWeightGeomMean, bc_val, pairStimCorr, singleCorr_center, singleCorr_periph)

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

        # ── Standard model suppression ──────────────────────────────────────────────
        # stdNormFits[n][10 + site] = alpha at peripheral site (a_1..a_4)
        # a0 (center) is always fixed = 1.0
        n1_std_fit = stdNormFits[n1] if (stdNormFits is not None and
                                         not np.any(np.isnan(stdNormFits[n1]))) else None
        n2_std_fit = stdNormFits[n2] if (stdNormFits is not None and
                                         not np.any(np.isnan(stdNormFits[n2]))) else None

        if n1_std_fit is not None and n2_std_fit is not None:
            alpha_site_n1 = float(n1_std_fit[STD_ALPHA_OFFSET + site])
            alpha_site_n2 = float(n2_std_fit[STD_ALPHA_OFFSET + site])

            n1_std_supp = standard_model_suppression(n1Sel, alpha_site_n1)
            n2_std_supp = standard_model_suppression(n2Sel, alpha_site_n2)

            if np.isfinite(n1_std_supp) and np.isfinite(n2_std_supp):
                pairStdSuppGeom = float(np.sqrt(n1_std_supp * n2_std_supp))
            else:
                pairStdSuppGeom = np.nan
        else:
            n1_std_supp = np.nan
            n2_std_supp = np.nan
            pairStdSuppGeom = np.nan

        # ---- Gate the four fields by R² ----
        deltaPD_keep = deltaPD if r2_ok else np.nan
        avg_pair_weight_2w_keep = avg_pair_weight_2w if r2_ok else np.nan
        bc_keep = bc_val if r2_ok else np.nan
        bd_keep = bd_val if r2_ok else np.nan

        pairVar1 = float(pairDCov[0][0])
        pairVar2 = float(pairDCov[1][1])

        singleCenterVar1 = float(singleCenterCov[0][0])
        singleCenterVar2 = float(singleCenterCov[1][1])

        singlePeriphVar1 = float(singlePeriCov[0][0])
        singlePeriphVar2 = float(singlePeriCov[1][1])

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
                   sig_corr=sig_corr_val,
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
                   pairMeanFR_single_periph_hz=pairMeanFR_single_periph_hz,
                   site=site,
                   center_single_idx=center_single_idx,
                   periph_single_idx=periph_single_idx,
                   single_center_isPref=single_center_isPref,
                   single_periph_isPref=single_periph_isPref,
                   paired_center_isPref=paired_center_isPref,
                   pairVar1=pairVar1,
                   pairVar2=pairVar2,
                   n1StdSupp=n1_std_supp,
                   n2StdSupp=n2_std_supp,
                   pairStdSuppGeom=pairStdSuppGeom,
                   singleCenterVar1=singleCenterVar1,
                   singleCenterVar2=singleCenterVar2,
                   singlePeriphVar1=singlePeriphVar1,
                   singlePeriphVar2=singlePeriphVar2)


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
        # For paired blocks: whether CENTER stimulus is preferred
        paired_center_isPref = (block == "centerP_offsetsPN")
        # For singles: infer preference from i parity if that’s how your indices are defined
        # (If you’re not 100% sure, set these to None for now and we’ll derive later from index maps)
        single_center_isPref = (i % 2 == 1)  # <-- ONLY if your single_indices_for_pair follows this convention
        single_periph_isPref = (i % 2 == 1)  # adjust if needed

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

        if (any(v < 0.0 for v in (pairStimCorr, singleCorr_center, singleCorr_periph))
                and pairWeightGeomMean > 0.5
                and bc_val < 0.5):
            print(comb, pairWeightGeomMean, bc_val, pairStimCorr, singleCorr_center, singleCorr_periph)

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

        # ── Standard model suppression ──────────────────────────────────────────────
        # stdNormFits[n][10 + site] = alpha at peripheral site (a_1..a_4)
        # a0 (center) is always fixed = 1.0
        n1_std_fit = stdNormFits[n1] if (stdNormFits is not None and
                                         not np.any(np.isnan(stdNormFits[n1]))) else None
        n2_std_fit = stdNormFits[n2] if (stdNormFits is not None and
                                         not np.any(np.isnan(stdNormFits[n2]))) else None

        if n1_std_fit is not None and n2_std_fit is not None:
            alpha_site_n1 = float(n1_std_fit[STD_ALPHA_OFFSET + site])
            alpha_site_n2 = float(n2_std_fit[STD_ALPHA_OFFSET + site])

            n1_std_supp = standard_model_suppression(n1Sel, alpha_site_n1)
            n2_std_supp = standard_model_suppression(n2Sel, alpha_site_n2)

            if np.isfinite(n1_std_supp) and np.isfinite(n2_std_supp):
                pairStdSuppGeom = float(np.sqrt(n1_std_supp * n2_std_supp))
            else:
                pairStdSuppGeom = np.nan
        else:
            n1_std_supp = np.nan
            n2_std_supp = np.nan
            pairStdSuppGeom = np.nan

        # ---- Gate the four fields by R² ----
        deltaPD_keep = deltaPD if r2_ok else np.nan
        avg_pair_weight_2w_keep = avg_pair_weight_2w if r2_ok else np.nan
        bc_keep = bc_val if r2_ok else np.nan
        bd_keep = bd_val if r2_ok else np.nan

        pairVar1 = float(pairDCov[0][0])
        pairVar2 = float(pairDCov[1][1])

        singleCenterVar1 = float(singleCenterCov[0][0])
        singleCenterVar2 = float(singleCenterCov[1][1])

        singlePeriphVar1 = float(singlePeriCov[0][0])
        singlePeriphVar2 = float(singlePeriCov[1][1])

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
                   sig_corr=sig_corr_val,
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
                   pairMeanFR_single_periph_hz=pairMeanFR_single_periph_hz,
                   site=site,
                   center_single_idx=center_single_idx,
                   periph_single_idx=periph_single_idx,
                   single_center_isPref=single_center_isPref,
                   single_periph_isPref=single_periph_isPref,
                   paired_center_isPref=paired_center_isPref,
                   pairVar1=pairVar1,
                   pairVar2=pairVar2,
                   n1StdSupp=n1_std_supp,
                   n2StdSupp=n2_std_supp,
                   pairStdSuppGeom=pairStdSuppGeom,
                   singleCenterVar1=singleCenterVar1,
                   singleCenterVar2=singleCenterVar2,
                   singlePeriphVar1=singlePeriphVar1,
                   singlePeriphVar2=singlePeriphVar2)


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
    "sigCorr": corr_df["sig_corr"].values,
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
    "pairVar1": corr_df["pairVar1"].values,
    "pairVar2": corr_df["pairVar2"].values,
    "pairDenom": corr_df["pairDenom"].values,
    "singleCenterVar1": corr_df["singleCenterVar1"].values,
    "singleCenterVar2": corr_df["singleCenterVar2"].values,
    "singleCenterDenom": corr_df["singleCenterDenom"].values,
    "singlePeriphVar1": corr_df["singlePeriphVar1"].values,
    "singlePeriphVar2": corr_df["singlePeriphVar2"].values,
    "singlePeriphDenom": corr_df["singlePeriphDenom"].values,
    "stdSuppGeom": corr_df["pairStdSuppGeom"].values,
    "n1StdSupp": corr_df["n1StdSupp"].values,
    "n2StdSupp": corr_df["n2StdSupp"].values,
})


############################## FIG 3 B C ##############################
# plot Von Mises Overlapping Tuning Curves with Signal Correlation
for fller in range(1):
    units = ['230627_71', '230627_114']
    locations = [0, 1, 2, 3, 4]
    base_dir = "/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN"
    eps = 1e-12

    # --- Infer session date and monkey from unit name ---
    sessionDate = units[0].split('_')[0]
    monkeyName = "Meetz" if sessionDate.startswith("23") else "Akshan"

    # --- Load dirData and SEM for this session ---
    sesh_dir = os.path.join(base_dir, monkeyName, sessionDate, "Direction Tuning")
    os.chdir(sesh_dir)
    dirData = np.load('unitsDirTuningMat.npy')  # already in spikes/s
    semData = np.load('unitsDirTuningSEM.npy')  # already in spikes/s
    os.chdir(base_dir)

    # --- Build unit -> row index mapping from masterAllUnits ---
    units_today = [u for u in masterAllUnits if u.startswith(sessionDate + "_")]
    unit_to_row = {u: i for i, u in enumerate(units_today)}

    # --- Angle grid ---
    n_dirs = dirData.shape[1]
    angles_deg = np.linspace(0, 360, n_dirs, endpoint=False)
    angles_rad = np.deg2rad(angles_deg)
    angles_rad_fine = np.linspace(0, 2 * np.pi, 1000)  # for smooth von Mises curve

    # --- Compute signal correlation on fitted values at original 6 angles ---
    fitted_vals = []
    for unit in units:
        idx = unit_to_row[unit]
        r = dirData[idx, :].astype(float)
        b, A, prefD, k = fit_von_mises_single(r, angles_rad)
        if np.isfinite(b) and np.isfinite(A) and np.isfinite(prefD) and np.isfinite(k):
            fitted_vals.append(von_mises_func(angles_rad, b, A, prefD, k))
        else:
            fitted_vals.append(r)  # fallback to raw if fit failed

    sig_corr = float(np.corrcoef(fitted_vals[0], fitted_vals[1])[0, 1])

    # --- Plot ---
    fig = plt.figure(figsize=(2.5, 2.5))

    ax1 = fig.add_subplot(2, 1, 1, projection='polar')
    for unit in units:
        idx = unit_to_row[unit]
        r = dirData[idx, :].astype(float)

        # Fit von Mises
        b, A, prefD, k = fit_von_mises_single(r, angles_rad)

        if np.isfinite(b) and np.isfinite(A) and np.isfinite(prefD) and np.isfinite(k):
            r2 = r2_score_safe(r, von_mises_func(angles_rad, b, A, prefD, k))
            vm_fit = von_mises_func(angles_rad_fine, b, A, prefD, k)
            ax1.plot(angles_rad_fine, vm_fit, linewidth=2,
                     label=f'Unit {unit}')
        else:
            ax1.plot([], [], label=f'Unit {unit} (fit failed)')

    ax1.set_theta_zero_location('W')
    ax1.set_theta_direction(1)
    ax1.set_title(f'Signal Corr = {sig_corr:.3f}', pad=15, fontsize=fontSize, **hfont)
    for label in ax1.get_xticklabels():
        label.set_fontname(hfont['fontname'])
        label.set_fontsize(fontSize)
    # ax1.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

    # 2nd unit comparison
    units = ['240611_37', '240611_42']
    base_dir = "/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN"
    eps = 1e-12

    # --- Infer session date and monkey from unit name ---
    sessionDate = units[0].split('_')[0]
    monkeyName = "Meetz" if sessionDate.startswith("23") else "Akshan"

    # --- Load dirData and SEM for this session ---
    sesh_dir = os.path.join(base_dir, monkeyName, sessionDate, "Direction Tuning")
    os.chdir(sesh_dir)
    dirData = np.load('unitsDirTuningMat.npy')  # already in spikes/s
    semData = np.load('unitsDirTuningSEM.npy')  # already in spikes/s
    os.chdir(base_dir)

    # --- Build unit -> row index mapping from masterAllUnits ---
    units_today = [u for u in masterAllUnits if u.startswith(sessionDate + "_")]
    unit_to_row = {u: i for i, u in enumerate(units_today)}

    # --- Angle grid ---
    n_dirs = dirData.shape[1]
    angles_deg = np.linspace(0, 360, n_dirs, endpoint=False)
    angles_rad = np.deg2rad(angles_deg)
    angles_rad_fine = np.linspace(0, 2 * np.pi, 1000)  # for smooth von Mises curve

    # --- Compute PMFs and BC ---
    # --- Compute signal correlation on fitted values at original 6 angles ---
    fitted_vals = []
    for unit in units:
        idx = unit_to_row[unit]
        r = dirData[idx, :].astype(float)
        b, A, prefD, k = fit_von_mises_single(r, angles_rad)
        if np.isfinite(b) and np.isfinite(A) and np.isfinite(prefD) and np.isfinite(k):
            fitted_vals.append(von_mises_func(angles_rad, b, A, prefD, k))
        else:
            fitted_vals.append(r)  # fallback to raw if fit failed

    sig_corr = float(np.corrcoef(fitted_vals[0], fitted_vals[1])[0, 1])

    ax2 = fig.add_subplot(2, 1, 2, projection='polar')
    for unit in units:
        idx = unit_to_row[unit]
        r = dirData[idx, :].astype(float)

        # Fit von Mises
        b, A, prefD, k = fit_von_mises_single(r, angles_rad)

        if np.isfinite(b) and np.isfinite(A) and np.isfinite(prefD) and np.isfinite(k):
            r2 = r2_score_safe(r, von_mises_func(angles_rad, b, A, prefD, k))
            vm_fit = von_mises_func(angles_rad_fine, b, A, prefD, k)
            ax2.plot(angles_rad_fine, vm_fit, linewidth=2,
                     label=f'Unit {unit}')
        else:
            ax2.plot([], [], label=f'Unit {unit} (fit failed)')

    ax2.set_theta_zero_location('W')
    ax2.set_theta_direction(1)
    ax2.set_title(f'Signal Corr = {sig_corr:.3f}', **hfont, fontsize=fontSize, pad=15)
    for label in ax2.get_xticklabels():
        label.set_fontname(hfont['fontname'])
        label.set_fontsize(fontSize)
    # ax2.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

    plt.show()

    plt.savefig(
        '/Users/chery/Documents/Grad School/Maunsell Lab/Manuscripts/Correlations/PythonMainFigures/exampleSigCorr.pdf',
        bbox_inches='tight', dpi=500)
    plt.close('all')


"""
correlation_analysis_pipeline.py
=================================
Single unified script covering:
  1. Heuristic weight bootstrap  (precision_bootstrap_v2)
  2. Standard suppression bootstrap (precision_bootstrap_stdsupp)
  3. Extended LME: 10-model comparison across original and bootstrapped predictors

Models tested (each run twice: original + bootstrapped predictors):
  M0:  BC only                                    (baseline)
  M1:  BC + weight
  M2:  BC + FR
  M3:  BC + stdSupp
  M4:  BC + weight + FR
  M5:  BC + weight + stdSupp
  M6:  BC + FR + stdSupp
  M7:  BC + weight + FR    (weight on top of FR)
  M8:  BC + weight + stdSupp (weight on top of stdSupp)
  M9:  BC + FR + stdSupp   (FR on top of stdSupp)
  M10: BC + weight + FR + stdSupp  (kitchen sink)

All models include quadratic BC term and BC interactions with each predictor.
Random structure: (1 + w_c | pair_id) with intercept-only fallback.

Execution order:
  corr_df  (built upstream in main loop, contains pairStdSuppGeom)
      ↓ run_precision_bootstrap_v2
  corr_df_subsampled  (+heuristic boot cols)
      ↓ run_stdsupp_bootstrap
  corr_df_full  (+stdsupp boot cols)
      ↓ run_lme_extended

Requires in scope from main pipeline:
  corr_df, masterGoodUnits, rfWeightNeuronFits, stdNormFits,
  prefNormalizedWithBase, nonprefNormalizedWithBase,
  pnNormalizedWithBase, npNormalizedWithBase,
  ppNormalizedWithBase, nnNormalizedWithBase,
  sponNormalizedToPrefCenter,
  contrast_centerSpon, contrast_peripherySpon,
  locationsSpon, stim_type_centerSpon, stim_type_peripherySpon,
  fullRFWeightedNormHeuristic, fullStandardNorm,
  shared_geom_share_rectified
"""
for filler in range(1):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from scipy import stats
    from scipy.stats import pearsonr, chi2
    import statsmodels.formula.api as smf
    import warnings
    warnings.filterwarnings("ignore")


    # ══════════════════════════════════════════════════════════════════════════════
    #  SHARED CONSTANTS
    # ══════════════════════════════════════════════════════════════════════════════

    VALID_COND_TYPES = [2, 3, 4, 5]
    COND_TYPE_NAMES  = {2: "pn", 3: "np", 4: "pp", 5: "nn"}
    STD_ALPHA_OFFSET = 10   # stdNormFits[i][10 + site] = a_{site+1}


    # ══════════════════════════════════════════════════════════════════════════════
    #  SECTION 1 — STIM SUBSET BUILDER (shared by both bootstraps)
    # ══════════════════════════════════════════════════════════════════════════════

    def build_stim_subsets(contrast_centerSpon, contrast_peripherySpon,
                           locationsSpon, stim_type_centerSpon,
                           stim_type_peripherySpon):
        """
        Split full stim vectors into per-condition subsets.
        Returns dict: cond_idx -> dict(cc, cp, loc, stc, stp) as float arrays.
        Only returns subsets for VALID_COND_TYPES (paired conditions with periph).
        Center-only conditions (0, 1) carry no W/alpha info and are excluded.
        """
        cc  = np.asarray(contrast_centerSpon,    dtype=float)
        cp  = np.asarray(contrast_peripherySpon, dtype=float)
        loc = np.asarray(locationsSpon,          dtype=float)
        stc = np.asarray(stim_type_centerSpon,   dtype=float)
        stp = np.asarray(stim_type_peripherySpon, dtype=float)

        has_periph = loc >= 0
        all_masks = [
            (~has_periph) & (stc == 1),              # 0: pref center only   — excluded
            (~has_periph) & (stc == 0),              # 1: nonpref center only — excluded
            has_periph & (stc == 1) & (stp == 0),   # 2: pn ✓
            has_periph & (stc == 0) & (stp == 1),   # 3: np ✓
            has_periph & (stc == 1) & (stp == 1),   # 4: pp ✓
            has_periph & (stc == 0) & (stp == 0),   # 5: nn ✓
        ]

        subsets = {}
        for ci in VALID_COND_TYPES:
            m = all_masks[ci]
            subsets[ci] = dict(cc=cc[m], cp=cp[m], loc=loc[m], stc=stc[m], stp=stp[m])
            print(f"  Cond {ci} ({COND_TYPE_NAMES[ci]}): {m.sum()} stim points")

        return subsets


    # ══════════════════════════════════════════════════════════════════════════════
    #  SECTION 2 — HEURISTIC WEIGHT BOOTSTRAP
    # ══════════════════════════════════════════════════════════════════════════════

    def _partial_r(y, x, controls):
        n = len(y)
        Z = np.column_stack([np.ones(n), controls])
        def resid(v):
            coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
            return v - Z @ coef
        r, p = stats.pearsonr(resid(y), resid(x))
        return dict(r=float(r), p=float(p), n=int(n))


    def _refit_W_one_cond(neuron_idx, cond_type_idx, all_resp_arrays,
                          stim_subsets, full_fit_params):
        """Refit W1-W4 holding Lp, Lnp, L0, b fixed."""
        sub  = stim_subsets[cond_type_idx]
        resp = all_resp_arrays[cond_type_idx][neuron_idx].astype(float)
        cc, cp, loc, stc, stp = sub["cc"], sub["cp"], sub["loc"], sub["stc"], sub["stp"]

        if len(resp) != len(cc) or len(resp) == 0:
            return None

        fp        = full_fit_params
        Lp_f      = float(fp[0])
        Lnp_f     = float(fp[1])
        L0_f      = float(fp[6])
        b_f       = float(fp[7])

        def model_W_only(data, W1, W2, W3, W4):
            cc_, cp_, loc_, stc_, stp_ = data
            return [fullRFWeightedNormHeuristic(
                        float(cc_[k]), float(cp_[k]), float(loc_[k]),
                        float(stc_[k]), float(stp_[k]),
                        Lp_f, Lnp_f, W1, W2, W3, W4, L0_f, b_f)
                    for k in range(len(cc_))]
        try:
            popt, _ = curve_fit(model_W_only, (cc, cp, loc, stc, stp), resp,
                                p0=[float(fp[2]), float(fp[3]),
                                    float(fp[4]), float(fp[5])],
                                bounds=([0.]*4, [np.inf]*4), maxfev=10000)
            return tuple(float(w) for w in popt)
        except Exception:
            return None


    def _recompute_pair_weight_geommean(W_n1, W_n2, site, sig1, sig2):
        w1, w2 = W_n1[site], W_n2[site]
        if w1 <= 0 or w2 <= 0:
            return np.nan
        return shared_geom_share_rectified(w1=w1, w2=w2, sigma1=sig1, sigma2=sig2)


    def _recompute_single_periph_weight(W_n1, W_n2, site, sig1, sig2):
        """singleWeightGeom_periph = sqrt((w/(w+sigma))_n1 * (w/(w+sigma))_n2)"""
        w1, w2 = W_n1[site], W_n2[site]
        if w1 <= 0 or w2 <= 0:
            return np.nan
        return float(np.sqrt((w1 / (w1 + sig1)) * (w2 / (w2 + sig2))))


    def run_precision_bootstrap_v2(corr_df, masterGoodUnits, rfWeightNeuronFits,
                                    all_resp_arrays, stim_subsets,
                                    n_boot=200, random_seed=42):
        """
        Bootstrap heuristic W1-W4 from single conditions.

        Adds to corr_df:
          pair_weight_geommean_1cond_{mean,std,median}
          singleWeightGeom_periph_1cond_{mean,std,median}

        Returns (boot_df, corr_df_subsampled)
        """
        rng        = np.random.default_rng(random_seed)
        good_index = {u: i for i, u in enumerate(masterGoodUnits)}

        required = ["corr", "pairMeanFR_paired_hz", "pair_weight_geommean",
                    "BC_dir_tuning", "unit_i", "unit_j", "site"]
        df_work = corr_df.dropna(subset=required).copy().reset_index(drop=True)
        df_work["Z_corr"] = np.arctanh(np.clip(df_work["corr"], -0.9999, 0.9999))

        weight_boot_mat        = np.full((len(df_work), n_boot), np.nan)
        single_periph_boot_mat = np.full((len(df_work), n_boot), np.nan)
        boot_records = []

        for boot in range(n_boot):
            cond_type = int(rng.choice(VALID_COND_TYPES))

            refitted_W = {}
            for u, idx in good_index.items():
                refitted_W[u] = _refit_W_one_cond(idx, cond_type, all_resp_arrays,
                                                   stim_subsets, rfWeightNeuronFits[idx])

            new_weights = np.full(len(df_work), np.nan)
            for row_num, (_, row) in enumerate(df_work.iterrows()):
                u1, u2   = row["unit_i"], row["unit_j"]
                site     = int(row["site"])
                W1, W2   = refitted_W.get(u1), refitted_W.get(u2)
                if W1 is None or W2 is None:
                    continue
                n1_idx   = good_index.get(u1)
                n2_idx   = good_index.get(u2)
                if n1_idx is None or n2_idx is None:
                    continue
                sig1 = float(rfWeightNeuronFits[n1_idx][6])
                sig2 = float(rfWeightNeuronFits[n2_idx][6])
                new_weights[row_num] = _recompute_pair_weight_geommean(W1, W2, site, sig1, sig2)
                single_periph_boot_mat[row_num, boot] = _recompute_single_periph_weight(
                    W1, W2, site, sig1, sig2)

            weight_boot_mat[:, boot] = new_weights

            # partial correlation stats
            valid = np.isfinite(new_weights)
            if valid.sum() < 50:
                continue
            y      = df_work.loc[valid, "Z_corr"].values
            rate   = df_work.loc[valid, "pairMeanFR_paired_hz"].values
            weight = new_weights[valid]
            bc     = df_work.loc[valid, "BC_dir_tuning"].values
            ok     = np.isfinite(rate) & np.isfinite(bc) & np.isfinite(weight)
            if ok.sum() < 50:
                continue

            res_rate   = _partial_r(y[ok], rate[ok],   np.column_stack([weight[ok], bc[ok]]))
            res_weight = _partial_r(y[ok], weight[ok], np.column_stack([rate[ok],   bc[ok]]))
            boot_records.append(dict(
                boot=boot, cond_type=cond_type, cond_name=COND_TYPE_NAMES[cond_type],
                n=ok.sum(),
                partial_r_rate=res_rate["r"],   p_rate=res_rate["p"],
                partial_r_weight=res_weight["r"], p_weight=res_weight["p"],
                diff_weight_minus_rate=res_weight["r"] - res_rate["r"],
            ))

            if (boot + 1) % 25 == 0:
                print(f"  [W boot {boot+1:>3}/{n_boot}] cond={COND_TYPE_NAMES[cond_type]}  "
                      f"r_rate={res_rate['r']:+.4f}  r_weight={res_weight['r']:+.4f}  n={ok.sum()}")

        # summarise
        for col, arr in [
            ("pair_weight_geommean_1cond_mean",      np.nanmean(weight_boot_mat,        axis=1)),
            ("pair_weight_geommean_1cond_std",       np.nanstd(weight_boot_mat,         axis=1)),
            ("pair_weight_geommean_1cond_median",    np.nanmedian(weight_boot_mat,      axis=1)),
            ("singleWeightGeom_periph_1cond_mean",   np.nanmean(single_periph_boot_mat, axis=1)),
            ("singleWeightGeom_periph_1cond_std",    np.nanstd(single_periph_boot_mat,  axis=1)),
            ("singleWeightGeom_periph_1cond_median", np.nanmedian(single_periph_boot_mat, axis=1)),
        ]:
            df_work[col] = arr

        new_cols = [c for c in df_work.columns if "1cond" in c]
        corr_df_out = corr_df.copy().merge(
            df_work[["unit_i", "unit_j", "block", "cond_index"] + new_cols],
            on=["unit_i", "unit_j", "block", "cond_index"], how="left"
        )

        print(f"\n  [W bootstrap done] "
              f"{corr_df_out['pair_weight_geommean_1cond_mean'].notna().sum():,} rows populated")
        return pd.DataFrame(boot_records), corr_df_out


    # ══════════════════════════════════════════════════════════════════════════════
    #  SECTION 3 — STANDARD SUPPRESSION BOOTSTRAP
    # ══════════════════════════════════════════════════════════════════════════════

    def _refit_alphas_one_cond(neuron_idx, cond_type_idx, all_resp_arrays,
                                stim_subsets, full_std_params):
        """Refit a_1..a_4 holding all other stdNormFits params fixed."""
        sub  = stim_subsets[cond_type_idx]
        resp = all_resp_arrays[cond_type_idx][neuron_idx].astype(float)
        cc, cp, loc, stc, stp = sub["cc"], sub["cp"], sub["loc"], sub["stc"], sub["stp"]

        if len(resp) != len(cc) or len(resp) == 0:
            return None

        fp = full_std_params
        Lp_0_f,  Lnp_0_f  = float(fp[0]),  float(fp[1])
        Lp_1_f,  Lnp_1_f  = float(fp[2]),  float(fp[3])
        Lp_2_f,  Lnp_2_f  = float(fp[4]),  float(fp[5])
        Lp_3_f,  Lnp_3_f  = float(fp[6]),  float(fp[7])
        Lp_4_f,  Lnp_4_f  = float(fp[8]),  float(fp[9])
        sigma_f            = float(fp[14])
        b_f                = float(fp[15])

        def model_alpha_only(data, a_1, a_2, a_3, a_4):
            cc_, cp_, loc_, stc_, stp_ = data
            return [fullStandardNorm(
                        float(cc_[k]), float(cp_[k]), float(loc_[k]),
                        float(stc_[k]), float(stp_[k]),
                        Lp_0_f, Lnp_0_f, Lp_1_f, Lnp_1_f, Lp_2_f, Lnp_2_f,
                        Lp_3_f, Lnp_3_f, Lp_4_f, Lnp_4_f,
                        a_1, a_2, a_3, a_4, sigma_f, b_f)
                    for k in range(len(cc_))]
        try:
            popt, _ = curve_fit(model_alpha_only, (cc, cp, loc, stc, stp), resp,
                                p0=[float(fp[10]), float(fp[11]),
                                    float(fp[12]), float(fp[13])],
                                bounds=([0.]*4, [np.inf]*4), maxfev=10000)
            return tuple(float(a) for a in popt)
        except Exception:
            return None


    def _recompute_pair_std_supp(n1Sel, n2Sel, alpha_site_n1, alpha_site_n2, a0=1.0):
        def single_supp(sel, alpha_site):
            if alpha_site <= 0 or not np.isfinite(alpha_site):
                return np.nan
            return alpha_site / (a0 + alpha_site) if sel >= 0 else a0 / (a0 + alpha_site)
        s1 = single_supp(n1Sel, alpha_site_n1)
        s2 = single_supp(n2Sel, alpha_site_n2)
        return float(np.sqrt(s1 * s2)) if np.isfinite(s1) and np.isfinite(s2) else np.nan


    def run_stdsupp_bootstrap(corr_df, masterGoodUnits, stdNormFits,
                               all_resp_arrays, stim_subsets,
                               n_boot=200, random_seed=42):
        """
        Bootstrap standard model a_1..a_4 from single conditions.

        Adds to corr_df:
          pairStdSuppGeom_1cond_{mean,std,median}

        Returns corr_df_full
        """
        rng        = np.random.default_rng(random_seed)
        good_index = {u: i for i, u in enumerate(masterGoodUnits)}

        def _std_valid(idx):
            fp = stdNormFits[idx]
            return fp is not None and np.all(np.isfinite(fp))

        required = ["unit_i", "unit_j", "site", "n1Sel", "n2Sel",
                    "corr", "pairMeanFR_paired_hz", "pair_weight_geommean",
                    "BC_dir_tuning"]
        df_work  = corr_df.dropna(subset=required).copy().reset_index(drop=True)

        supp_boot_mat = np.full((len(df_work), n_boot), np.nan)

        for boot in range(n_boot):
            cond_type = int(rng.choice(VALID_COND_TYPES))

            refitted_alphas = {}
            for u, idx in good_index.items():
                if not _std_valid(idx):
                    refitted_alphas[u] = None
                    continue
                refitted_alphas[u] = _refit_alphas_one_cond(
                    idx, cond_type, all_resp_arrays, stim_subsets, stdNormFits[idx])

            new_supp = np.full(len(df_work), np.nan)
            for row_num, row in df_work.iterrows():
                A1 = refitted_alphas.get(row["unit_i"])
                A2 = refitted_alphas.get(row["unit_j"])
                if A1 is None or A2 is None:
                    continue
                site = int(row["site"])
                new_supp[row_num] = _recompute_pair_std_supp(
                    row["n1Sel"], row["n2Sel"], A1[site], A2[site])

            supp_boot_mat[:, boot] = new_supp

            if (boot + 1) % 25 == 0:
                print(f"  [S boot {boot+1:>3}/{n_boot}] cond={COND_TYPE_NAMES[cond_type]}  "
                      f"valid_rows={np.isfinite(new_supp).sum()}")

        for col, arr in [
            ("pairStdSuppGeom_1cond_mean",   np.nanmean(supp_boot_mat,   axis=1)),
            ("pairStdSuppGeom_1cond_std",    np.nanstd(supp_boot_mat,    axis=1)),
            ("pairStdSuppGeom_1cond_median", np.nanmedian(supp_boot_mat, axis=1)),
        ]:
            df_work[col] = arr

        new_cols = ["pairStdSuppGeom_1cond_mean",
                    "pairStdSuppGeom_1cond_std",
                    "pairStdSuppGeom_1cond_median"]
        corr_df_out = corr_df.copy().merge(
            df_work[["unit_i", "unit_j", "block", "cond_index"] + new_cols],
            on=["unit_i", "unit_j", "block", "cond_index"], how="left"
        )

        print(f"\n  [S bootstrap done] "
              f"{corr_df_out['pairStdSuppGeom_1cond_mean'].notna().sum():,} rows populated")
        return corr_df_out


    # ══════════════════════════════════════════════════════════════════════════════
    #  SECTION 4 — LME HELPERS
    # ══════════════════════════════════════════════════════════════════════════════

    def _pseudo_r2(res):
        var_fe  = np.var(res.fittedvalues)
        var_res = res.scale
        var_re  = float(np.diag(res.cov_re).sum())
        r2_m = var_fe / (var_fe + var_re + var_res)
        r2_c = (var_fe + var_re) / (var_fe + var_re + var_res)
        return r2_m, r2_c


    def _fixed_pred(res, df):
        fe   = res.fe_params
        pred = pd.Series(0.0, index=df.index)
        for term, coef in fe.items():
            if term == "Intercept":
                pred += coef
            elif term in df.columns:
                pred += coef * df[term]
        return pred


    def _fit_one_lme(formula, df, label, re_formula="~ w_c", verbose=True):
        """Fit one LME with fallback to intercept-only RE."""
        for re_f in [re_formula, "~ 1"]:
            try:
                res = smf.mixedlm(formula, df, groups=df["pair_id"],
                                  re_formula=re_f).fit(reml=False, method="lbfgs",
                                                       maxiter=400)
                break
            except Exception as e:
                if re_f == "~ 1":
                    print(f"  WARNING: {label} — both RE structures failed: {e}")
                    return None

        r2_m, r2_c = _pseudo_r2(res)
        y_fe        = _fixed_pred(res, df)
        r_fe,  _    = pearsonr(df["rsc"], y_fe)
        r_full, _   = pearsonr(df["rsc"], res.fittedvalues)

        if verbose:
            print(f"\n{'═'*65}\n  {label}")
            print(f"  AIC={res.aic:.2f}  BIC={res.bic:.2f}  "
                  f"R²_marg={r2_m:.4f}  R²_cond={r2_c:.4f}")
            print(f"  Fixed-only r={r_fe:.4f}  Full r={r_full:.4f}")
            print(f"  {'term':<28} {'coef':>9} {'SE':>9} {'z':>7} {'p':>10}")
            for term in res.fe_params.index:
                z   = res.fe_params[term] / res.bse_fe[term]
                p   = res.pvalues[term]
                sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
                print(f"  {term:<28} {res.fe_params[term]:>9.4f} "
                      f"{res.bse_fe[term]:>9.4f} {z:>7.2f} {p:>10.4f} {sig}")

        return dict(res=res, r2m=r2_m, r2c=r2_c, y_fe=y_fe,
                    y_full=res.fittedvalues, r_fe=r_fe, r_full=r_full, label=label)


    def _run_lrts(fits, lrt_pairs, label=""):
        print(f"\n  LRTs [{label}]:")
        rows = []
        for null_k, alt_k, desc in lrt_pairs:
            if any(k not in fits or fits[k] is None for k in (null_k, alt_k)):
                continue
            null, alt = fits[null_k]["res"], fits[alt_k]["res"]
            df_d = len(alt.fe_params) - len(null.fe_params)
            if df_d <= 0:
                continue
            stat = 2 * (alt.llf - null.llf)
            p    = chi2.sf(stat, df_d)
            sig  = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
            print(f"    {desc:<42} χ²({df_d})={stat:7.2f}  p={p:.4e} {sig}")
            rows.append(dict(comparison=desc, chi2=stat, df=df_d, p=p, sig=sig))
        return pd.DataFrame(rows)


    def compute_standardized_betas(fits, df_stacked, round_label=""):
        """
        Compute standardized regression coefficients for w_c, fr_c, s_c
        across all models that contain each predictor.

        β_std = coef * (SD_predictor / SD_rsc)

        A 1 SD increase in predictor X corresponds to β_std SD change in rsc.
        """
        sd_rsc = df_stacked["rsc"].std()
        sd_w = df_stacked["w_c"].std()
        sd_fr = df_stacked["fr_c"].std()
        sd_s = df_stacked["s_c"].std()
        sd_sc = df_stacked["sc_c"].std() if "sc_c" in df_stacked.columns else np.nan  # <-- ADD
        sd_map = {"w_c": sd_w, "fr_c": sd_fr, "s_c": sd_s, "sc_c": sd_sc}  # <-- ADD sc_c

        print(f"\n{'═' * 70}")
        print(f"  STANDARDIZED BETAS [{round_label}]")
        print(f"  SD(rsc)={sd_rsc:.4f}  SD(w)={sd_w:.4f}  "
              f"SD(fr)={sd_fr:.4f}  SD(s)={sd_s:.4f}  SD(sc)={sd_sc:.4f}")  # <-- ADD SD(sc)
        print(f"{'═' * 70}")
        print(f"  {'model':<35} {'predictor':<10} {'coef':>9} "
              f"{'SD_pred':>9} {'β_std':>9} {'p':>10}")
        print(f"  {'-' * 68}")

        rows = []
        for key, _, label in _MODEL_SPECS:
            fit = fits.get(key)
            if fit is None:
                continue
            fe = fit["res"].fe_params
            pv = fit["res"].pvalues
            for cname in ["w_c", "fr_c", "s_c", "sc_c"]:  # <-- ADD sc_c
                if cname not in fe.index:
                    continue
                coef = float(fe[cname])
                sd_p = sd_map[cname]
                b_std = coef * sd_p / sd_rsc
                p = float(pv[cname])
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
                print(f"  {label:<35} {cname:<10} {coef:>9.4f} "
                      f"{sd_p:>9.4f} {b_std:>9.4f} {p:>10.4f} {sig}")
                rows.append(dict(model=label, predictor=cname,
                                 coef=coef, sd_pred=sd_p,
                                 beta_std=b_std, p=p))

        print(f"{'═' * 70}")
        return pd.DataFrame(rows)


    def compute_delta_r2(fits, round_label=""):
        """
        ΔR²_fe = r²_fe(larger model) - r²_fe(smaller model)
        for each LRT pair. Expresses how much additional variance
        in rsc the added predictor(s) explain at the population level.
        """
        # Detect which baseline key is present
        m0_key = "M0_bc" if "M0_bc" in fits else "M0_sc"

        lrt_pairs_for_delta = [
            (m0_key, "M1_weight", "weight alone"),
            (m0_key, "M2_fr", "FR alone"),
            (m0_key, "M3_supp", "stdSupp alone"),
            ("M2_fr", "M4_wt_fr", "weight | FR"),
            ("M3_supp", "M5_wt_supp", "weight | stdSupp"),
            ("M1_weight", "M4_wt_fr", "FR | weight"),
            ("M3_supp", "M6_fr_supp", "FR | stdSupp"),
            ("M1_weight", "M5_wt_supp", "stdSupp | weight"),
            ("M2_fr", "M6_fr_supp", "stdSupp | FR"),
            ("M6_fr_supp", "M7_all", "weight | FR+stdSupp"),
            ("M5_wt_supp", "M7_all", "FR | weight+stdSupp"),
            ("M4_wt_fr", "M7_all", "stdSupp | weight+FR"),
        ]

        print(f"\n{'═' * 65}")
        print(f"  ΔR²_fe (incremental fixed-effects variance) [{round_label}]")
        print(f"{'═' * 65}")
        print(f"  {'comparison':<35} {'R²_null':>8} {'R²_alt':>8} "
              f"{'ΔR²':>8} {'Δr_fe':>8}")
        print(f"  {'-' * 63}")

        rows = []
        for null_k, alt_k, desc in lrt_pairs_for_delta:
            null_fit = fits.get(null_k)
            alt_fit = fits.get(alt_k)
            if null_fit is None or alt_fit is None:
                continue

            r2_null = null_fit["r_fe"] ** 2
            r2_alt = alt_fit["r_fe"] ** 2
            delta_r2 = r2_alt - r2_null
            # also report the increment in correlation coefficient itself
            delta_r = alt_fit["r_fe"] - null_fit["r_fe"]

            print(f"  {desc:<35} {r2_null:>8.4f} {r2_alt:>8.4f} "
                  f"{delta_r2:>+8.4f} {delta_r:>+8.4f}")
            rows.append(dict(comparison=desc,
                             r2_null=r2_null, r2_alt=r2_alt,
                             delta_r2=delta_r2, delta_r_fe=delta_r))

        print(f"{'═' * 65}")
        return pd.DataFrame(rows)


    # ══════════════════════════════════════════════════════════════════════════════
    #  SECTION 5 — LME: 10 models, original vs bootstrapped predictors
    # ══════════════════════════════════════════════════════════════════════════════
    BC_COL    = "BC_dir_tuning"
    SIG_CORR_COL = "sig_corr"
    UNIT_I    = "unit_i"
    UNIT_J    = "unit_j"
    PAIRED_CORR   = "corr"
    SINGLE_C_CORR = "singleCorr_center"
    SINGLE_P_CORR = "singleCorr_periph"
    PAIRED_FR     = "pairMeanFR_paired_hz"
    SINGLE_C_FR   = "pairMeanFR_single_center_hz"
    SINGLE_P_FR   = "pairMeanFR_single_periph_hz"
    SC_WEIGHT     = "singleWeightGeom_center"

    _MODEL_SPECS = [
        # ── Single-predictor models ────────────────────────────────────────────────
        ("M0_bc",
         "rsc ~ bc_c + bc2",
         "M0: BC only"),

        ("M1_weight",
         "rsc ~ bc_c + bc2 + w_c + bc_w + w2",
         "M1: BC + weight"),

        ("M2_fr",
         "rsc ~ bc_c + bc2 + fr_c + bc_fr + fr2",
         "M2: BC + FR"),

        ("M3_supp",
         "rsc ~ bc_c + bc2 + s_c + bc_s + s2",
         "M3: BC + stdSupp"),

        # ── Two-predictor models ───────────────────────────────────────────────────
        ("M4_wt_fr",
         "rsc ~ bc_c + bc2 + w_c + fr_c + bc_w + bc_fr + w2 + fr2",
         "M4: BC + weight + FR"),

        ("M5_wt_supp",
         "rsc ~ bc_c + bc2 + w_c + s_c + bc_w + bc_s + w2 + s2",
         "M5: BC + weight + stdSupp"),

        ("M6_fr_supp",
         "rsc ~ bc_c + bc2 + fr_c + s_c + bc_fr + bc_s + fr2 + s2",
         "M6: BC + FR + stdSupp"),

        # ── Kitchen sink ──────────────────────────────────────────────────────────
        ("M7_all",
         "rsc ~ bc_c + bc2 + w_c + fr_c + s_c + bc_w + bc_fr + bc_s + w2 + fr2 + s2",
         "M7: BC + weight + FR + stdSupp"),
    ]

    _MODEL_SPECS_SC = [
        ("M0_sc",
         "rsc ~ sc_c + sc2",
         "M0: sigCorr only"),

        ("M1_weight",
         "rsc ~ sc_c + sc2 + w_c + sc_w + w2",
         "M1: sigCorr + weight"),

        ("M2_fr",
         "rsc ~ sc_c + sc2 + fr_c + sc_fr + fr2",
         "M2: sigCorr + FR"),

        ("M3_supp",
         "rsc ~ sc_c + sc2 + s_c + sc_s + s2",
         "M3: sigCorr + stdSupp"),

        ("M4_wt_fr",
         "rsc ~ sc_c + sc2 + w_c + fr_c + sc_w + sc_fr + w2 + fr2",
         "M4: sigCorr + weight + FR"),

        ("M5_wt_supp",
         "rsc ~ sc_c + sc2 + w_c + s_c + sc_w + sc_s + w2 + s2",
         "M5: sigCorr + weight + stdSupp"),

        ("M6_fr_supp",
         "rsc ~ sc_c + sc2 + fr_c + s_c + sc_fr + sc_s + fr2 + s2",
         "M6: sigCorr + FR + stdSupp"),

        ("M7_all",
         "rsc ~ sc_c + sc2 + w_c + fr_c + s_c + sc_w + sc_fr + sc_s + w2 + fr2 + s2",
         "M7: sigCorr + weight + FR + stdSupp"),
    ]

    _LRT_PAIRS = [
        # ── Each predictor alone vs baseline ──────────────────────────────────────
        ("M0_bc",    "M1_weight",  "M0→M1  (weight alone)"),
        ("M0_bc",    "M2_fr",      "M0→M2  (FR alone)"),
        ("M0_bc",    "M3_supp",    "M0→M3  (stdSupp alone)"),

        # ── Weight on top of each other predictor ─────────────────────────────────
        ("M2_fr",    "M4_wt_fr",   "M2→M4  (weight | FR)"),
        ("M3_supp",  "M5_wt_supp", "M3→M5  (weight | stdSupp)"),

        # ── FR on top of each other predictor ─────────────────────────────────────
        ("M1_weight","M4_wt_fr",   "M1→M4  (FR | weight)"),
        ("M3_supp",  "M6_fr_supp", "M3→M6  (FR | stdSupp)"),

        # ── stdSupp on top of each other predictor ────────────────────────────────
        ("M1_weight","M5_wt_supp", "M1→M5  (stdSupp | weight)"),
        ("M2_fr",    "M6_fr_supp", "M2→M6  (stdSupp | FR)"),

        # ── Each predictor on top of the other two ────────────────────────────────
        ("M6_fr_supp",  "M7_all",  "M6→M7  (weight | FR+stdSupp)"),
        ("M5_wt_supp",  "M7_all",  "M5→M7  (FR | weight+stdSupp)"),
        ("M4_wt_fr",    "M7_all",  "M4→M7  (stdSupp | weight+FR)"),
    ]

    _LRT_PAIRS_SC = [
        ("M0_sc", "M1_weight", "M0→M1  (weight alone)"),
        ("M0_sc", "M2_fr", "M0→M2  (FR alone)"),
        ("M0_sc", "M3_supp", "M0→M3  (stdSupp alone)"),
        ("M2_fr", "M4_wt_fr", "M2→M4  (weight | FR)"),
        ("M3_supp", "M5_wt_supp", "M3→M5  (weight | stdSupp)"),
        ("M1_weight", "M4_wt_fr", "M1→M4  (FR | weight)"),
        ("M3_supp", "M6_fr_supp", "M3→M6  (FR | stdSupp)"),
        ("M1_weight", "M5_wt_supp", "M1→M5  (stdSupp | weight)"),
        ("M2_fr", "M6_fr_supp", "M2→M6  (stdSupp | FR)"),
        ("M6_fr_supp", "M7_all", "M6→M7  (weight | FR+stdSupp)"),
        ("M5_wt_supp", "M7_all", "M5→M7  (FR | weight+stdSupp)"),
        ("M4_wt_fr", "M7_all", "M4→M7  (stdSupp | weight+FR)"),
    ]

    _MODEL_COLORS = {
        "M0_bc":       "#7f7f7f",
        "M0_sc":       "#7f7f7f",
        "M1_weight":   "#2ca02c",
        "M2_fr":       "#ff7f0e",
        "M3_supp":     "#d62728",
        "M4_wt_fr":    "#1f77b4",
        "M5_wt_supp":  "#9467bd",
        "M6_fr_supp":  "#e377c2",
        "M7_all":      "#8c564b",
    }


    def _build_stacked(corr_df, weight_col, supp_col, sp_weight_col, label=""):
        """
        Stack paired + single_center + single_periph.
        Weight varies by stim_type; stdSupp is structural (broadcast across all 3 rows).
        """
        # Step 1: supp lookup — one value per (pair × block × cond_index)
        supp_lookup = (corr_df[["unit_i","unit_j","block","cond_index",supp_col]]
                       .dropna(subset=[supp_col])
                       .drop_duplicates(subset=["unit_i","unit_j","block","cond_index"])
                       .copy())
        supp_lookup["pair_id"] = supp_lookup.apply(
            lambda r: "_".join(sorted([str(r["unit_i"]), str(r["unit_j"])])), axis=1)
        supp_lookup = supp_lookup[["pair_id","block","cond_index",supp_col]]

        # Step 2: base df
        needed = list(dict.fromkeys([BC_COL, SIG_CORR_COL, UNIT_I, UNIT_J, "block", "cond_index", "site",
                   PAIRED_CORR, SINGLE_C_CORR, SINGLE_P_CORR,
                   PAIRED_FR, SINGLE_C_FR, SINGLE_P_FR,
                   weight_col, SC_WEIGHT, sp_weight_col]))
        df = corr_df[needed].copy().dropna(subset=needed)
        df["pair_id"] = df.apply(
            lambda r: "_".join(sorted([str(r[UNIT_I]), str(r[UNIT_J])])), axis=1)

        # Step 3: stack stim types (weight varies; supp will be broadcast)
        # ── Paired rows: one row per cond_index, no deduplication needed ─────────
        frames = []
        tmp = df[["pair_id", "block", "cond_index", BC_COL, SIG_CORR_COL,
                  PAIRED_CORR, weight_col, PAIRED_FR]].copy()
        tmp.columns = ["pair_id", "block", "cond_index", "bc", 'sc', "rsc", "weight", "fr"]
        tmp["stim_type"] = "paired"
        frames.append(tmp)

        # ── single_center: deduplicate to 1 row per (pair, block) ────────────────
        # Keep one cond_index per unique (pair, block, singleCorr_center) —
        # the supp value from that cond_index will be broadcast via supp_lookup.
        # We keep the first occurrence which preserves a valid cond_index for the merge.
        sc_cols = ["pair_id", "block", "cond_index", BC_COL, SIG_CORR_COL,
                   SINGLE_C_CORR, SC_WEIGHT, SINGLE_C_FR]
        sc = (df[sc_cols]
              .drop_duplicates(subset=["pair_id","block",SINGLE_C_CORR])
              .copy())
        sc.columns = ["pair_id", "block", "cond_index", "bc", "sc", "rsc", "weight", "fr"]
        sc["stim_type"] = "single_center"
        frames.append(sc[["pair_id", "block", "cond_index",
                           "bc", "sc", "rsc", "weight", "fr", "stim_type"]])

        # ── single_periph: deduplicate to 1 row per (pair, block, unique rsc) ─────
        # Each unique periph measurement appears twice (once per center block type).
        # Keep first occurrence — preserves the cond_index so supp_lookup works.
        sp_cols = ["pair_id", "block", "cond_index", "site", BC_COL, SIG_CORR_COL,
                   SINGLE_P_CORR, sp_weight_col, SINGLE_P_FR]
        sp = (df[sp_cols]
              .drop_duplicates(subset=["pair_id", "site", SINGLE_P_CORR])
              .copy())
        sp.columns = ["pair_id", "block", "cond_index", "site", "bc", "sc", "rsc", "weight", "fr"]
        sp["stim_type"] = "single_periph"
        frames.append(sp[["pair_id", "block", "cond_index",
                           "bc", "sc", "rsc", "weight", "fr", "stim_type"]])

        stacked = pd.concat(frames, ignore_index=True)

        # Step 4: broadcast supp to all 3 rows per (pair × cond)
        stacked = (stacked
                   .merge(supp_lookup, on=["pair_id","block","cond_index"], how="left")
                   .rename(columns={supp_col: "supp"})
                   .dropna(subset=["bc","rsc","weight","supp","fr"]))

        # Step 5: center + quadratics + pairwise interactions with BC
        # Step 5: center all predictors first, then quadratics and interactions
        bc_mean = stacked["bc"].mean()
        w_mean = stacked["weight"].mean()
        s_mean = stacked["supp"].mean()
        fr_mean = stacked["fr"].mean()
        sc_mean = stacked["sc"].mean()

        stacked["bc_c"] = stacked["bc"] - bc_mean
        stacked["w_c"] = stacked["weight"] - w_mean
        stacked["s_c"] = stacked["supp"] - s_mean
        stacked["fr_c"] = stacked["fr"] - fr_mean
        stacked["sc_c"] = stacked["sc"] - sc_mean

        stacked["bc2"] = stacked["bc_c"] ** 2
        stacked["w2"] = stacked["w_c"] ** 2
        stacked["s2"] = stacked["s_c"] ** 2
        stacked["fr2"] = stacked["fr_c"] ** 2
        stacked["sc2"] = stacked["sc_c"] ** 2

        stacked["bc_w"] = stacked["bc_c"] * stacked["w_c"]
        stacked["bc_s"] = stacked["bc_c"] * stacked["s_c"]
        stacked["bc_fr"] = stacked["bc_c"] * stacked["fr_c"]
        stacked["sc_w"] = stacked["sc_c"] * stacked["w_c"]
        stacked["sc_s"] = stacked["sc_c"] * stacked["s_c"]
        stacked["sc_fr"] = stacked["sc_c"] * stacked["fr_c"]

        stacked["stim_type"] = pd.Categorical(
        stacked["stim_type"],
        categories=["paired", "single_center", "single_periph"])

        counts  = stacked["stim_type"].value_counts()
        n_pairs = stacked["pair_id"].nunique()
        print(f"  [{label}] rows={len(stacked):,}  pairs={n_pairs}  "
              f"paired={counts.get('paired',0):,}  "
              f"sc={counts.get('single_center',0):,}  "
              f"sp={counts.get('single_periph',0):,}")
        print(f"    mean_w={w_mean:.3f}  mean_supp={s_mean:.3f}  "
              f"mean_bc={bc_mean:.3f} mean_sc={sc_mean:.3f} mean_fr={fr_mean:.2f}Hz")
        return stacked


    def run_lme_extended(corr_df_full):
        """
        10-model LME: weight + FR + stdSupp and all combinations.
        Runs twice — original predictors and bootstrapped (precision-equated) predictors.
        Returns (all_fits, comp)
        """
        def _one_round(weight_col, supp_col, sp_weight_col, label):
            print(f"\n{'█'*72}\n  ROUND: {label}\n{'█'*72}")
            df   = _build_stacked(corr_df_full, weight_col, supp_col,
                                   sp_weight_col, label=label)
            fits = {k: _fit_one_lme(f, df, f"[{label}] {l}")
                    for k, f, l in _MODEL_SPECS}
            _run_lrts(fits, _LRT_PAIRS, label=label)
            compute_standardized_betas(fits, df, round_label=label)
            compute_delta_r2(fits, round_label=label)
            return fits, df

        fits_orig, df_orig = _one_round(
            "pair_weight_geommean",           "pairStdSuppGeom",
            "singleWeightGeom_periph",        "Original (all conditions)")
        fits_boot, df_boot = _one_round(
            "pair_weight_geommean_1cond_mean", "pairStdSuppGeom_1cond_mean",
            "singleWeightGeom_periph_1cond_mean", "Precision-equated (bootstrapped mean)")

        # ── Comparison table ──────────────────────────────────────────────────────
        rows = []
        for key, _, label in _MODEL_SPECS:
            o, b = fits_orig.get(key), fits_boot.get(key)
            if o is None or b is None:
                continue
            row = dict(model=label,
                       AIC_orig=o["res"].aic,  AIC_boot=b["res"].aic,
                       R2m_orig=o["r2m"],      R2m_boot=b["r2m"],
                       rfe_orig=o["r_fe"],     rfe_boot=b["r_fe"])
            for cname in ["w_c", "s_c", "fr_c"]:
                for suf, fits in [("_orig", fits_orig), ("_boot", fits_boot)]:
                    fe = fits[key]["res"].fe_params
                    row[f"{cname}{suf}"] = float(fe[cname]) if cname in fe.index else np.nan
            rows.append(row)
        comp = pd.DataFrame(rows)

        print(f"\n{'═'*100}\n  LME COMPARISON: Original vs Precision-Equated\n{'═'*100}")
        pd.set_option("display.float_format", "{:+.4f}".format)
        print(comp[["model", "AIC_orig", "AIC_boot",
                    "R2m_orig", "R2m_boot",
                    "w_c_orig", "w_c_boot",
                    "s_c_orig", "s_c_boot",
                    "fr_c_orig","fr_c_boot"]].to_string(index=False))
        pd.reset_option("display.float_format")

        _plot_extended_comparison(fits_orig, fits_boot)

        return {"orig": fits_orig, "boot": fits_boot}, comp, df_orig, df_boot


    def run_lme_sigcorr(corr_df_full):
        """
        Same 8-model LME as run_lme_extended but uses signal correlation
        as the base tuning-similarity predictor instead of BC.
        """
        import copy

        # Temporarily swap model specs
        global _MODEL_SPECS, _LRT_PAIRS
        _orig_specs = _MODEL_SPECS
        _orig_lrts = _LRT_PAIRS
        _MODEL_SPECS = _MODEL_SPECS_SC
        _LRT_PAIRS = _LRT_PAIRS_SC

        try:
            def _one_round(weight_col, supp_col, sp_weight_col, label):
                print(f"\n{'█' * 72}\n  ROUND: {label}\n{'█' * 72}")
                df = _build_stacked(corr_df_full, weight_col, supp_col,
                                    sp_weight_col, label=label)
                fits = {k: _fit_one_lme(f, df, f"[{label}] {l}")
                        for k, f, l in _MODEL_SPECS}
                _run_lrts(fits, _LRT_PAIRS, label=label)
                compute_standardized_betas(fits, df, round_label=label)
                compute_delta_r2(fits, round_label=label)
                return fits, df

            fits_orig, df_orig = _one_round(
                "pair_weight_geommean", "pairStdSuppGeom",
                "singleWeightGeom_periph", "sigCorr base — Original")
            fits_boot, df_boot = _one_round(
                "pair_weight_geommean_1cond_mean", "pairStdSuppGeom_1cond_mean",
                "singleWeightGeom_periph_1cond_mean", "sigCorr base — Bootstrapped")

        finally:
            # Always restore original specs
            _MODEL_SPECS = _orig_specs
            _LRT_PAIRS = _orig_lrts

        return {"orig": fits_orig, "boot": fits_boot}, df_orig, df_boot


    def _plot_extended_comparison(fits_orig, fits_boot):
        keys   = [k for k, _, _ in _MODEL_SPECS if fits_orig.get(k)]
        labels = [fits_orig[k]["label"].split("] ")[-1] for k in keys]
        colors = [_MODEL_COLORS[k] for k in keys]
        x, w   = np.arange(len(keys)), 0.32

        fig, axes = plt.subplots(1, 3, figsize=(20, 5))

        # Panel A: Marginal R²
        ax = axes[0]
        vals_o = [fits_orig[k]["r2m"] if fits_orig.get(k) else np.nan for k in keys]
        vals_b = [fits_boot[k]["r2m"] if fits_boot.get(k) else np.nan for k in keys]
        ax.bar(x - w/2, vals_o, w, color=colors, alpha=0.9,
               label="Original",     edgecolor="white")
        ax.bar(x + w/2, vals_b, w, color=colors, alpha=0.45,
               label="Bootstrapped", edgecolor="white", hatch="//")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=7)
        ax.set_ylabel("Marginal R²")
        ax.set_ylim(0, max(vals_o + vals_b) * 1.35)
        ax.set_title("Marginal R²\nacross all models")
        ax.legend(fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)

        # Panel B: Fixed-effects R²
        ax = axes[1]
        vals_o = [fits_orig[k]["r_fe"]**2 if fits_orig.get(k) else np.nan for k in keys]
        vals_b = [fits_boot[k]["r_fe"]**2 if fits_boot.get(k) else np.nan for k in keys]
        ax.bar(x - w/2, vals_o, w, color=colors, alpha=0.9,
               label="Original",     edgecolor="white")
        ax.bar(x + w/2, vals_b, w, color=colors, alpha=0.45,
               label="Bootstrapped", edgecolor="white", hatch="//")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=7)
        ax.set_ylabel("Fixed-effects R²")
        ax.set_ylim(0, max(vals_o + vals_b) * 1.35)
        ax.set_title("Fixed-effects-only R²\n(honest population fit)")
        ax.legend(fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)

        # Panel C: w_c, s_c, fr_c coefficients in models that contain all three
        ax = axes[2]
        cmp_keys = [k for k in ["M1_weight", "M2_fr", "M3_supp",
                                  "M4_wt_fr", "M5_wt_supp", "M6_fr_supp", "M7_all"]
                    if fits_orig.get(k)]
        xi         = np.arange(len(cmp_keys))
        cmp_labels = [fits_orig[k]["label"].split("] ")[-1].split(": ")[1]
                      for k in cmp_keys]

        offsets = [(-w, "w_c", "#2ca02c"), (0,  "s_c", "#d62728"), (+w, "fr_c", "#ff7f0e")]
        bar_w   = w * 0.9
        for offset, cname, color in offsets:
            for suf, fits_d, alpha in [("_orig", fits_orig, 0.9), ("_boot", fits_boot, 0.45)]:
                coefs = [fits_d[k]["res"].fe_params.get(cname, np.nan)
                         if fits_d.get(k) else np.nan for k in cmp_keys]
                ses   = [fits_d[k]["res"].bse_fe.get(cname, np.nan)
                         if fits_d.get(k) else np.nan for k in cmp_keys]
                lbl   = f"{'weight' if cname=='w_c' else 'supp' if cname=='s_c' else 'FR'}" \
                        f" {'orig' if suf=='_orig' else 'boot'}"
                ax.bar(xi + offset, coefs, bar_w, color=color, alpha=alpha,
                       label=lbl, edgecolor="white",
                       hatch="" if suf == "_orig" else "//")
                ax.errorbar(xi + offset, coefs, yerr=ses, fmt="none",
                            color="black" if suf == "_orig" else "#555",
                            capsize=2, lw=1.0)

        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_xticks(xi)
        ax.set_xticklabels(cmp_labels, rotation=25, ha="right", fontsize=7)
        ax.set_ylabel("Coefficient  (± SE)")
        ax.set_title("weight / stdSupp / FR coefficients\n(orig vs boot)")
        ax.legend(fontsize=6, frameon=False, ncol=3)
        ax.spines[["top", "right"]].set_visible(False)

        fig.suptitle(
            "LME: weight vs FR vs stdSupp → rsc  "
            "(8 models × original / bootstrapped predictors)",
            fontsize=12, fontweight="bold", y=1.02)
        plt.tight_layout()
        plt.savefig("lme_extended_comparison.png", dpi=150, bbox_inches="tight")
        plt.show()

# ══════════════════════════════════════════════════════════════════════════════
#  ENTRY POINT — run this block after corr_df is built
# ══════════════════════════════════════════════════════════════════════════════

# Build stim subsets (used by both bootstraps)
all_resp_arrays = [
    prefNormalizedWithBase,     # (n_neurons, 5)
    nonprefNormalizedWithBase,  # (n_neurons, 5)
    pnNormalizedWithBase,       # (n_neurons, 4)
    npNormalizedWithBase,       # (n_neurons, 4)
    ppNormalizedWithBase,       # (n_neurons, 4)
    nnNormalizedWithBase,       # (n_neurons, 4)
]
stim_subsets = build_stim_subsets(
    contrast_centerSpon, contrast_peripherySpon,
    locationsSpon, stim_type_centerSpon, stim_type_peripherySpon
)

# Step 1: heuristic weight bootstrap
boot_df, corr_df_subsampled = run_precision_bootstrap_v2(
    corr_df=corr_df, masterGoodUnits=masterGoodUnits,
    rfWeightNeuronFits=rfWeightNeuronFits,
    all_resp_arrays=all_resp_arrays, stim_subsets=stim_subsets, n_boot=100,
)

# Step 2: standard suppression bootstrap
corr_df_full = run_stdsupp_bootstrap(
    corr_df=corr_df_subsampled, masterGoodUnits=masterGoodUnits,
    stdNormFits=stdNormFits,
    all_resp_arrays=all_resp_arrays, stim_subsets=stim_subsets, n_boot=100,
)

# Verify all required columns present
expected_cols = ["pair_weight_geommean", "pairStdSuppGeom",
                 "pair_weight_geommean_1cond_mean",
                 "singleWeightGeom_periph_1cond_mean",
                 "pairStdSuppGeom_1cond_mean"]
for col in expected_cols:
    status = "✓" if col in corr_df_full.columns else "✗ MISSING"
    print(f"  {status}  {col}")

# Step 3: extended LME (8 models, weight + FR + stdSupp and all combinations)
all_fits, comp, df_orig, df_boot = run_lme_extended(corr_df_full)

fits_sc, df_sc_orig, df_sc_boot = run_lme_sigcorr(corr_df_full)

########################################################################################################################
########################################################## Steiger Analysis ############################################
########################################################################################################################

from scipy.stats import pearsonr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def meng_steiger_z(r_weight, r_rate, r_wr, n):
    """
    Meng, Rosenthal & Rubin (1992) test for difference between two
    dependent correlations sharing one variable (rsc).

    r_weight : partial r(rsc, weight | BC)
    r_rate   : partial r(rsc, rate   | BC)
    r_wr     : r(weight, rate) — correlation between the two predictors
    n        : number of observations
    """
    # Average r for determining R²
    r_bar = (r_weight ** 2 + r_rate ** 2) / 2
    # Covariance correction factor
    f = (1 - r_wr) / (2 * (1 - r_bar))
    h = (1 - f * r_bar) / (1 - r_bar)
    # Z-transformed correlations
    z_weight = np.arctanh(r_weight)
    z_rate = np.arctanh(r_rate)
    # Test statistic
    z_diff = (z_weight - z_rate) * np.sqrt((n - 3) / (2 * (1 - r_wr) * h))
    p = 2 * (1 - stats.norm.cdf(abs(z_diff)))
    return z_diff, p


def plot_steiger_scatter(pair_df, alpha=0.05):
    """
    Steiger scatter with two significance boundary lines.
    Points above upper line: weight significantly > rate
    Points between lines:    no significant difference
    Points below lower line: rate significantly > weight
    """
    from scipy import stats

    fig, ax = plt.subplots(figsize=(6, 6))

    z_w = pair_df["z_weight"].values
    z_r = pair_df["z_rate"].values
    n_vals = pair_df["n"].values

    # ── Classify each pair by Steiger significance ────────────────────────────
    colors = []
    for _, row in pair_df.iterrows():
        if row["p_steiger"] < alpha and row["z_diff"] > 0:
            colors.append("#d62728")   # red: weight wins
        elif row["p_steiger"] < alpha and row["z_diff"] < 0:
            colors.append("#1f77b4")   # blue: rate wins
        else:
            colors.append("#7f7f7f")   # gray: no significant difference

    colors = np.array(colors)
    n_weight_wins = (colors == "#d62728").sum()
    n_rate_wins   = (colors == "#1f77b4").sum()
    n_neither     = (colors == "#7f7f7f").sum()
    n_total       = len(pair_df)

    # Plot points
    for color, label in [("#d62728", f"weight > rate  (n={n_weight_wins}/{n_total})"),
                          ("#7f7f7f", f"no difference  (n={n_neither}/{n_total})"),
                          ("#1f77b4", f"rate > weight  (n={n_rate_wins}/{n_total})")]:
        mask = colors == color
        ax.scatter(z_r[mask], z_w[mask], s=30, color=color,
                   alpha=0.7, edgecolors="white", lw=0.4, label=label, zorder=3)

    # ── Significance boundary lines ───────────────────────────────────────────
    # The boundary is where Steiger Z = ±z_crit.
    # At the boundary: z_weight - z_rate = z_crit * SE
    # SE depends on n and r_wr (correlation between predictors).
    # Since r_wr varies per pair, we use the median r_wr for the boundary lines.
    z_crit   = stats.norm.ppf(1 - alpha / 2)
    r_wr_med = pair_df["r_wr"].median()
    n_med    = int(pair_df["n"].median())

    # From Meng et al: z_diff = (z_w - z_r) * sqrt((n-3) / (2*(1-r_wr)*h))
    # where h = (1 - f*r_bar) / (1 - r_bar), f = (1-r_wr)/(2*(1-r_bar))
    # At boundary z_diff = ±z_crit, solve for (z_w - z_r):
    # This requires r_bar which depends on both z values, so we approximate
    # using a representative r_bar from the data median
    r_bar_med = (pair_df["r_weight"].abs().median()**2 +
                 pair_df["r_rate"].abs().median()**2) / 2
    f_med  = (1 - r_wr_med) / (2 * (1 - r_bar_med))
    h_med  = (1 - f_med * r_bar_med) / (1 - r_bar_med)
    se_med = np.sqrt((2 * (1 - r_wr_med) * h_med) / (n_med - 3))
    delta  = z_crit * se_med   # critical difference in Z units

    # Boundary lines: z_weight = z_rate ± delta
    lims = [min(z_r.min(), z_w.min()) - 0.3,
            max(z_r.max(), z_w.max()) + 0.3]
    x_line = np.array(lims)

    ax.plot(x_line, x_line + delta, color="black", lw=1.2, ls="-",
            zorder=2)   # upper line: weight wins above this
    ax.plot(x_line, x_line - delta, color="black", lw=1.2, ls="-",
            zorder=2)   # lower line: rate wins below this

    # ── Annotations ───────────────────────────────────────────────────────────
    ax.set_xlabel("Z-transformed rate correlation (Z_rate | sigCorr)", fontsize=11)
    ax.set_ylabel("Z-transformed weight correlation (Z_weight | sigCorr)", fontsize=11)
    ax.set_xlim(lims);
    ax.set_ylim(lims)
    ax.set_aspect("equal")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)

    # Compute wilcoxon first, then set title once
    stat, p_wilcox = stats.wilcoxon(pair_df["z_diff"])
    ax.set_title(f"Steiger test: weight vs rate per pair (controlling for signal correlation)\n"
                 f"Wilcoxon on Z_diff: p={p_wilcox:.3e}", fontsize=10)

    plt.tight_layout()
    # plt.savefig("steiger_scatter_boundary.png", dpi=150, bbox_inches="tight")
    plt.show()
    return fig


def run_steiger_deduplicated(df_orig):
    from scipy import stats

    results = []
    pair_ids = df_orig["pair_id"].unique()

    for pid in pair_ids:
        sub = df_orig[df_orig["pair_id"] == pid].dropna(
            subset=["rsc", "w_c", "fr_c", "sc_c"])  # <-- sc_c

        if len(sub) < 5:
            continue

        rsc  = sub["rsc"].values
        w_c  = sub["w_c"].values
        fr_c = sub["fr_c"].values
        sc_c = sub["sc_c"].values                    # <-- sc_c
        n    = len(rsc)

        def resid_on_sc(y, sc):                      # <-- renamed
            X = np.column_stack([np.ones(len(sc)), sc, sc ** 2])
            coef, *_ = np.linalg.lstsq(X, y, rcond=None)
            return y - X @ coef

        rsc_r = resid_on_sc(rsc, sc_c)              # <-- sc_c
        w_r   = resid_on_sc(w_c,  sc_c)             # <-- sc_c
        fr_r  = resid_on_sc(fr_c, sc_c)             # <-- sc_c

        r_weight, _ = pearsonr(rsc_r, w_r)
        r_rate,   _ = pearsonr(rsc_r, fr_r)
        r_wr,     _ = pearsonr(w_r,   fr_r)

        z_w = np.arctanh(np.clip(r_weight, -0.9999, 0.9999))
        z_r = np.arctanh(np.clip(r_rate,   -0.9999, 0.9999))

        z_diff, p_steiger = meng_steiger_z(r_weight, r_rate, r_wr, n)

        results.append(dict(
            pair_id=pid,
            n=n,
            r_weight=r_weight,
            r_rate=r_rate,
            r_wr=r_wr,
            z_weight=z_w,
            z_rate=z_r,
            z_diff=z_diff,
            p_steiger=p_steiger,
        ))

    pair_df = pd.DataFrame(results)

    # ── Summary counts ────────────────────────────────────────────────────────
    n_total = len(pair_df)
    n_rate_wins = ((pair_df["p_steiger"] < 0.05) & (pair_df["z_diff"] < 0)).sum()
    n_weight_wins = ((pair_df["p_steiger"] < 0.05) & (pair_df["z_diff"] > 0)).sum()
    n_neither = ((pair_df["p_steiger"] >= 0.05)).sum()

    # FR exceeds weight (regardless of significance) — matches the 53.1% framing
    n_fr_exceeds = (pair_df["z_diff"] < 0).sum()  # z_diff = z_weight - z_rate, so negative = rate wins

    print(f"\n{'=' * 55}")
    print(f"  Steiger test summary (n={n_total} pairs)")
    print(f"{'=' * 55}")
    print(f"  FR > weight (z_diff < 0):        {n_fr_exceeds:>4}  ({100 * n_fr_exceeds / n_total:.1f}%)")
    print(
        f"  Weight > FR (z_diff > 0):        {n_total - n_fr_exceeds:>4}  ({100 * (n_total - n_fr_exceeds) / n_total:.1f}%)")
    print(f"  --- Significance breakdown ---")
    print(f"  FR significantly > weight:       {n_rate_wins:>4}  ({100 * n_rate_wins / n_total:.1f}%)")
    print(f"  Weight significantly > FR:       {n_weight_wins:>4}  ({100 * n_weight_wins / n_total:.1f}%)")
    print(f"  No significant difference:       {n_neither:>4}  ({100 * n_neither / n_total:.1f}%)")

    from scipy import stats
    stat, p_wilcox = stats.wilcoxon(pair_df["z_diff"])
    print(f"\n  Wilcoxon on z_diff: stat={stat:.3f}, p={p_wilcox:.4e}")
    print(f"{'=' * 55}")

    plot_steiger_scatter(pair_df, alpha=0.05)
    return pair_df


# ── Run it ────────────────────────────────────────────────────────────────────
pair_df = run_steiger_deduplicated(df_sc_boot)

########################################################################################################################
########################################################################################################################

