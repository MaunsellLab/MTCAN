# Imports
from usefulFns import *
from normalizationFunctions import *

p0 = time.time()

# Master list (both monkeys)
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

unitListCopy = np.copy(unitList)
unitListCopy = np.array(unitListCopy)

# population arrays
individualSpikeCountMat = []
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

    # add good units from each session to population arrays
    # for unit in sessionGoodUnits:
    corrUnits = []
    for unit in units:
        count = np.where(units == unit)[0][0]
        masterAllUnits.append(f'{seshDate}_{unit}')

        # add the spikeCountMat of good units to a master list so that individual trial spike counts
        # are accessible outside the individual file for loop (at the population level)
        individualSpikeCountMat.append(spikeCountMat[count, :-1, :])

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
    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
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
###################################### Direction Tuning and Correlation Analysis #######################################
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
        [ np.inf, np.inf, 2*np.pi, 1000.0]
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


# get the sessions with more than one unit and create a new fileList
# --- Extract dates from unitList ---
dates = [u.split('_')[0] for u in unitList]
# --- Count occurrences of each date ---
date_counts = Counter(dates)
# --- Filter dates with more than one unit ---
multi_unit_dates = [d for d, count in date_counts.items() if count > 1]
# --- Add prefixes depending on year ---
date_labels = []
for d in multi_unit_dates:
    if d.startswith("23"):
        date_labels.append(f"Meetz_{d}")
    elif d.startswith("24"):
        date_labels.append(f"Akshan_{d}")
    else:
        date_labels.append(d)  # fallback if another year shows up

# go through each session with multiple units and extract the direction tuning from the cells of interest
os.chdir('../')
similarity_records = []  # will accumulate dicts across all sessions
for date in date_labels:
    monkeyName, sessionDate = date.split('_')
    os.chdir(f'{monkeyName}/{sessionDate}/Direction Tuning')
    dirData = np.load('unitsDirTuningMat.npy')  # shape: (n_units_today, n_dirs)

    # Units present today (strings like '230607_77'), and the good ones
    units_today = [u for u in masterAllUnits if u.startswith(sessionDate + "_")]
    goodUnits_today = [u for u in masterGoodUnits if u.startswith(sessionDate + "_")]

    # Build a mapping from unit string -> row index in dirData
    # ASSUMPTION: rows of dirData correspond to units_today in order.
    unit_to_row = {u: i for i, u in enumerate(units_today)}

    # Keep only good units that actually exist in today's data
    goodUnits_today = [u for u in goodUnits_today if u in unit_to_row]

    # If fewer than 2 good units, nothing to compare
    if len(goodUnits_today) < 2:
        continue

    # Angle grid
    n_neurons_today, n_dirs = dirData.shape
    angles_deg = np.linspace(0, 360, n_dirs, endpoint=False)
    angles_rad = np.deg2rad(angles_deg)

    # Fit von Mises for good units only, store PDs aligned to goodUnits_today order
    pref_dirs_rad = np.zeros(len(goodUnits_today))
    for gi, u in enumerate(goodUnits_today):
        idx = unit_to_row[u]
        yi = dirData[idx, :]
        _, _, prefD, _ = fit_von_mises_single(yi, angles_rad)
        pref_dirs_rad[gi] = prefD

    pref_dirs_deg = np.rad2deg(pref_dirs_rad)

    # Similarity among good units only
    # cos(Δθ) in [-1,1]; 1 = identical, -1 = 180° apart. SI01 in [0,1].
    c = np.cos(pref_dirs_rad).reshape(-1, 1)
    s = np.sin(pref_dirs_rad).reshape(-1, 1)
    SI = (c @ c.T) + (s @ s.T)            # (ngood x ngood), includes self-diagonal=1
    SI01 = 0.5 * (SI + 1.0)

    # Record only unique pairs (i < j), skip self-pairs
    for i, j in combinations(range(len(goodUnits_today)), 2):
        similarity_records.append({
            "date": sessionDate,
            "monkey": monkeyName,
            "unit_i": goodUnits_today[i],
            "unit_j": goodUnits_today[j],
            "PD_i_deg": float(pref_dirs_deg[i]),
            "PD_j_deg": float(pref_dirs_deg[j]),
            "SI": float(SI[i, j]),
            "SI01": float(SI01[i, j]),
        })

    # change directory back
    os.chdir('../../../')


########################################################################################################################
#################################################### Correlation Analysis ##############################################
########################################################################################################################
def shared_geom_share_rectified(w1, w2, sigma1, sigma2):
    eps = 1e-9
    w1p = max(w1, 0.0)
    w2p = max(w2, 0.0)
    s1 = w1p / (1 + w1p + sigma1 + eps)   # in [0,1)
    s2 = w2p / (1 + w2p + sigma2 + eps)
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


# === NEW: map paired index i to single indices (center/periphery) and site ===
def single_indices_for_pair(i):
    """
    Return (center_single_idx, periph_single_idx, site)
    - For i in 10..17 (center=N paired): center_single_idx=9
    - For i in 19..26 (center=P paired): center_single_idx=18
    - periph_single_idx alternates N/P by within-site parity
      and yields {1,2,3,4,5,6,7,8} depending on site and parity.
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


# Helper to add a record (extended to include singles)
def add_record(block, i, comb, pairStimCorr, pairWeightGeomMean, n1Sel, n2Sel,
               sim_SI=None, sim_SI01=None, delta_PD_abs_deg=None, avg_pair_weight_2w=None,
               monkey=None, singleCorr_center=None, singleCorr_periph=None,
               singleWeightGeom_center=None, singleWeightGeom_periph=None):
    u1, u2 = comb
    date = u1.split("_")[0]  # both units should share the same date
    corr_records.append({
        "date": date,
        "unit_i": u1,
        "unit_j": u2,
        "block": block,              # e.g., "centerN_offsetsPN" or "centerP_offsetsPN"
        "cond_index": i,             # 10..17 or 19..26
        "corr": float(pairStimCorr),
        "pair_weight_geommean": float(pairWeightGeomMean),
        "n1Sel": float(n1Sel),
        "n2Sel": float(n2Sel),
        "pairSel_geom": float(np.sign(n1Sel)*np.sign(n2Sel)*np.sqrt(abs(n1Sel)*abs(n2Sel))),
        "similarity_SI": (None if sim_SI is None else float(sim_SI)),
        "similarity_SI01": (None if sim_SI01 is None else float(sim_SI01)),
        "delta_PD_abs_deg": (None if delta_PD_abs_deg is None else float(delta_PD_abs_deg)),
        "avg_pair_weight_2w": (None if avg_pair_weight_2w is None else float(avg_pair_weight_2w)),
        "monkey": (None if monkey is None else str(monkey)),
        # NEW singles:
        "singleCorr_center": (None if singleCorr_center is None else float(singleCorr_center)),
        "singleCorr_periph": (None if singleCorr_periph is None else float(singleCorr_periph)),
        "singleWeightGeom_center": (None if singleWeightGeom_center is None else float(singleWeightGeom_center)),
        "singleWeightGeom_periph": (None if singleWeightGeom_periph is None else float(singleWeightGeom_periph)),
    })


# stim index mat (not strictly needed now, but kept for reference)
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

# Faster index maps for master arrays to avoid repeated np.where calls
good_index = {u: i for i, u in enumerate(masterGoodUnits)}
corr_index = {u: i for i, u in enumerate(masterCorrUnits)}

# Accumulators
corr_records = []   # list of dict rows
masterCorr = []
masterPairWeight = []
masterPairSel = []
sameSelCorr = []
sameSelWeight = []
oppositeSelCorr = []
oppositeSelWeight = []

for comb in all_combs:
    if comb[0] not in good_index or comb[1] not in good_index:
        continue

    n1 = good_index[comb[0]]
    n2 = good_index[comb[1]]

    # similarity/PD lookup (order-agnostic)
    key = tuple(sorted((comb[0], comb[1])))
    monkey_name = monkey_lookup.get(key)
    sim_val = sim_SI_lookup.get(key)
    sim01_val = sim_SI01_lookup.get(key)
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
        pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

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

        # Paired weight geom-mean (existing)
        pairWeightGeomMean = shared_geom_share_rectified(w1=n1w, w2=n2w, sigma1=sig1, sigma2=sig2)
        masterPairWeight.append(pairWeightGeomMean)
        masterCorr.append(pairStimCorr)

        # --- NEW: single-stim correlations ---
        n1_single_center = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, center_single_idx]
        n2_single_center = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, center_single_idx]
        singleCorr_center, _, _ = pairCorrExclude3SD(n1_single_center, n2_single_center)

        n1_single_periph = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, periph_single_idx]
        n2_single_periph = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, periph_single_idx]
        singleCorr_periph, _, _ = pairCorrExclude3SD(n1_single_periph, n2_single_periph)

        # --- NEW: single-stim weight geom-means ---
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

        # Add record with NEW singles fields
        add_record(block="centerN_offsetsPN", i=i, comb=comb,
                   pairStimCorr=pairStimCorr,
                   pairWeightGeomMean=pairWeightGeomMean,
                   n1Sel=n1Sel, n2Sel=n2Sel,
                   sim_SI=sim_val, sim_SI01=sim01_val,
                   delta_PD_abs_deg=deltaPD, avg_pair_weight_2w=avg_pair_weight_2w,
                   monkey=monkey_name,
                   singleCorr_center=singleCorr_center,
                   singleCorr_periph=singleCorr_periph,
                   singleWeightGeom_center=singleWeightGeom_center,
                   singleWeightGeom_periph=singleWeightGeom_periph)

    # =========================
    # Block 2: i = 19..26 (center=P offsets P or N)
    # =========================
    for i in range(19, 27):
        n1SpkMatIdx = corr_index[comb[0]]
        n2SpkMatIdx = corr_index[comb[1]]
        n1SpikeMat = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, i]
        n2SpikeMat = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, i]
        pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

        # Correct single indices and site
        center_single_idx, periph_single_idx, site = single_indices_for_pair(i)

        w0 = 1.0
        n1w = rfWeightNeuronFits[n1][2 + site]
        sig1 = rfWeightNeuronFits[n1][6]
        n2w = rfWeightNeuronFits[n2][2 + site]
        sig2 = rfWeightNeuronFits[n2][6]
        if n1w <= 0 or n2w <= 0:
            continue

        pairWeightGeomMean = shared_geom_share_rectified(w1=n1w, w2=n2w, sigma1=sig1, sigma2=sig2)
        masterPairWeight.append(pairWeightGeomMean)
        masterCorr.append(pairStimCorr)

        # --- NEW: single-stim correlations ---
        n1_single_center = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, center_single_idx]
        n2_single_center = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, center_single_idx]
        singleCorr_center, _, _ = pairCorrExclude3SD(n1_single_center, n2_single_center)

        n1_single_periph = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, periph_single_idx]
        n2_single_periph = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, periph_single_idx]
        singleCorr_periph, _, _ = pairCorrExclude3SD(n1_single_periph, n2_single_periph)

        # --- NEW: single-stim weight geom-means ---
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

        add_record(block="centerP_offsetsPN", i=i, comb=comb,
                   pairStimCorr=pairStimCorr,
                   pairWeightGeomMean=pairWeightGeomMean,
                   n1Sel=n1Sel, n2Sel=n2Sel,
                   sim_SI=sim_val, sim_SI01=sim01_val,
                   delta_PD_abs_deg=deltaPD, avg_pair_weight_2w=avg_pair_weight_2w,
                   monkey=monkey_name,
                   singleCorr_center=singleCorr_center,
                   singleCorr_periph=singleCorr_periph,
                   singleWeightGeom_center=singleWeightGeom_center,
                   singleWeightGeom_periph=singleWeightGeom_periph)

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
    "nonPrefSupp": corr_df['pair_weight_geommean'].values,
    "selectivity": corr_df['delta_PD_abs_deg'].values,
    "pairedCorr": corr_df['corr'].values,
    "singleCorrLoc0": corr_df['singleCorr_center'].values,
    "singleCorrLoc1": corr_df['singleCorr_periph'].values,
    "singleCenterNonPrefSupp": corr_df['singleWeightGeom_center'].values,
    "singlePeriphNonPrefSupp": corr_df['singleWeightGeom_periph'].values,
    "monkey": corr_df['monkey'].values.astype('U')  # save as strings
})

# paired
x = corr_df['pair_weight_geommean'].values
y = corr_df['delta_PD_abs_deg'].values
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


# binned plot
for filler in range(1):
    # 1) filter to 0–2 on x
    mask = (masterPairWeight >= 0) & (masterPairWeight <= 2)
    x = masterPairWeight[mask]
    y = masterCorr[mask]

    # x = sameSelWeight
    # y = sameSelCorr

    # x = oppositeSelWeight
    # y = oppositeSelCorr

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
    # x = masterPairWeight[mask]
    # y = masterCorr[mask]

    # x = sameSelWeight
    # y = sameSelCorr

    x = oppositeSelWeight
    y = oppositeSelCorr

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

########################################################################################################################
################################################# extra code ###########################################################
########################################################################################################################

    # ------- Block 3: center single stim N or P (i=9 or i=18) -------
    for i in [9, 18]:
        n1SpkMatIdx = corr_index[comb[0]]
        n2SpkMatIdx = corr_index[comb[1]]
        n1SpikeMat = masterCorrUnitsSpikeCountMat[n1SpkMatIdx][:, i]
        n2SpikeMat = masterCorrUnitsSpikeCountMat[n2SpkMatIdx][:, i]
        pairSingleCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)

        if i == 9:
            indx = 0
        if i == 18:
            indx = 1
        w0 = 1.0
        n1w = 0
        n2w = 0
        sig1 = rfWeightNeuronFits[n1][6]
        sig2 = rfWeightNeuronFits[n2][6]

        pairWeightGeomMean = shared_geom_share_rectified(w1=n1w, w2=n2w, sigma1=sig1, sigma2=sig2)
        masterPairWeight.append(pairWeightGeomMean)
        masterCorr.append(pairSingleCorr)

        # --- NEW: avg of each neuron's two weights (nW, w0=1), then avg across neurons ---
        avg_n1_2w = 0.5 * (float(n1w) + w0)
        avg_n2_2w = 0.5 * (float(n2w) + w0)
        avg_pair_weight_2w = 0.5 * (avg_n1_2w + avg_n2_2w)

        # Add fully-labeled record (now with delta_PD_abs_deg and avg_pair_weight_2w)
        add_record(block="centerP_or_N", i=i, comb=comb,
                   pairSingleCorr=pairSingleCorr,
                   pairWeightGeomMean=pairWeightGeomMean,
                   delta_PD_abs_deg=deltaPD, avg_pair_weight_2w=avg_pair_weight_2w,
                   monkey=monkey_name)

