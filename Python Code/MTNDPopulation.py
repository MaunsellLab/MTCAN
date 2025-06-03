"""
MTND population analysis script. This script will take a list of 'good units'
across days/sessions and add their normalized responses to a population array.
From there, this script plots how normalization changes with distance at a population
level.

Chery - June 2024
"""

# Imports
from usefulFns import *
from normalizationFunctions import *

p0 = time.time()


## null resp > baseline and non-pref for Meetz + non-pref for Akshan

# # # file list (MEETZ) (null and non-pref with pref)
# fileList = ['Meetz_230607', 'Meetz_230608', 'Meetz_230612', 'Meetz_230613', 'Meetz_230615',
#             'Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628', 'Meetz_230630',
#             'Meetz_230707', 'Meetz_230710', 'Meetz_230711', 'Meetz_230718', 'Meetz_230719',
#             'Meetz_230720']
#
# # # unit list with null resp > baseline (MEETZ)
# unitList = ['230607_147', '230607_146', '230607_144', '230607_126','230607_125',
#             '230608_143', '230608_145', '230608_149', '230612_213', '230613_187',
#             '230615_83', '230615_166', '230620_58', '230626_131', '230627_147',
#             '230627_141', '230627_106', '230627_71', '230627_70', '230627_14',
#             '230628_143', '230630_50', '230707_140', '230710_137', '230711_45',
#             '230711_50', '230718_178', '230719_147', '230719_156', '230720_149']

# # file list, only non-pref + pref (Akshan)
# fileList = ['Akshan_240530', 'Akshan_240603', 'Akshan_240606', 'Akshan_240607',
#             'Akshan_240611', 'Akshan_240613', 'Akshan_240628', 'Akshan_240703',
#             'Akshan_240704', 'Akshan_240705', 'Akshan_240709', 'Akshan_241014',
#             'Akshan_241024']
#
# unitList = ['240530_55', '240603_167', '240606_176', '240607_137', '240607_170',
#             '240607_179', '240611_136', '240613_138', '240628_19', '240703_19',
#             '240703_91', '240704_86', '240705_74', '240705_162', '240709_53',
#             '240709_108', '241014_3','241014_61', '241014_71', '241024_2',
#             '241024_21', '241024_25']

# # Master list (both monkeys)
fileList = ['Meetz_230607', 'Meetz_230608', 'Meetz_230612', 'Meetz_230613', 'Meetz_230615',
            'Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628', 'Meetz_230630',
            'Meetz_230707', 'Meetz_230710', 'Meetz_230711', 'Meetz_230718', 'Meetz_230719',
            'Meetz_230720', 'Akshan_240530', 'Akshan_240603', 'Akshan_240606', 'Akshan_240607',
            'Akshan_240610', 'Akshan_240611', 'Akshan_240613', 'Akshan_240628', 'Akshan_240703',
            'Akshan_240704', 'Akshan_240705', 'Akshan_240709', 'Akshan_241014', 'Akshan_241024']

# # original
# unitList = ['230607_147', '230607_146', '230607_144', '230607_126', '230607_125',
#             '230608_143', '230608_145', '230608_149', '230612_213', '230613_187',
#             '230615_83', '230615_166', '230620_58', '230626_131', '230627_147',
#             '230627_141', '230627_106', '230627_71', '230627_70', '230627_14',
#             '230628_143', '230630_50', '230707_140', '230710_137', '230711_45',
#             '230711_50', '230718_178', '230719_147', '230719_156', '230720_149',
#             '240530_55', '240603_167', '240606_176', '240607_137', '240607_170',
#             '240607_179', '240611_136', '240613_138', '240628_19', '240703_19',
#             '240703_91', '240704_86', '240705_74', '240705_162', '240709_53',
#             '240709_108', '241014_3', '241014_61', '241014_71', '241024_2',
#             '241024_21', '241024_25']

# some new units added
# unitList = ['230607_147', '230607_146', '230607_144', '230607_126', '230607_125',
#             '230607_80', '230607_77', '230608_112', '230608_117', '230608_140',
#             '230608_143', '230608_145', '230608_149', '230612_213', '230613_187',
#             '230615_83', '230615_166', '230615_164', '230615_74', '230620_58',
#             '230626_131', '230626_58', '230626_65', '230626_64', '230627_147',
#             '230627_141', '230627_106', '230627_71', '230627_70', '230627_14',
#             '230627_79', '230627_148', '230627_112', '230627_93', '230628_143',
#             '230630_50', '230630_138', '230707_140', '230707_23', '230707_144',
#             '230707_143', '230710_137', '230711_45', '230711_50', '230711_55',
#             '230718_178', '230719_147', '230719_156', '230719_155', '230720_149',
#             '240530_55', '240603_167', '240606_176', '240607_137', '240607_170',
#             '240607_179', '240611_136', '240613_138', '240628_19', '240703_19',
#             '240703_91', '240704_86', '240705_74', '240705_162', '240709_53',
#             '240709_108', '241014_3', '241014_61', '241014_71', '241024_2',
#             '241024_21', '241024_25']


# # unitList with only one optimized cell/day ignoring the ancillary cells
# unitList = ['230607_144', '230612_213', '230613_187', '230615_83', '230620_58',
#             '230626_131', '230627_147', '230627_141', '230628_143', '230630_50',
#             '230707_140', '230710_137', '230711_50', '230718_178', '230719_147',
#             '230720_149', '240530_55', '240603_167', '240606_176', '240607_179',
#             '240611_136', '240613_138', '240628_19', '240703_91', '240704_86',
#             '240705_74', '240709_108', '241014_61',  '241024_2']

# some new units added + most of Akshan's sessions run through KS4
unitList = ['230607_147', '230607_146', '230607_144', '230607_126', '230607_125',
            '230607_80', '230607_77', '230608_112', '230608_117', '230608_140',
            '230608_143', '230608_145', '230608_149', '230612_213', '230613_187',
            '230615_83', '230615_166', '230615_164', '230615_74', '230620_58',
            '230626_131', '230626_58', '230626_65', '230626_64', '230627_147',
            '230627_141', '230627_106', '230627_71', '230627_70', '230627_14',
            '230627_79', '230627_148', '230627_112', '230627_93', '230628_143',
            '230630_50', '230630_138', '230707_140', '230707_23', '230707_144',
            '230707_143', '230710_137', '230711_45', '230711_50', '230711_55',
            '230718_178', '230719_147', '230719_156', '230719_155', '230720_149',
            '240530_54', '240603_118', '240606_46', '240606_52', '240606_78',
            '240606_79', '240606_80', '240606_99', '240607_137', '240607_170',
            '240607_179', '240611_21', '240611_28', '240611_29', '240611_36',
            '240611_37', '240611_46', '240611_81', '240611_86', '240611_97',
            '240611_101', '240611_103', '240611_110', '240613_12', '240628_49',
            '240703_7', '240703_57', '240703_66', '240703_72', '240704_24',
            '240705_56', '240705_60', '240709_53', '240709_108', '241014_3',
            '241014_61', '241014_71', '241024_2', '241024_21', '241024_25']


########################################################################################################################

################################################ PREF + NON PREF ONLY ##################################################

########################################################################################################################

# # # pref + nonpref (Meetz) (n=18)
# fileList = ['Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628',
#             'Meetz_230630', 'Meetz_230707', 'Meetz_230710', 'Meetz_230711',
#             'Meetz_230718', 'Meetz_230719', 'Meetz_230720']
#
# unitList = ['230620_58', '230626_131', '230627_147', '230627_141', '230627_106',
#             '230627_71', '230627_70', '230627_14', '230628_143', '230630_50',
#             '230707_140', '230710_137', '230711_45', '230711_50', '230718_178',
#             '230719_147', '230719_156', '230720_149']

# pref + nonpref (Akshan)
# fileList = ['Akshan_240530', 'Akshan_240603', 'Akshan_240606', 'Akshan_240607',
#             'Akshan_240611', 'Akshan_240613', 'Akshan_240628', 'Akshan_240703',
#             'Akshan_240704', 'Akshan_240705', 'Akshan_240709', 'Akshan_241014',
#             'Akshan_241024']
#
# unitList = ['240530_55', '240603_167', '240606_176', '240607_137', '240607_170',
#             '240607_179', '240611_136', '240613_138', '240628_19', '240703_19',
#             '240703_91', '240704_86', '240705_74', '240705_162', '240709_53',
#             '240709_108', '241014_3','241014_61', '241014_71', '241024_2',
#             '241024_21', '241024_25']

## Akshan File list with most sessions run through KS4
# unitList = ['240530_54', '240603_118', '240606_46', '240606_52', '240606_78',
#             '240606_79', '240606_80', '240606_99', '240607_137', '240607_170',
#             '240607_179', '240611_21', '240611_28', '240611_29', '240611_36',
#             '240611_37', '240611_46',  '240611_81', '240611_86', '240611_97',
#             '240611_101', '240611_103', '240611_110', '240613_12', '240628_49',
#             '240703_7', '240703_57', '240703_66','240703_72', '240704_24',
#             '240705_56', '240705_60', '240709_53', '240709_108', '241014_3',
#             '241014_61', '241014_71', '241024_2', '241024_21', '241024_25']


# # both monkeys (Meetz + Akshan)
# fileList = ['Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628',
#             'Meetz_230630', 'Meetz_230707', 'Meetz_230710', 'Meetz_230711',
#             'Meetz_230718', 'Meetz_230719', 'Meetz_230720', 'Akshan_240530',
#             'Akshan_240603', 'Akshan_240606', 'Akshan_240607', 'Akshan_240611',
#             'Akshan_240613', 'Akshan_240628', 'Akshan_240703', 'Akshan_240704',
#             'Akshan_240705', 'Akshan_240709', 'Akshan_241014', 'Akshan_241024']
#
# unitList = ['230620_58', '230626_131', '230627_147', '230627_141', '230627_106',
#             '230627_71', '230627_70', '230627_14', '230628_143', '230630_50',
#             '230707_140', '230710_137', '230711_45', '230711_50', '230718_178',
#             '230719_147', '230719_156', '230720_149', '240530_55', '240603_167',
#             '240606_176', '240607_137', '240607_170', '240607_179', '240611_136',
#             '240613_138', '240628_19', '240703_19', '240703_91', '240704_86',
#             '240705_74', '240705_162', '240709_53', '240709_108', '241014_3',
#             '241014_61', '241014_71', '241024_2', '241024_21', '241024_25']

## Combined file list with most of Akshan's sessions run through KS4
# unitList = ['230620_58', '230626_131', '230627_147', '230627_141', '230627_106',
#             '230627_71', '230627_70', '230627_14', '230628_143', '230630_50',
#             '230707_140', '230710_137', '230711_45', '230711_50', '230718_178',
#             '230719_147', '230719_156', '230720_149', '240530_54', '240603_118',
#             '240606_46', '240606_52', '240606_78', '240606_79', '240606_80',
#             '240606_99', '240607_137', '240607_170', '240607_179', '240611_21',
#             '240611_28', '240611_29', '240611_36', '240611_37', '240611_46',
#              '240611_81', '240611_86', '240611_97', '240611_101', '240611_103',
#              '240611_110', '240613_12', '240628_49', '240703_7', '240703_57',
#              '240703_66','240703_72', '240704_24', '240705_56', '240705_60',
#              '240709_53', '240709_108', '241014_3', '241014_61', '241014_71',
#              '241024_2', '241024_21', '241024_25']

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
pnNMIPop = []
npNMIPop = []
ppNMIPop = []
nnNMIPop = []
totNMIPop = []
off1CorrPop = []
off2CorrPop = []
off3CorrPop = []
off4CorrPop = []
totCorrPop = []
singleCorrPop = []
pairedCorrPop = []
sponCorrPop = []
rfWeightPop = []
pairNMIPop = []
rfWeightPairedPop = []
rfWeightSinglePop = []
normalizedRespWeightPop = []
totRFWeightPop = []
offsetDegSepNormPop = []
rfGaborSigmaPop = []
adaptationMat = np.zeros((10000, 4, 100))
adaptationMat[:] = np.nan
adaptC = 0
allBlocksDone = []
masterGoodUnits = []
totUnits = []
corrLists = [[] for _ in range(4)]
popGaborSD = []
popGaborOffsetDeg = []


for file in fileList:
    # Load relevant file here with pyMat reader
    monkeyName, seshDate = file.split('_')

    fileName = f'{monkeyName}_{seshDate}_MTND_Spikes.mat'
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

# make new folder a level above to combine data from both monkeys
    # if not os.path.exists('Normalization'):
    #     os.makedirs('Normalization')
    # os.chdir('Normalization/')

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
            pnNMIPop.append(pCentNMI)
            npNMIPop.append(nCentNMI)
            ppNMIPop.append(ppNMI)
            nnNMIPop.append(nnNMI)
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

    if file == 'Meetz_230719':
        count = np.where(units == 147)[0][0]
        exampleUnitHist = spikeHists[count] * 1000 / blocksDone

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
transectNormalized = transectNormalized - sponNormalized[:, np.newaxis]
transectNormalized = transectNormalized / normVal

pnNMIPop = np.array(pnNMIPop)
npNMIPop = np.array(npNMIPop)
ppNMIPop = np.array(ppNMIPop)
nnNMIPop = np.array(nnNMIPop)

# sort by pref response strength with decreaing response from position 1 to 4.
# for filler in range(1):
#     reordering = []
#     # Loop over each row and get the sorting indices for columns 1 to 4
#     for i in range(prefNormalized.shape[0]):
#         # Get the indices that would sort columns 1 to 4 in descending order
#         sorted_indices = np.argsort(prefNormalized[i, 1:])[::-1]
#         sorted_indices = sorted_indices + 1
#         sorted_indices = np.insert(sorted_indices, 0, 0)
#         reordering.append(sorted_indices)
#
#     reordering = np.array(reordering)
#
#     row_indices = np.arange(prefNormalized.shape[0])[:, None]
#     # Apply reordering to each row at once using advanced indexing
#     prefNormalized1 = prefNormalized[row_indices, reordering]
#     nonprefNormalized1 = nonprefNormalized[row_indices, reordering]
#     pnNormalized1 = pnNormalized[row_indices, (reordering[:, 1:] - 1)]
#     npNormalized1 = npNormalized[row_indices, (reordering[:, 1:] - 1)]
#     ppNormalized1 = ppNormalized[row_indices, (reordering[:, 1:] - 1)]
#     nnNormalized1 = nnNormalized[row_indices, (reordering[:, 1:] - 1)]
#
#     prefMean = np.mean(prefNormalized1, axis=0)
#     nonprefMean = np.mean(nonprefNormalized1, axis=0)
#     pnMean = np.mean(pnNormalized1, axis=0)
#     npMean = np.mean(npNormalized1, axis=0)
#     ppMean = np.mean(ppNormalized1, axis=0)
#     nnMean = np.mean(nnNormalized1, axis=0)
#     prefSEM = np.std(prefNormalized1, axis=0) / np.sqrt(numUnits)
#     nonprefSEM = np.std(nonprefNormalized1, axis=0) / np.sqrt(numUnits)
#     pnSEM = np.std(pnNormalized1, axis=0) / np.sqrt(numUnits)
#     npSEM = np.std(npNormalized1, axis=0) / np.sqrt(numUnits)


# # fit RF weighted normalization to the data
# # R = (Ccent*[Lp, Lnp]) + Cperi*[w1,w2,w3,w4]*[Lp, Lnp]) / (Ccent + [w1,w2,w3,w4]*Cperi + sigma)
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
x = np.arange(0, numSteps + 1)
xNorm = np.arange(1, numSteps + 1)
transX = np.concatenate((np.arange(-numSteps, 0),
                         np.arange(0, numSteps + 1)),
                        axis=0)

# known variables
contrast_center = [1, 0, 0, 0, 0,
                   1, 0, 0, 0, 0,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1]  # First response has contrast in the center, second in periphery
contrast_periphery = [0, 1, 1, 1, 1,
                      0, 1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1]  # First has no periphery, second has periphery with contrast
locations = np.array([-1, 0, 1, 2, 3,
                      -1, 0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3])  # First has no peripheral stimulus, second is at location 1
stim_type_center = np.array([1, -1, -1, -1, -1,
                             0, -1, -1, -1, -1,
                             1, 1, 1, 1,
                             0, 0, 0, 0,
                             1, 1, 1, 1,
                             0, 0, 0, 0])
stim_type_periphery = np.array([-1, 1, 1, 1, 1,
                                -1, 0, 0, 0, 0,
                                0, 0, 0, 0,
                                1, 1, 1, 1,
                                1, 1, 1, 1,
                                0, 0, 0, 0])

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
###################################################### Figure 1c #######################################################
########################################################################################################################

# fit gaussian to average responses
params = gaussFit(transX, transectMean)
xFull = np.linspace(transX[0], transX[-1], 1000)
respFull = gauss(xFull, *params)

# plot transect response
fig, ax1 = plt.subplots(figsize=(8, 6), layout='constrained')
ax1.plot(xFull, respFull, color='black', linestyle='--', label='Gauss Fit')
ax1.plot(transX, transectMean, color='black', label='Pref')
ax1.errorbar(transX, transectMean, yerr=transectSEM, fmt='o', ecolor='black',
             color='black', markersize=7)
ax1.plot(x, nonprefMean, color='grey', label='NO')
ax1.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
             color='grey', markersize=7)

# # individual neuron fits
# for i in transectNormalized:
#     ax1.plot(transX, i, color='purple', alpha=0.2)

ax1.set_ylabel('Normalized Response', fontsize=25, **hfont)
ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
ax1.set_ylim([0, 1.2])
ax1.set_title('transect')

plt.show()

########################################################################################################################
############################################## Example Unit Histogram Figure 2A/B ######################################
########################################################################################################################

# # example unit 230719_147, defined in the main fileList for loop

# # remove spon rate
# exampleUnitHist = exampleUnitHist - exampleUnitHist[0]

# figure
fig, axs = plt.subplots(1, 2, figsize=(14, 6), layout='constrained')
yMax = 0
smoothValue = 5

# P center, N center
pIndx = stimCountIndex[2, 0]
nIndx = stimCountIndex[1, 0]
pHist = exampleUnitHist[pIndx]
nHist = exampleUnitHist[nIndx]
gaussSmoothP = gaussian_filter1d(pHist, smoothValue)
gaussSmoothN = gaussian_filter1d(nHist, smoothValue)
if np.max([gaussSmoothP, gaussSmoothN]) > yMax:
    yMax = np.max([gaussSmoothP, gaussSmoothN])

# PN paired: p, pn loc 1, pn loc 4
for axCount, i in enumerate([1, 2]):
    nOffIndx = stimCountIndex[0, i * 2 - 1]
    pnIndx = stimCountIndex[2, i * 2 - 1]
    nOffHist = exampleUnitHist[nOffIndx]
    pnHist = exampleUnitHist[pnIndx]
    gaussSmoothNOff = gaussian_filter1d(nOffHist, smoothValue)
    gaussSmoothPN = gaussian_filter1d(pnHist, smoothValue)
    gaussSmoothPNPred = (gaussSmoothP + gaussSmoothNOff) / 2
    if np.max([gaussSmoothNOff, gaussSmoothPN]) > yMax:
        yMax = np.max([gaussSmoothNOff, gaussSmoothPN])
    axs[axCount].plot(gaussSmoothP, label='P0', color='black')
    axs[axCount].plot(gaussSmoothNOff, label='N Offset', color='grey')
    axs[axCount].plot(gaussSmoothPN, label='PN', color='green')
    axs[axCount].plot(gaussSmoothPNPred, label='PN Pred', color='green', linestyle='--')
    axs[axCount].set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                             2 * histPrePostMS + trueStimDurMS])
    axs[axCount].set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                                  trueStimDurMS + histPrePostMS], fontsize=20)
    axs[axCount].axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
               color='grey', alpha=0.2)
    axs[axCount].set_xlabel('Time (ms)', fontsize=25)
    axs[axCount].set_ylabel('Firing Rate (spikes/s)', fontsize=25)
    axs[axCount].set_title(f'P Center, N Offset {i}')
    axs[axCount].legend(fontsize=5)

for ax in axs:
    ax.set_ylim([0, yMax * 1.1])

plt.show()

########################################################################################################################
################################################## Bar Plot of Simple Normalization ####################################
########################################################################################################################

# bar plot for individual neurons (PN): Figure 2a
for filler in range(1):
    hfont = {'fontname': 'Arial'}

    # find unit 230626_71 (reference unit)
    j = np.where(masterGoodUnits == '230627_71')[0][0]

    prefCenterBar = prefNormalized[j][0]
    barResponsesPN = {'PN measured': [pnNormalized[j][0], pnNormalized[j][3]],
                      'PN predicted': [(prefNormalized[j][0] + nonprefNormalized[j][1]) / 2,
                                       (prefNormalized[j][0] + nonprefNormalized[j][4]) / 2],
                      'Non-pref offset': [nonprefNormalized[j][1], nonprefNormalized[j][3]]}

    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')

    width = 0.25
    colorList = ['tab:green', 'tab:green', 'grey']
    cluster_positions = [0, 1, 3]
    x_single = [cluster_positions[0]]  # For Position 0
    x_multiple = [cluster_positions[1], cluster_positions[2]]  # For Positions 1 and 4

    ax1.bar(x_single, prefCenterBar, width, label="Pref Center", color="black")
    multiplier = 0
    for attribute, measurement in barResponsesPN.items():
        for i, x_pos in enumerate(x_multiple):
            offset = width * multiplier
            if multiplier == 1:
                rects = ax1.bar(x_pos + offset, measurement[i], width, label=attribute if i == 0 else "",
                                color=colorList[multiplier], hatch='//', alpha=0.5)
            else:
                rects = ax1.bar(x_pos + offset, measurement[i], width, label=attribute if i == 0 else "",
                                color=colorList[multiplier])
        multiplier += 1

    # Calculate the center of each cluster
    tick_positions = [cluster_positions[0],
                      cluster_positions[1] + width,
                      cluster_positions[2] + width]

    # Set custom x-axis tick labels at the center of each cluster
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(['0', '1', '4'], fontsize=20, **hfont)
    ax1.set_xlabel('Receptive Field Location', **hfont, fontsize=25)
    ax1.set_ylim([0, 1.6])
    ax1.set_ylabel('Response Rate (spikes/s)', **hfont, fontsize=25)
    ax1.legend()
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        f'{round(normVal[j][0])}' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4])
    ax1.tick_params(axis='both', width=2, length=8)

    plt.show()
    # # Save as a PDF with tight bounding box and high DPI
    # filePath = f'/Users/chery/Documents/Grad School/Maunsell Lab/Manuscripts/Python Figures/{unitList[j]}.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')


########################################################################################################################
################################# Fig 2c: simple average of representative unit with baseline ##########################
########################################################################################################################

# i = np.where(masterGoodUnits == '230627_71')[0][0]
i = np.where(masterGoodUnits == '230719_147')[0][0]

fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')

mSize = 70
x = np.arange(0, numSteps + 1)
pnFitPred = (prefNormalizedWithBase[i][0] + nonprefNormalizedWithBase[i][1:]) / 2

ax1.scatter(x, prefNormalizedWithBase[i], s=mSize,
            facecolor='white', edgecolor='black')
ax1.plot(x, prefNormalizedWithBase[i], color='black', label='PO', linewidth=1.5)
ax1.scatter(x, nonprefNormalizedWithBase[i], s=mSize,
            facecolor='white', edgecolor='grey')
ax1.plot(x, nonprefNormalizedWithBase[i], color='grey', label='NO', linewidth=1.5)
ax1.plot(xNorm, pnNormalizedWithBase[i], color='green', label='P0 N1', linewidth=1.5)
ax1.scatter(xNorm, pnNormalizedWithBase[i], color='green', s=mSize)
ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
ax1.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=mSize)
ax1.axhline(y=sponNormalizedToPrefCenter[i], color='blue')
ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
ax1.set_ylabel('Response Rate (spikes/s)', fontsize=25, **hfont)
ax1.set_ylim([0, 1.6])
y_ticks = ax1.get_yticks()
# Create a list of y-tick labels with only 0 and 1.0 labeled
y_tick_labels = [
    '0' if tick == 0 else
    f'{round(normVal[i][0])}' if tick == 1.0 else
    ''  # blank for other ticks
    for tick in y_ticks
]
ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
ax1.set_xticks([0, 1, 2, 3, 4])
ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
ax1.set_xlim([-0.5, 4.5])
ax1.tick_params(axis='both', width=2, length=8)

plt.show()


########################################################################################################################
####################################### RF Weighted Normalization with baseline ########################################
########################################################################################################################

# individual neurons
r2ScoresRFWeightFullWithBase = []
for i in range(len(prefNormalized)):

    # population fits
    resp = np.concatenate((prefNormalizedWithBase[i], nonprefNormalizedWithBase[i],
                           pnNormalizedWithBase[i], npNormalizedWithBase[i],
                           ppNormalizedWithBase[i], nnNormalizedWithBase[i],
                           sponNormalizedToPrefCenter[i]), axis=0)
    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                     np.float(sponNormalizedToPrefCenter[i])]  # Lp, Lnp, W1, W2, W3, W4, sigma, b

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

    # # figure
    # fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    # popt_str = (f"Lp = {popt[0]:.2f}\n"
    #             f"Lnp = {popt[1]:.2f}\n"
    #             f"W1 = {popt[2]:.2f}\n"
    #             f"W2 = {popt[3]:.2f}\n"
    #             f"W3 = {popt[4]:.2f}\n"
    #             f"W4 = {popt[5]:.2f}\n"
    #             f"sigma = {popt[6]:.2f}\n"
    #             f"spon = {popt[7]:.2f}\n"
    #             f"r2= {r2:.3f}")
    #
    # mSize = 70
    # x = np.arange(0, numSteps + 1)
    #
    # ax1.scatter(x, prefNormalizedWithBase[i], s=mSize,
    #             facecolor='white', edgecolor='black')
    # ax1.plot(x, prefNormalizedWithBase[i], color='black', label='PO', linewidth=1.5)
    # ax1.scatter(x, nonprefNormalizedWithBase[i], s=mSize,
    #             facecolor='white', edgecolor='grey')
    # ax1.plot(x, nonprefNormalizedWithBase[i], color='grey', label='NO', linewidth=1.5)
    # ax1.plot(xNorm, pnNormalizedWithBase[i], color='green', label='P0 N1', linewidth=1.5)
    # ax1.scatter(xNorm, pnNormalizedWithBase[i], color='green', s=mSize)
    # ax1.plot(xNorm, npNormalizedWithBase[i], color='magenta', label='N0 P1', linewidth=1.5)
    # ax1.scatter(xNorm, npNormalizedWithBase[i], color='magenta', s=mSize)
    # ax1.plot(xNorm, ppNormalizedWithBase[i], color='black', label='P0 P1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, ppNormalizedWithBase[i], color='black', s=mSize)
    # ax1.plot(xNorm, nnNormalizedWithBase[i], color='grey', label='N0 N1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, nnNormalizedWithBase[i], color='grey', s=mSize)
    # ax1.axhline(y=sponNormalizedToPrefCenter[i], color='blue', alpha=1)
    # ax1.axhline(y=y_pred1[-1], linestyle='--', color='blue', alpha=0.8)
    # ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=mSize)
    # ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, npFitPred, facecolors='magenta', edgecolors='m', s=mSize)
    # ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=mSize)
    # ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, nnFitPred, facecolors='gray', edgecolors='gray', s=mSize)
    # ax1.scatter(x, pSingleFitPred, facecolors='white', edgecolors='black', s=mSize)
    # ax1.plot(x, pSingleFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(x, npSingleFitPred, facecolors='white', edgecolors='gray', s=mSize)
    # ax1.plot(x, npSingleFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    # ax1.set_ylabel('Response Rate (spikes/s)', fontsize=25, **hfont)
    # ax1.set_ylim([0, 1.6])
    # ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
    #          bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
    #          transform=ax1.transAxes)
    # y_ticks = ax1.get_yticks()
    # # Create a list of y-tick labels with only 0 and 1.0 labeled
    # y_tick_labels = [
    #     '0' if tick == 0 else
    #     f'{round(normValBase[i][0])}' if tick == 1.0 else
    #     ''  # blank for other ticks
    #     for tick in y_ticks]
    #
    # ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    # ax1.set_xticks([0, 1, 2, 3, 4])
    # ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['bottom'].set_linewidth(2)
    # ax1.spines['left'].set_linewidth(2)
    # ax1.set_xlim([-0.5, 4.5])
    # ax1.tick_params(axis='both', width=2, length=8)
    #
    # # Save as a PDF with tight bounding box and high DPI
    # filePath = f'/Users/chery/Desktop/Project 2 Figs/individualNeuronsFits/{masterGoodUnits[i]}_RFWeightWithBase.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')

r2ScoresRFWeightFullWithBase = np.array(r2ScoresRFWeightFullWithBase)
np.median(r2ScoresRFWeightFullWithBase)

# population
resp = np.concatenate((prefMeanWithBase, nonprefMeanWithBase, pnMeanWithBase,
                       npMeanWithBase, ppMeanWithBase, nnMeanWithBase,
                       [np.mean(sponNormalizedToPrefCenter)]), axis=0)
initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                 np.mean(sponNormalizedToPrefCenter)]  # Lp, Lnp, W1, W2, W3, W4, sigma, b


popt, pcov = curve_fit(model_wrapperWithBase, (contrast_centerSpon,
                                               contrast_peripherySpon,
                                               locationsSpon,
                                               stim_type_centerSpon,
                                               stim_type_peripherySpon),
                       resp,
                       p0=initial_guess)
y_pred1 = apply_fitted_modelWithBase(contrast_centerSpon, contrast_peripherySpon,
                                     locationsSpon, stim_type_centerSpon,
                                     stim_type_peripherySpon, *popt)
r2 = r2_score(resp, y_pred1)
print(r2)
print(popt)

pnFitPred = [fullRFWeightedNormWithBase(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [fullRFWeightedNormWithBase(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullRFWeightedNormWithBase(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullRFWeightedNormWithBase(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullRFWeightedNormWithBase(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullRFWeightedNormWithBase(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullRFWeightedNormWithBase(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullRFWeightedNormWithBase(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# RF weighted with baseline prediction figure
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    popt_str = (f"Lp = {popt[0]:.2f}\n"
                f"Lnp = {popt[1]:.2f}\n"
                f"W1 = {popt[2]:.2f}\n"
                f"W2 = {popt[3]:.2f}\n"
                f"W3 = {popt[4]:.2f}\n"
                f"W4 = {popt[5]:.2f}\n"
                f"sigma = {popt[6]:.2f}")
    ax1.plot(x, prefMeanWithBase, color='black', label='PO', linewidth=1.5)
    ax1.errorbar(x, prefMeanWithBase, yerr=prefWithBaseSEM, fmt='o', ecolor='black',
                 markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax1.plot(x, nonprefMeanWithBase, color='gray', label='NO', linewidth=1.5)
    ax1.errorbar(x, nonprefMeanWithBase, yerr=nonprefWithBaseSEM, fmt='o', ecolor='grey',
                 markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax1.plot(xNorm, pnMeanWithBase, color='green', label='P0 N1', linewidth=1.5)
    ax1.errorbar(xNorm, pnMeanWithBase, yerr=pnWithBaseSEM, fmt='o', ecolor='green',
                 color='green', markersize=7)
    ax1.plot(xNorm, npMeanWithBase, color='magenta', label='N0 P1', linewidth=1.5)
    ax1.errorbar(xNorm, npMeanWithBase, yerr=npWithBaseSEM, fmt='o', ecolor='magenta',
                 color='magenta', markersize=7)
    ax1.plot(xNorm, ppMeanWithBase, color='black', label='P0 P1', linewidth=1.5)
    ax1.errorbar(xNorm, ppMeanWithBase, yerr=ppWithBaseSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax1.plot(xNorm, nnMeanWithBase, color='grey', label='N0 N1', linewidth=1.5)
    ax1.errorbar(xNorm, nnMeanWithBase, yerr=nnWithBaseSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    ax1.axhline(y=np.mean(sponNormalizedToPrefCenter), color='blue')
    ax1.axhline(y=y_pred1[-1], linestyle='--', color='blue', alpha=0.8)
    ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, pnFitPred, facecolors='g', edgecolors='g', s=50)
    ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, npFitPred, facecolors='m', edgecolors='m', s=50)
    ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=50)
    ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, nnFitPred, facecolors='grey', edgecolors='gray', s=50)
    ax1.plot(x, pSingleFitPred, color='black', linestyle='--')
    ax1.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
    ax1.plot(x, npSingleFitPred, color='grey', linestyle='--')
    ax1.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
    ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    ax1.set_ylabel('Normalized Response Rate', fontsize=25, **hfont)
    ax1.set_ylim([-0.05, 1.5])
    # Add a textbox with the popt values to the plot
    ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
             transform=ax1.transAxes)
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.set_xticks([0, 1, 2, 3, 4])
    ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4.5])
    ax1.tick_params(axis='both', width=2, length=8)

    plt.show()


########################################################################################################################
######################################## Classic Heeger with baseline ##################################################
########################################################################################################################
# individual neurons
r2ScoresClassicHeegerWithBase = []
for i in range(len(prefNormalized)):

    resp = np.concatenate((prefNormalizedWithBase[i], nonprefNormalizedWithBase[i],
                           pnNormalizedWithBase[i], npNormalizedWithBase[i],
                           ppNormalizedWithBase[i], nnNormalizedWithBase[i],
                           sponNormalizedToPrefCenter[i]), axis=0)
    initial_guess = [1.0, 0.5, .10,
                     np.float(sponNormalizedToPrefCenter[i])]  # Lp, Lnp, sigma, b

    # classic norm with baseline
    popt, pcov = curve_fit(model_wrapperClassicWithBase, (contrast_centerSpon,
                                                          contrast_peripherySpon,
                                                          locationsSpon,
                                                          stim_type_centerSpon,
                                                          stim_type_peripherySpon), resp,
                           p0=initial_guess)
    y_pred1 = apply_fitted_modelClassicWithBase(contrast_centerSpon, contrast_peripherySpon,
                                                locationsSpon, stim_type_centerSpon,
                                                stim_type_peripherySpon, *popt)
    r2 = r2_score(resp, y_pred1)
    r2ScoresClassicHeegerWithBase.append(r2)

    # single presentation prediction
    pCenterSinglePred = np.array(fullClassicNormWithBase(1, 0, -1, 1, -1, *popt))
    pCenterSinglePred = pCenterSinglePred[np.newaxis]
    pPeriSinglePred = np.array([fullClassicNormWithBase(0, 1, j, -1, 1, *popt) for j in range(4)])
    pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

    npCenterSinglePred = np.array(fullClassicNormWithBase(1, 0, -1, 0, -1, *popt))
    npCenterSinglePred = npCenterSinglePred[np.newaxis]
    npPeriSinglePred = np.array([fullClassicNormWithBase(0, 1, j, -1, 0, *popt) for j in range(4)])
    npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

    # paired presentation prediction
    pnFitPred = [fullClassicNormWithBase(1, 1, j, 1, 0, *popt) for j in range(4)]
    npFitPred = [fullClassicNormWithBase(1, 1, j, 0, 1, *popt) for j in range(4)]
    ppFitPred = [fullClassicNormWithBase(1, 1, j, 1, 1, *popt) for j in range(4)]
    nnFitPred = [fullClassicNormWithBase(1, 1, j, 0, 0, *popt) for j in range(4)]

    predResp = np.concatenate((pSingleFitPred, npSingleFitPred, pnFitPred,
                               npFitPred, ppFitPred, nnFitPred), axis=0)

    # # figure
    # fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    # popt_str = (f"Lp = {popt[0]:.2f}\n"
    #             f"Lnp = {popt[1]:.2f}\n"
    #             f"sigma = {popt[2]:.2f}\n"
    #             f"spon = {popt[3]:.2f}\n"
    #             f"r2= {r2:.3f}")
    #
    # mSize = 70
    # x = np.arange(0, numSteps + 1)
    #
    # ax1.scatter(x, prefNormalizedWithBase[i], s=mSize,
    #             facecolor='white', edgecolor='black')
    # ax1.plot(x, prefNormalizedWithBase[i], color='black', label='PO', linewidth=1.5)
    # ax1.scatter(x, nonprefNormalizedWithBase[i], s=mSize,
    #             facecolor='white', edgecolor='grey')
    # ax1.plot(x, nonprefNormalizedWithBase[i], color='grey', label='NO', linewidth=1.5)
    # ax1.plot(xNorm, pnNormalizedWithBase[i], color='green', label='P0 N1', linewidth=1.5)
    # ax1.scatter(xNorm, pnNormalizedWithBase[i], color='green', s=mSize)
    # ax1.plot(xNorm, npNormalizedWithBase[i], color='magenta', label='N0 P1', linewidth=1.5)
    # ax1.scatter(xNorm, npNormalizedWithBase[i], color='magenta', s=mSize)
    # ax1.plot(xNorm, ppNormalizedWithBase[i], color='black', label='P0 P1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, ppNormalizedWithBase[i], color='black', s=mSize)
    # ax1.plot(xNorm, nnNormalizedWithBase[i], color='grey', label='N0 N1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, nnNormalizedWithBase[i], color='grey', s=mSize)
    # ax1.axhline(y=sponNormalizedToPrefCenter[i], color='blue', alpha=1)
    # ax1.axhline(y=y_pred1[-1], linestyle='--', color='blue', alpha=0.8)
    # ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=mSize)
    # ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, npFitPred, facecolors='magenta', edgecolors='m', s=mSize)
    # ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=mSize)
    # ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, nnFitPred, facecolors='gray', edgecolors='gray', s=mSize)
    # ax1.scatter(x, pSingleFitPred, facecolors='white', edgecolors='black', s=mSize)
    # ax1.plot(x, pSingleFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(x, npSingleFitPred, facecolors='white', edgecolors='gray', s=mSize)
    # ax1.plot(x, npSingleFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    # ax1.set_ylabel('Response Rate (spikes/s)', fontsize=25, **hfont)
    # ax1.set_ylim([0, 1.6])
    # ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
    #          bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
    #          transform=ax1.transAxes)
    # y_ticks = ax1.get_yticks()
    # # Create a list of y-tick labels with only 0 and 1.0 labeled
    # y_tick_labels = [
    #     '0' if tick == 0 else
    #     f'{round(normValBase[i][0])}' if tick == 1.0 else
    #     ''  # blank for other ticks
    #     for tick in y_ticks]
    #
    # ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    # ax1.set_xticks([0, 1, 2, 3, 4])
    # ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['bottom'].set_linewidth(2)
    # ax1.spines['left'].set_linewidth(2)
    # ax1.set_xlim([-0.5, 4.5])
    # ax1.tick_params(axis='both', width=2, length=8)
    #
    # # Save as a PDF with tight bounding box and high DPI
    # filePath = f'/Users/chery/Desktop/Project 2 Figs/individualNeuronsFits/{masterGoodUnits[i]}_classicHeegerWithBase.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')

r2ScoresClassicHeegerWithBase = np.array(r2ScoresClassicHeegerWithBase)
np.median(r2ScoresClassicHeegerWithBase)

# population
resp = np.concatenate((prefMeanWithBase, nonprefMeanWithBase, pnMeanWithBase,
                       npMeanWithBase, ppMeanWithBase, nnMeanWithBase,
                       [np.mean(sponNormalizedToPrefCenter)]), axis=0)
initial_guess = [1.0, 0.5, 0.10,
                 np.mean(sponNormalizedToPrefCenter)]  # Lp, Lnp, sigma, b


popt, pcov = curve_fit(model_wrapperClassicWithBase, (contrast_centerSpon,
                                                      contrast_peripherySpon,
                                                      locationsSpon,
                                                      stim_type_centerSpon,
                                                      stim_type_peripherySpon),
                       resp,
                       p0=initial_guess)
y_pred1 = apply_fitted_modelClassicWithBase(contrast_centerSpon, contrast_peripherySpon,
                                            locationsSpon, stim_type_centerSpon,
                                            stim_type_peripherySpon, *popt)
r2 = r2_score(resp, y_pred1)
print(r2)
print(popt)

pnFitPred = [fullClassicNormWithBase(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [fullClassicNormWithBase(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullClassicNormWithBase(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullClassicNormWithBase(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullClassicNormWithBase(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullClassicNormWithBase(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullClassicNormWithBase(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullClassicNormWithBase(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)


########################################################################################################################
####################################### Weighted Heeger with baseline ##################################################
######################################### weight in numerator only #####################################################
########################################################################################################################

# individual neurons
r2ScoresWeightedHeegerWithBase = []
for i in range(len(prefNormalized)):

    resp = np.concatenate((prefNormalizedWithBase[i], nonprefNormalizedWithBase[i],
                           pnNormalizedWithBase[i], npNormalizedWithBase[i],
                           ppNormalizedWithBase[i], nnNormalizedWithBase[i],
                           sponNormalizedToPrefCenter[i]), axis=0)
    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                     np.float(sponNormalizedToPrefCenter[i])]  # Lp, Lnp, W1, W2, W3, W4, sigma, b

    # weighted norm (numerator only) with baseline
    popt, pcov = curve_fit(model_wrapperWeightedHeegerWithBase, (contrast_centerSpon,
                                                                 contrast_peripherySpon,
                                                                 locationsSpon,
                                                                 stim_type_centerSpon,
                                                                 stim_type_peripherySpon), resp,
                           p0=initial_guess)
    y_pred1 = apply_fitted_modelWeightedHeegerWithBase(contrast_centerSpon, contrast_peripherySpon,
                                                       locationsSpon, stim_type_centerSpon,
                                                       stim_type_peripherySpon, *popt)
    r2 = r2_score(resp, y_pred1)
    r2ScoresWeightedHeegerWithBase.append(r2)

    # single presentation prediction
    pCenterSinglePred = np.array(weightedHeegerNormWithBase(1, 0, -1, 1, -1, *popt))
    pCenterSinglePred = pCenterSinglePred[np.newaxis]
    pPeriSinglePred = np.array([weightedHeegerNormWithBase(0, 1, j, -1, 1, *popt) for j in range(4)])
    pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

    npCenterSinglePred = np.array(weightedHeegerNormWithBase(1, 0, -1, 0, -1, *popt))
    npCenterSinglePred = npCenterSinglePred[np.newaxis]
    npPeriSinglePred = np.array([weightedHeegerNormWithBase(0, 1, j, -1, 0, *popt) for j in range(4)])
    npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

    # paired presentation prediction
    pnFitPred = [weightedHeegerNormWithBase(1, 1, j, 1, 0, *popt) for j in range(4)]
    npFitPred = [weightedHeegerNormWithBase(1, 1, j, 0, 1, *popt) for j in range(4)]
    ppFitPred = [weightedHeegerNormWithBase(1, 1, j, 1, 1, *popt) for j in range(4)]
    nnFitPred = [weightedHeegerNormWithBase(1, 1, j, 0, 0, *popt) for j in range(4)]

    predResp = np.concatenate((pSingleFitPred, npSingleFitPred, pnFitPred,
                               npFitPred, ppFitPred, nnFitPred), axis=0)

    # # figure
    # fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    # popt_str = (f"Lp = {popt[0]:.2f}\n"
    #             f"Lnp = {popt[1]:.2f}\n"
    #             f"W1 = {popt[2]:.2f}\n"
    #             f"W2 = {popt[3]:.2f}\n"
    #             f"W3 = {popt[4]:.2f}\n"
    #             f"W4 = {popt[5]:.2f}\n"
    #             f"sigma = {popt[6]:.2f}\n"
    #             f"spon = {popt[7]:.2f}\n"
    #             f"r2= {r2:.3f}")
    #
    # mSize = 70
    # x = np.arange(0, numSteps + 1)
    #
    # ax1.scatter(x, prefNormalizedWithBase[i], s=mSize,
    #             facecolor='white', edgecolor='black')
    # ax1.plot(x, prefNormalizedWithBase[i], color='black', label='PO', linewidth=1.5)
    # ax1.scatter(x, nonprefNormalizedWithBase[i], s=mSize,
    #             facecolor='white', edgecolor='grey')
    # ax1.plot(x, nonprefNormalizedWithBase[i], color='grey', label='NO', linewidth=1.5)
    # ax1.plot(xNorm, pnNormalizedWithBase[i], color='green', label='P0 N1', linewidth=1.5)
    # ax1.scatter(xNorm, pnNormalizedWithBase[i], color='green', s=mSize)
    # ax1.plot(xNorm, npNormalizedWithBase[i], color='magenta', label='N0 P1', linewidth=1.5)
    # ax1.scatter(xNorm, npNormalizedWithBase[i], color='magenta', s=mSize)
    # ax1.plot(xNorm, ppNormalizedWithBase[i], color='black', label='P0 P1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, ppNormalizedWithBase[i], color='black', s=mSize)
    # ax1.plot(xNorm, nnNormalizedWithBase[i], color='grey', label='N0 N1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, nnNormalizedWithBase[i], color='grey', s=mSize)
    # ax1.axhline(y=sponNormalizedToPrefCenter[i], color='blue', alpha=1)
    # ax1.axhline(y=y_pred1[-1], linestyle='--', color='blue', alpha=0.8)
    # ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=mSize)
    # ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, npFitPred, facecolors='magenta', edgecolors='m', s=mSize)
    # ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=mSize)
    # ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, nnFitPred, facecolors='gray', edgecolors='gray', s=mSize)
    # ax1.scatter(x, pSingleFitPred, facecolors='white', edgecolors='black', s=mSize)
    # ax1.plot(x, pSingleFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(x, npSingleFitPred, facecolors='white', edgecolors='gray', s=mSize)
    # ax1.plot(x, npSingleFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    # ax1.set_ylabel('Response Rate (spikes/s)', fontsize=25, **hfont)
    # ax1.set_ylim([0, 1.6])
    # ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
    #          bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
    #          transform=ax1.transAxes)
    # y_ticks = ax1.get_yticks()
    # # Create a list of y-tick labels with only 0 and 1.0 labeled
    # y_tick_labels = [
    #     '0' if tick == 0 else
    #     f'{round(normValBase[i][0])}' if tick == 1.0 else
    #     ''  # blank for other ticks
    #     for tick in y_ticks]
    #
    # ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    # ax1.set_xticks([0, 1, 2, 3, 4])
    # ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['bottom'].set_linewidth(2)
    # ax1.spines['left'].set_linewidth(2)
    # ax1.set_xlim([-0.5, 4.5])
    # ax1.tick_params(axis='both', width=2, length=8)
    #
    # # Save as a PDF with tight bounding box and high DPI
    # filePath = f'/Users/chery/Desktop/Project 2 Figs/individualNeuronsFits/{masterGoodUnits[i]}_WeightedHeegerWithBase.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')

r2ScoresWeightedHeegerWithBase = np.array(r2ScoresWeightedHeegerWithBase)
np.median(r2ScoresWeightedHeegerWithBase)

# population
resp = np.concatenate((prefMeanWithBase, nonprefMeanWithBase, pnMeanWithBase,
                       npMeanWithBase, ppMeanWithBase, nnMeanWithBase,
                       [np.mean(sponNormalizedToPrefCenter)]), axis=0)
initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                 np.mean(sponNormalizedToPrefCenter)]  # Lp, Lnp, W1, W2, W3, W4, sigma, b


popt, pcov = curve_fit(model_wrapperWeightedHeegerWithBase, (contrast_centerSpon,
                                                             contrast_peripherySpon,
                                                             locationsSpon,
                                                             stim_type_centerSpon,
                                                             stim_type_peripherySpon),
                       resp,
                       p0=initial_guess)
y_pred1 = apply_fitted_modelWeightedHeegerWithBase(contrast_centerSpon, contrast_peripherySpon,
                                                   locationsSpon, stim_type_centerSpon,
                                                   stim_type_peripherySpon, *popt)
r2 = r2_score(resp, y_pred1)
print(r2)
print(popt)

pnFitPred = [weightedHeegerNormWithBase(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [weightedHeegerNormWithBase(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [weightedHeegerNormWithBase(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [weightedHeegerNormWithBase(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(weightedHeegerNormWithBase(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([weightedHeegerNormWithBase(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(weightedHeegerNormWithBase(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([weightedHeegerNormWithBase(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# weighted heeger with baseline prediction figure
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    popt_str = (f"Lp = {popt[0]:.2f}\n"
                f"Lnp = {popt[1]:.2f}\n"
                f"W1 = {popt[2]:.2f}\n"
                f"W2 = {popt[3]:.2f}\n"
                f"W3 = {popt[4]:.2f}\n"
                f"W4 = {popt[5]:.2f}\n"
                f"sigma = {popt[6]:.2f}")
    ax1.plot(x, prefMeanWithBase, color='black', label='PO', linewidth=1.5)
    ax1.errorbar(x, prefMeanWithBase, yerr=prefWithBaseSEM, fmt='o', ecolor='black',
                 markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax1.plot(x, nonprefMeanWithBase, color='gray', label='NO', linewidth=1.5)
    ax1.errorbar(x, nonprefMeanWithBase, yerr=nonprefWithBaseSEM, fmt='o', ecolor='grey',
                 markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax1.plot(xNorm, pnMeanWithBase, color='green', label='P0 N1', linewidth=1.5)
    ax1.errorbar(xNorm, pnMeanWithBase, yerr=pnWithBaseSEM, fmt='o', ecolor='green',
                 color='green', markersize=7)
    ax1.plot(xNorm, npMeanWithBase, color='magenta', label='N0 P1', linewidth=1.5)
    ax1.errorbar(xNorm, npMeanWithBase, yerr=npWithBaseSEM, fmt='o', ecolor='magenta',
                 color='magenta', markersize=7)
    ax1.plot(xNorm, ppMeanWithBase, color='black', label='P0 P1', linewidth=1.5)
    ax1.errorbar(xNorm, ppMeanWithBase, yerr=ppWithBaseSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax1.plot(xNorm, nnMeanWithBase, color='grey', label='N0 N1', linewidth=1.5)
    ax1.errorbar(xNorm, nnMeanWithBase, yerr=nnWithBaseSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    ax1.axhline(y=np.mean(sponNormalizedToPrefCenter), color='blue')
    ax1.axhline(y=y_pred1[-1], linestyle='--', color='blue', alpha=0.8)
    ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, pnFitPred, facecolors='g', edgecolors='g', s=50)
    ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, npFitPred, facecolors='m', edgecolors='m', s=50)
    ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=50)
    ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, nnFitPred, facecolors='grey', edgecolors='gray', s=50)
    ax1.plot(x, pSingleFitPred, color='black', linestyle='--')
    ax1.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
    ax1.plot(x, npSingleFitPred, color='grey', linestyle='--')
    ax1.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
    ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    ax1.set_ylabel('Normalized Response Rate', fontsize=25, **hfont)
    ax1.set_ylim([-0.05, 1.5])
    # Add a textbox with the popt values to the plot
    ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
             transform=ax1.transAxes)
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.set_xticks([0, 1, 2, 3, 4])
    ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4.5])
    ax1.tick_params(axis='both', width=2, length=8)

    plt.show()


########################################################################################################################
###################################################### Work Space 2 ####################################################
################################################# Explainable Variance #################################################
################################# In this version I simulate the raw rates with the baseline ###########################
########################################################################################################################
explainable_variances = []
goodI = []
for i in range(len(prefNormalized)):
    # Skip if any entry is zero-dimensional or empty
    if (prefNormalized[i].ndim == 0 or nonprefNormalized[i].ndim == 0 or
        pnNormalized[i].ndim == 0 or npNormalized[i].ndim == 0 or
        ppNormalized[i].ndim == 0 or nnNormalized[i].ndim == 0):
        continue

    # Convert normalized rates to raw rates
    raw_spon = sponNormalizedToPrefCenter[i] * normVal[i]
    # add spon rate back in
    raw_pref = (prefNormalized[i] * normVal[i]) + raw_spon
    raw_nonpref = (nonprefNormalized[i] * normVal[i]) + raw_spon
    raw_pn = (pnNormalized[i] * normVal[i]) + raw_spon
    raw_np = (npNormalized[i] * normVal[i]) + raw_spon
    raw_pp = (ppNormalized[i] * normVal[i]) + raw_spon
    raw_nn = (nnNormalized[i] * normVal[i]) + raw_spon

    # Concatenate into a single response vector
    raw_resp = np.concatenate((raw_pref, raw_nonpref, raw_pn, raw_np, raw_pp, raw_nn, raw_spon), axis=0)

    # # Fit model to raw responses
    # initial_guess = [1.0 * float(normVal[i]), 0.5 * float(normVal[i]),
    #                  1.0, 0.75, 0.50, 0.10, 0.10, raw_spon]
    initial_guess = [raw_pref[0], raw_nonpref[0], 1.0, 0.75, 0.50, 0.10, 0.10,
                     np.float(raw_spon)]

    try:
        popt, _ = curve_fit(model_wrapperWithBase,
                            (contrast_centerSpon, contrast_peripherySpon, locationsSpon,
                             stim_type_centerSpon, stim_type_peripherySpon),
                            raw_resp, p0=initial_guess)
    except RuntimeError:
        continue  # Skip neuron if fit fails

    # Predicted rates in spikes/s
    predicted_rates = np.array(apply_fitted_modelWithBase(contrast_centerSpon,
                                                          contrast_peripherySpon,
                                                          locationsSpon,
                                                          stim_type_centerSpon,
                                                          stim_type_peripherySpon,
                                                          *popt))

    # Simulate trials using Poisson spike counts (convert to spikes per 250 ms)
    n_repeats = allBlocksDone[i]
    rate_per_trial = predicted_rates * 0.25

    # Ensure all λ values are valid for Poisson
    rate_per_trial = np.nan_to_num(rate_per_trial, nan=0.0, posinf=0.0, neginf=0.0)
    rate_per_trial[rate_per_trial < 0] = 0.0

    simulated_trials = np.random.poisson(rate_per_trial[:, None], size=(len(rate_per_trial), n_repeats))

    # Average simulated responses
    simulated_means = simulated_trials.mean(axis=1)

    # # # Convert back to spikes/s
    simulated_means = simulated_means / 0.25

    # Refit model to simulated averages
    try:
        popt_sim, _ = curve_fit(model_wrapperWithBase,
                                (contrast_centerSpon, contrast_peripherySpon, locationsSpon,
                                 stim_type_centerSpon, stim_type_peripherySpon),
                                simulated_means, p0=initial_guess)
    except RuntimeError:
        continue

    refit_preds = np.array(apply_fitted_modelWithBase(contrast_centerSpon, contrast_peripherySpon,
                                                      locationsSpon, stim_type_centerSpon,
                                                      stim_type_peripherySpon, *popt_sim))

    # Compute explainable variance (R²)
    ss_resid = np.sum((simulated_means - refit_preds) ** 2)
    ss_total = np.sum((simulated_means - np.mean(simulated_means)) ** 2)
    R2 = 1 - (ss_resid / ss_total)
    explainable_variances.append(R2)
    goodI.append(i)

explainable_variances = np.array(explainable_variances)
# Population-level explainable variance
population_explainable_variance = np.median(explainable_variances)
print("Population-level explainable variance:", population_explainable_variance)

# Scatter plot
plt.figure(figsize=(5, 5))
plt.scatter(explainable_variances, r2ScoresRFWeightFullWithBase[goodI], color='blue', alpha=0.7)

minVal = 0.7
maxVal = 1.03
# Axis limits
plt.xlim(minVal, maxVal)
plt.ylim(minVal, maxVal)

# Unity line
plt.plot([minVal, maxVal], [minVal, maxVal], 'k--')

# Labels and title
plt.xlabel("Explainable variance (R² from simulated data)")
plt.ylabel("Model fit R² (original data)")
plt.title("Explainable variance vs. Original model R²")
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Compute ratio
ratio = r2ScoresRFWeightFullWithBase / explainable_variances

# Plot using matplotlib
plt.figure(figsize=(8, 5))
plt.hist(ratio, bins=20, color="royalblue", edgecolor="black")
plt.axvline(1.0, color='red', linestyle='--', linewidth=2, label='Unity (1.0)')
plt.xlabel("Model R² / Explainable Variance", fontsize=12)
plt.ylabel("Count", fontsize=12)
plt.title("Distribution of Model R² as a Fraction of Explainable Variance", fontsize=14)
plt.legend()
plt.tight_layout()
plt.show()

########################################################################################################################
############# Evaluate Weighted Normalization vs Weighted Normalization Heuristic Model with baseline ##################
########################################################################################################################

# individual neurons
r2ScoresRFWeightHeuristic = []
for i in range(len(prefNormalized)):

    resp = np.concatenate((prefNormalizedWithBase[i], nonprefNormalizedWithBase[i],
                           pnNormalizedWithBase[i], npNormalizedWithBase[i],
                           ppNormalizedWithBase[i], nnNormalizedWithBase[i],
                           sponNormalizedToPrefCenter[i]), axis=0)
    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                     np.float(sponNormalizedToPrefCenter[i])]

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

r2ScoresRFWeightHeuristic = np.array(r2ScoresRFWeightHeuristic)
np.median(r2ScoresRFWeightHeuristic)

# population
# weighted norm with baseline (heuristic model)
resp = np.concatenate((prefMeanWithBase, nonprefMeanWithBase, pnMeanWithBase,
                       npMeanWithBase, ppMeanWithBase, nnMeanWithBase,
                       [np.mean(sponNormalizedToPrefCenter)]), axis=0)
initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                 np.mean(sponNormalizedToPrefCenter)]  # Lp, Lnp, W1, W2, W3, W4, sigma, b

popt, pcov = curve_fit(model_wrapperHeuristic, (contrast_centerSpon, contrast_peripherySpon,
                                                locationsSpon, stim_type_centerSpon,
                                                stim_type_peripherySpon), resp,
                       p0=initial_guess)
y_pred1 = apply_fitted_modelHeuristic(contrast_centerSpon, contrast_peripherySpon,
                                      locationsSpon, stim_type_centerSpon,
                                      stim_type_peripherySpon, *popt)
r2 = r2_score(resp, y_pred1)
print(r2)
print(popt)

# scatter plot of individual neuron fit r2 scores for heuristic weighted vs rf weighted
for filler in range(1):
    fix, ax1 = plt.subplots()
    ax1.scatter(r2ScoresRFWeightFullWithBase, r2ScoresRFWeightHeuristic)
    ax1.set_xlabel('r2 Scores weighted normalization')
    ax1.set_ylabel('r2 Scores heuristic weighted norm')

    x_min, x_max = ax1.get_xlim()
    y_min, y_max = ax1.get_ylim()
    min_val = min(x_min, y_min)
    max_val = max(x_max, y_max)
    ax1.plot([min_val, max_val], [min_val, max_val], 'k--', label='Unity Line', zorder=0)

    # Reset limits to ensure the unity line fits perfectly
    ax1.set_xlim([min_val, max_val])
    ax1.set_ylim([min_val, max_val])
    stat, p_value = wilcoxon(r2ScoresRFWeightFullWithBase, r2ScoresRFWeightHeuristic, alternative='two-sided')
    print(p_value)
    plt.show()


########################################################################################################################
########################################################################################################################
########################################################################################################################


########################################################################################################################
####################################### The analyses below use models without the baseline term ########################
########################################################################################################################


########################################################################################################################
########################################################################################################################
########################################################################################################################


########################################################################################################################
####################################### Fig 2c: simple average of representative unit without base #####################
########################################################################################################################

# i = np.where(masterGoodUnits == '230627_71')[0][0]
i = np.where(masterGoodUnits == '230719_147')[0][0]

fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')

mSize = 70
x = np.arange(0, numSteps + 1)
pnFitPred = (prefNormalized[i][0] + nonprefNormalized[i][1:]) / 2

ax1.scatter(x, prefNormalized[i], s=mSize,
            facecolor='white', edgecolor='black')
ax1.plot(x, prefNormalized[i], color='black', label='PO', linewidth=1.5)
ax1.scatter(x, nonprefNormalized[i], s=mSize,
            facecolor='white', edgecolor='grey')
ax1.plot(x, nonprefNormalized[i], color='grey', label='NO', linewidth=1.5)
ax1.plot(xNorm, pnNormalized[i], color='green', label='P0 N1', linewidth=1.5)
ax1.scatter(xNorm, pnNormalized[i], color='green', s=mSize)
ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
ax1.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=mSize)
ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
ax1.set_ylabel('Response Rate (spikes/s)', fontsize=25, **hfont)
ax1.set_ylim([0, 1.6])
y_ticks = ax1.get_yticks()
# Create a list of y-tick labels with only 0 and 1.0 labeled
y_tick_labels = [
    '0' if tick == 0 else
    f'{round(normVal[i][0])}' if tick == 1.0 else
    ''  # blank for other ticks
    for tick in y_ticks
]
ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
ax1.set_xticks([0, 1, 2, 3, 4])
ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
ax1.set_xlim([-0.5, 4.5])
ax1.tick_params(axis='both', width=2, length=8)

plt.show()

########################################################################################################################
###################################### Classic Heeger Normalization Weighted Numerator #################################
########################################################################################################################

# individual neurons fit: Figure 2d
fullClassicNormPred = []
measuredResp = []
r2ScoresClassicNormFull = []
for i in range(len(prefNormalized)):
    resp = np.concatenate((prefNormalized[i], nonprefNormalized[i], pnNormalized[i],
                           npNormalized[i], ppNormalized[i], nnNormalized[i]), axis=0)

    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10]  # Lp, Lnp, W1, W2, W3, W4, sigma

    popt, pcov = curve_fit(model_wrapperClassicWeightedNorm, (contrast_center, contrast_periphery, locations, stim_type_center,
                                                              stim_type_periphery),
                           resp, p0=initial_guess, maxfev=10000000)
    y_pred1 = apply_fitted_modelClassicWeighted(contrast_center, contrast_periphery, locations,
                                                stim_type_center, stim_type_periphery, *popt)
    r2 = r2_score(resp, y_pred1)
    r2ScoresClassicNormFull.append(r2)

    # single presentation prediction
    pCenterSinglePred = np.array(fullClassicWeightedNorm(1, 0, -1, 1, -1, *popt))
    pCenterSinglePred = pCenterSinglePred[np.newaxis]
    pPeriSinglePred = np.array([fullClassicWeightedNorm(0, 1, j, -1, 1, *popt) for j in range(4)])
    pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

    npCenterSinglePred = np.array(fullClassicWeightedNorm(1, 0, -1, 0, -1, *popt))
    npCenterSinglePred = npCenterSinglePred[np.newaxis]
    npPeriSinglePred = np.array([fullClassicWeightedNorm(0, 1, j, -1, 0, *popt) for j in range(4)])
    npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

    # paired presentation prediction
    pnFitPred = [fullClassicWeightedNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
    npFitPred = [fullClassicWeightedNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
    ppFitPred = [fullClassicWeightedNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
    nnFitPred = [fullClassicWeightedNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

    predResp = np.concatenate((pSingleFitPred, npSingleFitPred, pnFitPred,
                               npFitPred, ppFitPred, nnFitPred), axis=0)

    for j in range(len(resp)):
        fullClassicNormPred.append(predResp[j])
        measuredResp.append(resp[j])

    # fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    # popt_str = (f"Lp = {popt[0]:.2f}\n"
    #             f"Lnp = {popt[1]:.2f}\n"
    #             f"W1 = {popt[2]:.2f}\n"
    #             f"W2 = {popt[3]:.2f}\n"
    #             f"W3 = {popt[4]:.2f}\n"
    #             f"W4 = {popt[5]:.2f}\n"
    #             f"sigma = {popt[6]:.2f}\n"
    #             f"r2={r2:.3f}")
    # mSize = 70
    # x = np.arange(0, numSteps + 1)
    #
    # ax1.plot(x, prefNormalized[i], color='black', label='PO', linewidth=1.5)
    # ax1.scatter(x, prefNormalized[i], s=mSize,
    #             facecolor='white', edgecolor='black')
    # ax1.plot(x, nonprefNormalized[i], color='grey', label='NO', linewidth=1.5)
    # ax1.scatter(x, nonprefNormalized[i], color='grey', s=mSize,
    #             facecolor='white', edgecolor='gray')
    # ax1.plot(xNorm, pnNormalized[i], color='green', label='P0 N1', linewidth=1.5)
    # ax1.scatter(xNorm, pnNormalized[i], color='green', s=mSize)
    # ax1.plot(xNorm, npNormalized[i], color='magenta', label='N0 P1', linewidth=1.5)
    # ax1.scatter(xNorm, npNormalized[i], color='magenta', s=mSize)
    # ax1.plot(xNorm, ppNormalized[i], color='black', label='P0 P1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, ppNormalized[i], color='black', s=mSize)
    # ax1.plot(xNorm, nnNormalized[i], color='grey', label='N0 N1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, nnNormalized[i], color='grey', s=mSize)
    # ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=mSize)
    # ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, npFitPred, facecolors='magenta', edgecolors='m', s=mSize)
    # ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=mSize)
    # ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, nnFitPred, facecolors='gray', edgecolors='gray', s=mSize)
    # ax1.plot(x, pSingleFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(x, pSingleFitPred, facecolors='white', edgecolors='black', s=mSize)
    # ax1.plot(x, npSingleFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.scatter(x, npSingleFitPred, facecolors='white', edgecolors='gray', s=mSize)
    # ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    # ax1.set_ylabel('Response Rate (spikes/s)', fontsize=25, **hfont)
    # ax1.set_ylim([0, 1.6])
    # ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
    #          bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
    #          transform=ax1.transAxes)
    # y_ticks = ax1.get_yticks()
    # # Create a list of y-tick labels with only 0 and 1.0 labeled
    # y_tick_labels = [
    #     '0' if tick == 0 else
    #     f'{round(normVal[i][0])}' if tick == 1.0 else
    #     ''  # blank for other ticks
    #     for tick in y_ticks
    # ]
    # ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    # ax1.set_xticks([0, 1, 2, 3, 4])
    # ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['bottom'].set_linewidth(2)
    # ax1.spines['left'].set_linewidth(2)
    # ax1.set_xlim([-0.5, 4.5])
    # ax1.tick_params(axis='both', width=2, length=8)
    # ax1.tick_params(axis='both', width=2, length=8)
    #
    # # Save the figure to the specified path
    # filePath = f'/Users/chery/Desktop/Project 2 Figs/individualNeuronsFits/{masterGoodUnits[i]}_weightedHeeger.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')

# population fits
resp = np.concatenate((prefMean, nonprefMean, pnMean, npMean, ppMean, nnMean), axis=0)
initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10]  # Lp, Lnp, W1, W2, W3, W4, sigma
popt, pcov = curve_fit(model_wrapperClassicWeightedNorm, (contrast_center, contrast_periphery,
                                                          locations, stim_type_center,
                                                          stim_type_periphery),
                       resp, p0=initial_guess)

y_pred1 = apply_fitted_modelClassicWeighted(contrast_center, contrast_periphery, locations,
                                            stim_type_center, stim_type_periphery, *popt)
r2 = r2_score(resp, y_pred1)

pnFitPred = [fullClassicWeightedNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [fullClassicWeightedNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullClassicWeightedNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullClassicWeightedNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullClassicWeightedNorm(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullClassicWeightedNorm(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullClassicWeightedNorm(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullClassicWeightedNorm(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# main heeger weighed plot: figure 3a
hfont = {'fontname': 'Arial'}
for filler in range(1):
    x = np.arange(0, numSteps + 1)
    fig, ax2 = plt.subplots(figsize=(7, 7), layout='constrained')

    ax2.plot(x, prefMean, color='black', label='PO', linewidth=1.5)
    # ax2.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
    #              markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax2.errorbar(x, prefMean, yerr=0, fmt='o', ecolor='black', markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax2.plot(x, nonprefMean, color='grey', label='NO', linewidth=1.5)
    # ax2.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
    #              markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax2.errorbar(x, nonprefMean, yerr=0, fmt='o', ecolor='grey', markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax2.plot(xNorm, pnMean, color='green', label='P0 N1', linewidth=1.5)
    # ax2.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
    #              color='green', markersize=7)
    ax2.errorbar(xNorm, pnMean, yerr=0, fmt='o', ecolor='green', color='green', markersize=7)
    ax2.plot(xNorm, npMean, color='magenta', label='N0 P1', linewidth=1.5)
    # ax2.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='magenta',
    #              color='magenta', markersize=7)
    ax2.errorbar(xNorm, npMean, yerr=0, fmt='o', ecolor='magenta', color='magenta', markersize=7)
    ax2.plot(xNorm, ppMean, color='black', label='P0 P1', linewidth=1.5)
    # ax2.errorbar(xNorm, ppMean, yerr=ppSEM, fmt='o', ecolor='black',
    #              color='black', markersize=7)
    ax2.errorbar(xNorm, ppMean, yerr=0, fmt='o', ecolor='black', color='black', markersize=7)
    ax2.plot(xNorm, nnMean, color='grey', label='N0 N1', linewidth=1.5)
    # ax2.errorbar(xNorm, nnMean, yerr=nnSEM, fmt='o', ecolor='grey',
    #              color='grey', markersize=7)
    ax2.errorbar(xNorm, nnMean, yerr=0, fmt='o', ecolor='grey', color='grey', markersize=7)
    ax2.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    ax2.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=50)
    ax2.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    ax2.scatter(xNorm, npFitPred, facecolors='magenta', edgecolors='m', s=50)
    ax2.plot(xNorm, ppFitPred, color='black', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax2.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=50)
    ax2.plot(xNorm, nnFitPred, color='grey', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax2.scatter(xNorm, nnFitPred, facecolors='gray', edgecolors='gray', s=50)
    ax2.plot(x, pSingleFitPred, color='black', linestyle='--', linewidth=1.5)
    ax2.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
    ax2.plot(x, npSingleFitPred, color='grey', linestyle='--', linewidth=1.5)
    ax2.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
    ax2.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    ax2.set_ylabel('Normalized Response Rate', fontsize=25, **hfont)
    ax2.set_ylim([0, 1.50])
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_ticks = ax2.get_yticks()
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax2.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax2.set_xticks([0, 1, 2, 3, 4])
    ax2.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['left'].set_linewidth(2)
    ax2.set_xlim([-0.5, 4.5])
    ax2.tick_params(axis='both', width=2, length=8)

    plt.show()
    # # Save the figure to the specified path
    # filePath = f'/Users/chery/Documents/Grad School/Maunsell Lab/Manuscripts/Python Figures/populationWeightedHeeger.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')

########################################################################################################################
########################################### Receptive Field Weighted Normalization #####################################
########################################################################################################################

# individual neurons fit: Figure 2e
fullRFWeightPred = []
measuredResp = []
r2ScoresRFWeightFull = []
modelPredictedVariance = []
fitSigmaVals = []
for i in range(len(prefNormalized)):
    resp = np.concatenate((prefNormalized[i], nonprefNormalized[i], pnNormalized[i],
                           npNormalized[i], ppNormalized[i], nnNormalized[i]), axis=0)

    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10]  # Lp, Lnp, W1, W2, W3, W4, sigma

    popt, pcov = curve_fit(model_wrapper, (contrast_center, contrast_periphery, locations, stim_type_center,
                                           stim_type_periphery), resp, p0=initial_guess)
    y_pred1 = apply_fitted_model(contrast_center, contrast_periphery, locations,
                                 stim_type_center, stim_type_periphery, *popt)
    r2 = r2_score(resp, y_pred1)
    r2ScoresRFWeightFull.append(r2)

    # single presentation prediction
    pCenterSinglePred = np.array(fullRFWeightedNorm(1, 0, -1, 1, -1, *popt))
    pCenterSinglePred = pCenterSinglePred[np.newaxis]
    pPeriSinglePred = np.array([fullRFWeightedNorm(0, 1, j, -1, 1, *popt) for j in range(4)])
    pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

    npCenterSinglePred = np.array(fullRFWeightedNorm(1, 0, -1, 0, -1, *popt))
    npCenterSinglePred = npCenterSinglePred[np.newaxis]
    npPeriSinglePred = np.array([fullRFWeightedNorm(0, 1, j, -1, 0, *popt) for j in range(4)])
    npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

    # paired presentation prediction
    pnFitPred = [fullRFWeightedNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
    npFitPred = [fullRFWeightedNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
    ppFitPred = [fullRFWeightedNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
    nnFitPred = [fullRFWeightedNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

    predResp = np.concatenate((pSingleFitPred, npSingleFitPred, pnFitPred,
                               npFitPred, ppFitPred, nnFitPred), axis=0)

    for j in range(len(resp)):
        fullRFWeightPred.append(predResp[j])
        measuredResp.append(resp[j])

    modelPredictedVariance.append(np.var(predResp))
    fitSigmaVals.append(popt[6])

    # fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    # popt_str = (f"Lp = {popt[0]:.2f}\n"
    #             f"Lnp = {popt[1]:.2f}\n"
    #             f"W1 = {popt[2]:.2f}\n"
    #             f"W2 = {popt[3]:.2f}\n"
    #             f"W3 = {popt[4]:.2f}\n"
    #             f"W4 = {popt[5]:.2f}\n"
    #             f"sigma = {popt[6]:.2f}\n"
    #             f"r2={r2:.3f}")
    # mSize = 70
    # x = np.arange(0, numSteps + 1)
    #
    # ax1.scatter(x, prefNormalized[i], s=mSize,
    #             facecolor='white', edgecolor='black')
    # ax1.plot(x, prefNormalized[i], color='black', label='PO', linewidth=1.5)
    # ax1.scatter(x, nonprefNormalized[i], s=mSize,
    #             facecolor='white', edgecolor='grey')
    # ax1.plot(x, nonprefNormalized[i], color='grey', label='NO', linewidth=1.5)
    # ax1.plot(xNorm, pnNormalized[i], color='green', label='P0 N1', linewidth=1.5)
    # ax1.scatter(xNorm, pnNormalized[i], color='green', s=mSize)
    # ax1.plot(xNorm, npNormalized[i], color='magenta', label='N0 P1', linewidth=1.5)
    # ax1.scatter(xNorm, npNormalized[i], color='magenta', s=mSize)
    # ax1.plot(xNorm, ppNormalized[i], color='black', label='P0 P1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, ppNormalized[i], color='black', s=mSize)
    # ax1.plot(xNorm, nnNormalized[i], color='grey', label='N0 N1',
    #          linewidth=1.5)
    # ax1.scatter(xNorm, nnNormalized[i], color='grey', s=mSize)
    # ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=mSize)
    # ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, npFitPred, facecolors='magenta', edgecolors='m', s=mSize)
    # ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=mSize)
    # ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.scatter(xNorm, nnFitPred, facecolors='gray', edgecolors='gray', s=mSize)
    # ax1.scatter(x, pSingleFitPred, facecolors='white', edgecolors='black', s=mSize)
    # ax1.plot(x, pSingleFitPred, color='black', linestyle='--', linewidth=1.5)
    # ax1.scatter(x, npSingleFitPred, facecolors='white', edgecolors='gray', s=mSize)
    # ax1.plot(x, npSingleFitPred, color='grey', linestyle='--', linewidth=1.5)
    # ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    # ax1.set_ylabel('Response Rate (spikes/s)', fontsize=25, **hfont)
    # ax1.set_ylim([0, 1.6])
    # ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
    #          bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
    #          transform=ax1.transAxes)
    # y_ticks = ax1.get_yticks()
    # # Create a list of y-tick labels with only 0 and 1.0 labeled
    # y_tick_labels = [
    #     '0' if tick == 0 else
    #     f'{round(normVal[i][0])}' if tick == 1.0 else
    #     ''  # blank for other ticks
    #     for tick in y_ticks
    # ]
    # ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    # ax1.set_xticks([0, 1, 2, 3, 4])
    # ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['bottom'].set_linewidth(2)
    # ax1.spines['left'].set_linewidth(2)
    # ax1.set_xlim([-0.5, 4.5])
    # ax1.tick_params(axis='both', width=2, length=8)
    #
    # # Save as a PDF with tight bounding box and high DPI
    # filePath = f'/Users/chery/Desktop/Project 2 Figs/individualNeuronsFits/{masterGoodUnits[i]}.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')

r2ScoresRFWeightFull = np.array(r2ScoresRFWeightFull)

# population fits: Figure 3b
resp = np.concatenate((prefMean, nonprefMean, pnMean, npMean, ppMean, nnMean), axis=0)
initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10]  # Lp, Lnp, W1, W2, W3, W4, sigma
popt, pcov = curve_fit(model_wrapper, (contrast_center, contrast_periphery, locations, stim_type_center,
                                       stim_type_periphery), resp,
                       p0=initial_guess)


y_pred1 = apply_fitted_model(contrast_center, contrast_periphery, locations,
                             stim_type_center, stim_type_periphery, *popt)
r2 = r2_score(resp, y_pred1)

pnFitPred = [fullRFWeightedNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [fullRFWeightedNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullRFWeightedNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullRFWeightedNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullRFWeightedNorm(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullRFWeightedNorm(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullRFWeightedNorm(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullRFWeightedNorm(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# simple average prediction
pnAvgPred = (prefMean[0] + nonprefMean[1:]) / 2

#### Figure 3b main RF weighted prediction figure
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    popt_str = (f"Lp = {popt[0]:.2f}\n"
                f"Lnp = {popt[1]:.2f}\n"
                f"W1 = {popt[2]:.2f}\n"
                f"W2 = {popt[3]:.2f}\n"
                f"W3 = {popt[4]:.2f}\n"
                f"W4 = {popt[5]:.2f}\n"
                f"sigma = {popt[6]:.2f}")
    ax1.plot(x, prefMean, color='black', label='PO', linewidth=1.5)
    # ax1.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
    #              markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax1.errorbar(x, prefMean, yerr=0, fmt='o', ecolor='black', markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax1.plot(x, nonprefMean, color='gray', label='NO', linewidth=1.5)
    # ax1.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
    #              markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax1.errorbar(x, nonprefMean, yerr=0, fmt='o', ecolor='grey', markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax1.plot(xNorm, pnMean, color='green', label='P0 N1', linewidth=1.5)
    # ax1.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
    #              color='green', markersize=7)
    ax1.errorbar(xNorm, pnMean, yerr=0, fmt='o', ecolor='green', color='green', markersize=7)
    ax1.plot(xNorm, npMean, color='magenta', label='N0 P1', linewidth=1.5)
    # ax1.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='magenta',
    #              color='magenta', markersize=7)
    ax1.errorbar(xNorm, npMean, yerr=0, fmt='o', ecolor='magenta', color='magenta', markersize=7)
    ax1.plot(xNorm, ppMean, color='black', label='P0 P1', linewidth=1.5)
    # ax1.errorbar(xNorm, ppMean, yerr=ppSEM, fmt='o', ecolor='black',
    #              color='black', markersize=7)
    ax1.errorbar(xNorm, ppMean, yerr=0, fmt='o', ecolor='black', color='black', markersize=7)
    ax1.plot(xNorm, nnMean, color='grey', label='N0 N1', linewidth=1.5)
    # ax1.errorbar(xNorm, nnMean, yerr=nnSEM, fmt='o', ecolor='grey',
    #              color='grey', markersize=7)
    ax1.errorbar(xNorm, nnMean, yerr=0, fmt='o', ecolor='grey', color='grey', markersize=7)
    ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, pnFitPred, facecolors='g', edgecolors='g', s=50)
    ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, npFitPred, facecolors='m', edgecolors='m', s=50)
    ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=50)
    ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, nnFitPred, facecolors='grey', edgecolors='gray', s=50)
    ax1.plot(x, pSingleFitPred, color='black', linestyle='--')
    ax1.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
    ax1.plot(x, npSingleFitPred, color='grey', linestyle='--')
    ax1.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
    ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    ax1.set_ylabel('Normalized Response Rate', fontsize=25, **hfont)
    ax1.set_ylim([0, 1.5])
    # Add a textbox with the popt values to the plot
    # ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
    #          bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
    #          transform=ax1.transAxes)
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.set_xticks([0, 1, 2, 3, 4])
    ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4.5])
    ax1.tick_params(axis='both', width=2, length=8)

    plt.show()

    # # Save the figure to the specified path
    # filePath = f'/Users/chery/Documents/Grad School/Maunsell Lab/Manuscripts/Python Figures/populationRFWeighted.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')


#### Figure 3c scatter comparing R2 for weighted heeger vs RF weighted norm
for filler in range(1):
    fig, ax2 = plt.subplots(figsize=(7, 7), layout='constrained')
    ax2.scatter(r2ScoresClassicNormFull, r2ScoresRFWeightFull, color='black',
                facecolors='none', s=100)
    ax2.set_ylim([0.70, 1])
    ax2.set_xlim([0.70, 1])
    ax2.set_xlabel(r'$r^2$ Scores Weighted Classic', fontsize=25, **hfont)
    ax2.set_ylabel(r'$r^2$ Scores RF Weighted', fontsize=25, **hfont)
    ax2.spines['top'].set_linewidth(2)
    ax2.spines['right'].set_linewidth(2)
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['left'].set_linewidth(2)
    y_ticks = ax2.get_yticks()
    tolerance = 0.01
    y_tick_labels = [
        '0.70' if abs(tick - 0.7) < tolerance else
        '1.00' if abs(tick - 1.0) < tolerance else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax2.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax2.set_xticklabels(y_tick_labels, fontsize=20, **hfont)

    line = lines.Line2D([0, 1], [0, 1], color='black')
    transform = ax2.transAxes
    line.set_transform(transform)
    ax2.add_line(line)

    # stat test to compare models and to see if RF weight does better than classic weight
    stat, p_value = wilcoxon(r2ScoresRFWeightFull, r2ScoresClassicNormFull, alternative='greater')
    ax2.text(0.02, 0.98, f'Wilcoxon p-value: {p_value:.2e}',
             transform=ax2.transAxes, fontsize=15, verticalalignment='top',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'))

    plt.show()
    # # Save the figure to the specified path
    # filePath = f'/Users/chery/Documents/Grad School/Maunsell Lab/Manuscripts/Python Figures/r2ScoreComparison.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')

########################################################################################################################
####################### measured responses vs intensity weighted average or contrast weighted average ##################
##################################### and then observable variance from RF weighted model ##############################
########################################################################################################################

# Initialize storage for all pair types
from sklearn.metrics import r2_score
from scipy.stats import wilcoxon
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

allMeasured = []
allRFWeightedPred = []
allContrastWeightedPred = []

for i in range(len(prefNormalized)):
    # Fit model and get predictions (as you've already done)
    resp = np.concatenate((prefNormalized[i], nonprefNormalized[i], pnNormalized[i],
                           npNormalized[i], ppNormalized[i], nnNormalized[i]), axis=0)
    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10]
    popt, pcov = curve_fit(model_wrapper, (contrast_center, contrast_periphery, locations,
                                           stim_type_center, stim_type_periphery), resp, p0=initial_guess)

    # Predicted single-stim responses
    pCenter = fullRFWeightedNorm(1, 0, -1, 1, -1, *popt)
    pPeri = [fullRFWeightedNorm(0, 1, j, -1, 1, *popt) for j in range(4)]
    npCenter = fullRFWeightedNorm(1, 0, -1, 0, -1, *popt)
    npPeri = [fullRFWeightedNorm(0, 1, j, -1, 0, *popt) for j in range(4)]

    # Paired responses
    pnFitPred = [fullRFWeightedNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
    npFitPred = [fullRFWeightedNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
    ppFitPred = [fullRFWeightedNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
    nnFitPred = [fullRFWeightedNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

    # Compare each condition
    for j in range(2):
        # PN condition: pref center, non-pref periphery
        allMeasured.append(pnNormalized[i][j])
        allRFWeightedPred.append(pnFitPred[j])
        allContrastWeightedPred.append(0.5 * (prefNormalized[i][0] + nonprefNormalized[i][j]))

        # NP condition: non-pref center, pref periphery
        allMeasured.append(npNormalized[i][j])
        allRFWeightedPred.append(npFitPred[j])
        allContrastWeightedPred.append(0.5 * (nonprefNormalized[i][0] + prefNormalized[i][j]))

        # # PP condition: pref center, pref periphery
        # allMeasured.append(ppNormalized[i][j])
        # allRFWeightedPred.append(ppFitPred[j])
        # allContrastWeightedPred.append(0.5 * (prefNormalized[i][0] + prefNormalized[i][j]))
        #
        # # NN condition: non-pref center, non-pref periphery
        # allMeasured.append(nnNormalized[i][j])
        # allRFWeightedPred.append(nnFitPred[j])
        # allContrastWeightedPred.append(0.5 * (nonprefNormalized[i][0] + nonprefNormalized[i][j]))

# Convert to arrays
measured = np.array(allMeasured)
rf_pred = np.array(allRFWeightedPred)
contrast_pred = np.array(allContrastWeightedPred)

# R² scores
r2_rf = r2_score(measured, rf_pred)
r2_contrast = r2_score(measured, contrast_pred)
print(f"RF-weighted R²: {r2_rf:.3f}")
print(f"Contrast-weighted R²: {r2_contrast:.3f}")
print(f"ΔR² (RF - Contrast): {r2_rf - r2_contrast:.3f}")

# Error comparison
rf_error = np.abs(measured - rf_pred)
contrast_error = np.abs(measured - contrast_pred)

stat, pval = wilcoxon(rf_error, contrast_error)
print(f"Wilcoxon signed-rank test p = {pval:.4g}")

# Assume measured, rf_pred, contrast_pred are already defined
total_var = np.var(measured)

resid_var_rf = np.var(measured - rf_pred)
resid_var_contrast = np.var(measured - contrast_pred)

explained_var_rf = 1 - resid_var_rf / total_var
explained_var_contrast = 1 - resid_var_contrast / total_var

gain_in_explained_var = explained_var_rf - explained_var_contrast

print(f"Explained Variance (RF-weighted): {explained_var_rf:.3f}")
print(f"Explained Variance (Contrast-weighted): {explained_var_contrast:.3f}")
print(f"Gain in Explained Variance: {gain_in_explained_var:.3f}")

# Set up the figure
plt.figure(figsize=(8, 6))
sns.set(style="whitegrid")

# Jitter for visualization
jitter = 0.1 * np.random.randn(len(rf_error))

# Plot individual residuals
plt.scatter(np.zeros_like(rf_error) + jitter, rf_error, color='royalblue', alpha=0.6, label='RF-weighted model')
plt.scatter(np.ones_like(contrast_error) + jitter, contrast_error, color='darkorange', alpha=0.6, label='Contrast-avg model')

# Plot means with error bars
means = [np.mean(rf_error), np.mean(contrast_error)]
sems = [np.std(rf_error)/np.sqrt(len(rf_error)), np.std(contrast_error)/np.sqrt(len(contrast_error))]
plt.errorbar([0, 1], means, yerr=sems, fmt='o', color='black', capsize=5)

# Beautify plot
plt.xticks([0, 1], ['RF-weighted model', 'Contrast-avg model'])
plt.ylabel('Absolute Residual (|Measured - Predicted|)')
plt.title('Residuals Comparison: RF-weighted vs. Contrast-weighted')
plt.legend()
plt.tight_layout()
plt.show()

########################################################################################################################
### variance of RF weight vs variance of simple average for normalization
########################################################################################################################

# PN, NP
pnNMI = ((prefNormalized[:, 0] + nonprefNormalized[:, 1]) - pnNormalized[:, 0]) / (
        (prefNormalized[:, 0] + nonprefNormalized[:, 1]) + pnNormalized[:, 0])

npNMI = ((nonprefNormalized[:, 0] + prefNormalized[:, 1]) - npNormalized[:, 0]) / (
        (nonprefNormalized[:, 0] + prefNormalized[:, 1]) + npNormalized[:, 0])

pnNMI2 = ((prefNormalized[:, 0] + nonprefNormalized[:, 2]) - pnNormalized[:, 1]) / (
        (prefNormalized[:, 0] + nonprefNormalized[:, 2]) + pnNormalized[:, 1])

npNMI2 = ((nonprefNormalized[:, 0] + prefNormalized[:, 2]) - npNormalized[:, 1]) / (
        (nonprefNormalized[:, 0] + prefNormalized[:, 2]) + npNormalized[:, 1])

pnNMI3 = ((prefNormalized[:, 0] + nonprefNormalized[:, 3]) - pnNormalized[:, 2]) / (
        (prefNormalized[:, 0] + nonprefNormalized[:, 3]) + pnNormalized[:, 2])

npNMI3 = ((nonprefNormalized[:, 0] + prefNormalized[:, 3]) - npNormalized[:, 2]) / (
        (nonprefNormalized[:, 0] + prefNormalized[:, 3]) + npNormalized[:, 2])

pnNMI4 = ((prefNormalized[:, 0] + nonprefNormalized[:, 4]) - pnNormalized[:, 3]) / (
        (prefNormalized[:, 0] + nonprefNormalized[:, 4]) + pnNormalized[:, 3])

npNMI4 = ((nonprefNormalized[:, 0] + prefNormalized[:, 4]) - npNormalized[:, 3]) / (
        (nonprefNormalized[:, 0] + prefNormalized[:, 4]) + npNormalized[:, 3])


plt.hist(np.concatenate((pnNMI4, npNMI4), axis=0), bins=20)
plt.xlim(-0.2, 0.8)
plt.show()

# population

# bar plots for figures
for filler in range(1):
    conditions = [1, 2, 3, 4]
    barResponsesPN = {'P Center': [prefMean[0]] * 4,
                      'PN predicted': pnFitPred,
                      'PN measured': pnMean,
                      'Average prediction': (([prefMean[0]] * 4) + nonprefMean[1:]) / 2,
                      'N offset': nonprefMean[1:]}

    barResponsesPP = {'P Center': [prefMean[0]] * 4,
                      'PP predicted': ppFitPred,
                      'PP measured': ppMean,
                      'Average prediction': (([prefMean[0]] * 4) + prefMean[1:]) / 2,
                      'N Offset': prefMean[1:]}

    barResponsesNP = {'N Center': [nonprefMean[0]] * 4,
                      'NP predicted': npFitPred,
                      'NP measured': npMean,
                      'Average prediction': (([nonprefMean[0]] * 4) + prefMean[1:]) / 2,
                      'P Offset': prefMean[1:]}

    barResponsesNN = {'N Center': [nonprefMean[0]] * 4,
                      'NN predicted': nnFitPred,
                      'NN (offset) measured': nnMean,
                      'Average prediction': (([nonprefMean[0]] * 4) + nonprefMean[1:]) / 2,
                      'N Offset': nonprefMean[1:]}

    fig, axs = plt.subplots(2, 2, layout='constrained')
    fig.set_size_inches(10, 8)
    ax1, ax2, ax3, ax4 = axs.flatten()

    colorList = ['black', 'tab:green', 'tab:red', 'tab:blue', 'gray']

    x = np.arange(len(conditions)) * 1.5  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0
    for attribute, measurement in barResponsesPN.items():
        offset = width * multiplier
        if multiplier == 0:
            rects = ax1.bar(x + offset, measurement, width, label=attribute,
                            facecolor='none', edgecolor='black')
        elif multiplier == 4:
            rects = ax1.bar(x + offset, measurement, width, label=attribute,
                            facecolor='none', edgecolor='black', hatch='//')
        else:
            rects = ax1.bar(x + offset, measurement, width, label=attribute,
                            color=colorList[multiplier])
        multiplier += 1

    ax1.set_ylabel('Normalized Responses', fontsize=15)
    ax1.set_xticks(x + width, conditions)
    # ax1.legend(loc='upper left', ncols=5)
    ax1.set_ylim(0, 1.25)

    multiplier = 0
    for attribute, measurement in barResponsesPP.items():
        offset = width * multiplier
        if multiplier == 0:
            rects = ax2.bar(x + offset, measurement, width, label=attribute,
                            facecolor='none', edgecolor='black')
        elif multiplier == 4:
            rects = ax2.bar(x + offset, measurement, width, label=attribute,
                            facecolor='none', edgecolor='black', hatch='//')
        else:
            rects = ax2.bar(x + offset, measurement, width, label=attribute,
                            color=colorList[multiplier])
        multiplier += 1

    ax2.set_ylabel('Normalized Responses', fontsize=15)
    ax2.set_xticks(x + width, conditions)
    # ax2.legend(loc='upper left', ncols=5)
    ax2.set_ylim(0, 1.25)

    multiplier = 0
    for attribute, measurement in barResponsesNP.items():
        offset = width * multiplier
        if multiplier == 0:
            rects = ax3.bar(x + offset, measurement, width, label=attribute,
                            facecolor='none', edgecolor='black')
        elif multiplier == 4:
            rects = ax3.bar(x + offset, measurement, width, label=attribute,
                            facecolor='none', edgecolor='black', hatch='//')
        else:
            rects = ax3.bar(x + offset, measurement, width, label=attribute,
                            color=colorList[multiplier])
        multiplier += 1

    ax3.set_ylabel('Normalized Responses', fontsize=15)
    ax3.set_xticks(x + width, conditions)
    # ax3.legend(loc='upper left', ncols=5)
    ax3.set_ylim(0, 1.25)

    multiplier = 0
    for attribute, measurement in barResponsesNN.items():
        offset = width * multiplier
        if multiplier == 0:
            rects = ax4.bar(x + offset, measurement, width, label=attribute,
                            facecolor='none', edgecolor='black')
        elif multiplier == 4:
            rects = ax4.bar(x + offset, measurement, width, label=attribute,
                            facecolor='none', edgecolor='black', hatch='//')
        else:
            rects = ax4.bar(x + offset, measurement, width, label=attribute,
                            color=colorList[multiplier])
        multiplier += 1

    ax4.set_ylabel('Normalized Responses', fontsize=15)
    ax4.set_xticks(x + width, conditions)
    ax4.legend(loc='upper left', ncols=2)
    ax4.set_ylim(0, 1.25)

    plt.tight_layout()
    plt.show()


# figure/plotting
for filler in range(1):
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
    fig.set_size_inches(20, 6)
    popt_str = (f"Lp = {popt[0]:.2f}\n"
                f"Lnp = {popt[1]:.2f}\n"
                f"W1 = {popt[2]:.2f}\n"
                f"W2 = {popt[3]:.2f}\n"
                f"W3 = {popt[4]:.2f}\n"
                f"W4 = {popt[5]:.2f}\n"
                f"sigma = {popt[6]:.2f}")
    ax1.plot(x, prefMean, color='black', label='PO')
    ax1.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax1.scatter(x, nonprefMean, color='grey', label='NO', hatch='//')
    ax1.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    ax1.plot(xNorm, pnMean, color='green', label='P0 N1')
    ax1.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
                 color='green', markersize=7)
    ax1.plot(xNorm, npMean, color='red', label='N0 P1')
    ax1.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='red',
                 color='red', markersize=7)
    ax1.plot(xNorm, ppMean, color='black', label='P0 P1')
    ax1.errorbar(xNorm, ppMean, yerr=ppSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax1.plot(xNorm, nnMean, color='grey', label='N0 N1')
    ax1.errorbar(xNorm, nnMean, yerr=nnSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    ax1.plot(xNorm, pnFitPred, color='green', linestyle='--')
    ax1.scatter(xNorm, pnFitPred, facecolors='none', edgecolors='g', s=50)
    ax1.plot(xNorm, npFitPred, color='red', linestyle='--')
    ax1.scatter(xNorm, npFitPred, facecolors='none', edgecolors='r', s=50)
    ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', alpha=0.8)
    ax1.scatter(xNorm, ppFitPred, facecolors='none', edgecolors='black', s=50)
    ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', alpha=0.8)
    ax1.scatter(xNorm, nnFitPred, facecolors='none', edgecolors='gray', s=50)
    ax1.plot(x, pSingleFitPred, color='black', linestyle='--')
    ax1.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
    ax1.plot(x, npSingleFitPred, color='grey', linestyle='--')
    ax1.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
    ax1.set_xlabel('Receptive Field Location', fontsize=25)
    ax1.set_ylabel('Normalized Response Rate', fontsize=25)
    ax1.set_ylim([0, 1.7])
    # ax1.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
    # ax1.fill_between(x=x, y1=meanSpon + semSpon, y2=meanSpon - semSpon,
    #                  color='blue', alpha=0.2)
    # Add a textbox with the popt values to the plot
    ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
             transform=ax1.transAxes)
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.set_xticks([0, 1, 2, 3, 4])
    ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4.5])

    ax2.scatter(measuredResp, fullRFWeightPred, color='black', facecolors= 'none')
    ax2.set_ylim([0, 2])
    ax2.set_xlim([0, 2])
    ax2.set_xlabel('Measured Response', fontsize=15, **hfont)
    ax2.set_ylabel('Predicted Response', fontsize=15, **hfont)
    line = lines.Line2D([0, 1], [0, 1], color='black')
    transform = ax2.transAxes
    line.set_transform(transform)
    ax2.add_line(line)

    ax3.hist(r2ScoresRFWeightFull, bins=5)
    ax3.set_xlim([0, 1])
    ax3.set_xlabel('Explained R2 RF Weighted Model', fontsize=15, **hfont)
    ax3.set_ylabel('Frequency', fontsize=15, **hfont)

    ax4.hist(r2ScoresSimpleAvgFull, bins=50)
    ax4.set_xlim([0, 1])
    ax4.set_ylim([0, 35])
    ax4.set_xlabel('Explained R2 Simple Avg Model', fontsize=15, **hfont)
    ax4.set_ylabel('Frequency', fontsize=15, **hfont)

    plt.tight_layout()
    plt.show()


########################################################################################################################
############# Evaluate Weighted Normalization vs Weighted Normalization Heuristic Model with baseline ##################
########################################################################################################################

# known variables
contrast_center = [1, 0, 0, 0, 0,
                   1, 0, 0, 0, 0,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   0]  # First response has contrast in the center, second in periphery
contrast_periphery = [0, 1, 1, 1, 1,
                      0, 1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      0]  # First has no periphery, second has periphery with contrast
locations = np.array([-1, 0, 1, 2, 3,
                      -1, 0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      -1])  # First has no peripheral stimulus, second is at location 1
stim_type_center = np.array([1, -1, -1, -1, -1,
                             0, -1, -1, -1, -1,
                             1, 1, 1, 1,
                             0, 0, 0, 0,
                             1, 1, 1, 1,
                             0, 0, 0, 0,
                             1])
stim_type_periphery = np.array([-1, 1, 1, 1, 1,
                                -1, 0, 0, 0, 0,
                                0, 0, 0, 0,
                                1, 1, 1, 1,
                                1, 1, 1, 1,
                                0, 0, 0, 0,
                                1])


# scatter plot of individual neuron fit r2 scores for heuristic weighted vs rf weighted
for filler in range(1):
    fix, ax1 = plt.subplots()
    ax1.scatter(r2ScoresRFWeightBase, r2ScoresRFWeightHeuristic)
    ax1.set_xlabel('r2 Scores weighted normalization')
    ax1.set_ylabel('r2 Scores heuristic weighted norm')

    x_min, x_max = ax1.get_xlim()
    y_min, y_max = ax1.get_ylim()
    min_val = min(x_min, y_min)
    max_val = max(x_max, y_max)
    ax1.plot([min_val, max_val], [min_val, max_val], 'k--', label='Unity Line', zorder=0)

    # Reset limits to ensure the unity line fits perfectly
    ax1.set_xlim([min_val, max_val])
    ax1.set_ylim([min_val, max_val])
    stat, p_value = wilcoxon(r2ScoresRFWeightBase, r2ScoresRFWeightHeuristic, alternative='two-sided')
    print(p_value)
    plt.show()


# population fits
resp = np.concatenate((prefMeanWithBase, nonprefMeanWithBase, pnMeanWithBase,
                       npMeanWithBase, ppMeanWithBase, nnMeanWithBase,
                       [np.mean(sponNormalizedToPrefCenter)]), axis=0)
initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                 np.mean(sponNormalizedToPrefCenter)]  # Lp, Lnp, W1, W2, W3, W4, sigma, b

# weighted norm heuristic model
popt, pcov = curve_fit(model_wrapperHeuristic, (contrast_center, contrast_periphery, locations, stim_type_center,
                                                stim_type_periphery), resp,
                       p0=initial_guess)
y_pred1 = apply_fitted_modelHeuristic(contrast_center, contrast_periphery, locations,
                                      stim_type_center, stim_type_periphery, *popt)
r2 = r2_score(resp, y_pred1)
print(r2)
print(popt)

# weighted norm with baseline
popt, pcov = curve_fit(model_wrapperWithBase, (contrast_center, contrast_periphery, locations, stim_type_center,
                                               stim_type_periphery), resp,
                       p0=initial_guess)
y_pred1 = apply_fitted_modelWithBase(contrast_center, contrast_periphery, locations,
                                     stim_type_center, stim_type_periphery, *popt)
r2 = r2_score(resp, y_pred1)
print(r2)
print(popt)

pnFitPred = [fullRFWeightedNormWithBase(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [fullRFWeightedNormWithBase(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullRFWeightedNormWithBase(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullRFWeightedNormWithBase(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullRFWeightedNormWithBase(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullRFWeightedNormWithBase(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullRFWeightedNormWithBase(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullRFWeightedNormWithBase(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# RF weighted with baseline prediction figure
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    popt_str = (f"Lp = {popt[0]:.2f}\n"
                f"Lnp = {popt[1]:.2f}\n"
                f"W1 = {popt[2]:.2f}\n"
                f"W2 = {popt[3]:.2f}\n"
                f"W3 = {popt[4]:.2f}\n"
                f"W4 = {popt[5]:.2f}\n"
                f"sigma = {popt[6]:.2f}")
    ax1.plot(x, prefMeanWithBase, color='black', label='PO', linewidth=1.5)
    ax1.errorbar(x, prefMeanWithBase, yerr=prefWithBaseSEM, fmt='o', ecolor='black',
                 markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax1.plot(x, nonprefMeanWithBase, color='gray', label='NO', linewidth=1.5)
    ax1.errorbar(x, nonprefMeanWithBase, yerr=nonprefWithBaseSEM, fmt='o', ecolor='grey',
                 markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax1.plot(xNorm, pnMeanWithBase, color='green', label='P0 N1', linewidth=1.5)
    ax1.errorbar(xNorm, pnMeanWithBase, yerr=pnWithBaseSEM, fmt='o', ecolor='green',
                 color='green', markersize=7)
    ax1.plot(xNorm, npMeanWithBase, color='magenta', label='N0 P1', linewidth=1.5)
    ax1.errorbar(xNorm, npMeanWithBase, yerr=npWithBaseSEM, fmt='o', ecolor='magenta',
                 color='magenta', markersize=7)
    ax1.plot(xNorm, ppMeanWithBase, color='black', label='P0 P1', linewidth=1.5)
    ax1.errorbar(xNorm, ppMeanWithBase, yerr=ppWithBaseSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax1.plot(xNorm, nnMeanWithBase, color='grey', label='N0 N1', linewidth=1.5)
    ax1.errorbar(xNorm, nnMeanWithBase, yerr=nnWithBaseSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, pnFitPred, facecolors='g', edgecolors='g', s=50)
    ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, npFitPred, facecolors='m', edgecolors='m', s=50)
    ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=50)
    ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, nnFitPred, facecolors='grey', edgecolors='gray', s=50)
    ax1.plot(x, pSingleFitPred, color='black', linestyle='--')
    ax1.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
    ax1.plot(x, npSingleFitPred, color='grey', linestyle='--')
    ax1.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
    ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    ax1.set_ylabel('Normalized Response Rate', fontsize=25, **hfont)
    ax1.set_ylim([-0.2, 1.5])
    # Add a textbox with the popt values to the plot
    ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
             transform=ax1.transAxes)
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.set_xticks([0, 1, 2, 3, 4])
    ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4.5])
    ax1.tick_params(axis='both', width=2, length=8)

    plt.show()


# weighted norm with baseline (heuristic model)
resp = np.concatenate((prefMeanWithBase, nonprefMeanWithBase, pnMeanWithBase,
                       npMeanWithBase, ppMeanWithBase, nnMeanWithBase,
                       [np.mean(sponNormalizedToPrefCenter)]), axis=0)
initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10,
                 np.mean(sponNormalizedToPrefCenter)]  # Lp, Lnp, W1, W2, W3, W4, sigma, b

popt, pcov = curve_fit(model_wrapperHeuristic, (contrast_center, contrast_periphery, locations, stim_type_center,
                                                stim_type_periphery), resp,
                       p0=initial_guess)
y_pred1 = apply_fitted_modelHeuristic(contrast_center, contrast_periphery, locations,
                                     stim_type_center, stim_type_periphery, *popt)
r2 = r2_score(resp, y_pred1)
print(r2)
print(popt)

pnFitPred = [fullRFWeightedNormHeuristic(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [fullRFWeightedNormHeuristic(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullRFWeightedNormHeuristic(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullRFWeightedNormHeuristic(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullRFWeightedNormHeuristic(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullRFWeightedNormHeuristic(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullRFWeightedNormHeuristic(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullRFWeightedNormHeuristic(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# RF weighted with baseline prediction figure (heuristic model)
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    popt_str = (f"Lp = {popt[0]:.2f}\n"
                f"Lnp = {popt[1]:.2f}\n"
                f"W1 = {popt[2]:.2f}\n"
                f"W2 = {popt[3]:.2f}\n"
                f"W3 = {popt[4]:.2f}\n"
                f"W4 = {popt[5]:.2f}\n"
                f"sigma = {popt[6]:.2f}")
    ax1.plot(x, prefMeanWithBase, color='black', label='PO', linewidth=1.5)
    ax1.errorbar(x, prefMeanWithBase, yerr=prefWithBaseSEM, fmt='o', ecolor='black',
                 markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax1.plot(x, nonprefMeanWithBase, color='gray', label='NO', linewidth=1.5)
    ax1.errorbar(x, nonprefMeanWithBase, yerr=nonprefWithBaseSEM, fmt='o', ecolor='grey',
                 markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax1.plot(xNorm, pnMeanWithBase, color='green', label='P0 N1', linewidth=1.5)
    ax1.errorbar(xNorm, pnMeanWithBase, yerr=pnWithBaseSEM, fmt='o', ecolor='green',
                 color='green', markersize=7)
    ax1.plot(xNorm, npMeanWithBase, color='magenta', label='N0 P1', linewidth=1.5)
    ax1.errorbar(xNorm, npMeanWithBase, yerr=npWithBaseSEM, fmt='o', ecolor='magenta',
                 color='magenta', markersize=7)
    ax1.plot(xNorm, ppMeanWithBase, color='black', label='P0 P1', linewidth=1.5)
    ax1.errorbar(xNorm, ppMeanWithBase, yerr=ppWithBaseSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax1.plot(xNorm, nnMeanWithBase, color='grey', label='N0 N1', linewidth=1.5)
    ax1.errorbar(xNorm, nnMeanWithBase, yerr=nnWithBaseSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, pnFitPred, facecolors='g', edgecolors='g', s=50)
    ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, npFitPred, facecolors='m', edgecolors='m', s=50)
    ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=50)
    ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax1.scatter(xNorm, nnFitPred, facecolors='grey', edgecolors='gray', s=50)
    ax1.plot(x, pSingleFitPred, color='black', linestyle='--')
    ax1.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
    ax1.plot(x, npSingleFitPred, color='grey', linestyle='--')
    ax1.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
    ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    ax1.set_ylabel('Normalized Response Rate', fontsize=25, **hfont)
    ax1.set_ylim([-0.2, 1.5])
    # Add a textbox with the popt values to the plot
    ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
             transform=ax1.transAxes)
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.set_xticks([0, 1, 2, 3, 4])
    ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4.5])
    ax1.tick_params(axis='both', width=2, length=8)

    plt.show()


########################################################################################################################
###################################################### Work Space  #####################################################
################################################# Explainable Variance #################################################
########################################################################################################################
explainable_variances = []
goodI = []

for i in range(len(prefNormalized)):
    # Skip if any entry is zero-dimensional or empty
    if (prefNormalized[i].ndim == 0 or nonprefNormalized[i].ndim == 0 or
        pnNormalized[i].ndim == 0 or npNormalized[i].ndim == 0 or
        ppNormalized[i].ndim == 0 or nnNormalized[i].ndim == 0):
        continue

    # Convert normalized rates to raw rates
    raw_pref = prefNormalized[i] * normVal[i]
    raw_nonpref = nonprefNormalized[i] * normVal[i]
    raw_pn = pnNormalized[i] * normVal[i]
    raw_np = npNormalized[i] * normVal[i]
    raw_pp = ppNormalized[i] * normVal[i]
    raw_nn = nnNormalized[i] * normVal[i]

    # Concatenate into a single response vector
    raw_resp = np.concatenate((raw_pref, raw_nonpref, raw_pn, raw_np, raw_pp, raw_nn), axis=0)

    # Fit model to raw responses
    initial_guess = [1.0 * float(normVal[i]), 0.5 * float(normVal[i]),
                     1.0, 0.75, 0.50, 0.10, 0.10]
    try:
        popt, _ = curve_fit(model_wrapper,
                            (contrast_center, contrast_periphery, locations, stim_type_center, stim_type_periphery),
                            raw_resp, p0=initial_guess)
    except RuntimeError:
        continue  # Skip neuron if fit fails

    # Predicted rates in spikes/s
    predicted_rates = np.array(apply_fitted_model(contrast_center, contrast_periphery, locations,
                                                  stim_type_center, stim_type_periphery, *popt))

    # Simulate trials using Poisson spike counts (convert to spikes per 250 ms)
    n_repeats = allBlocksDone[i]
    rate_per_trial = predicted_rates * 0.25

    # Ensure all λ values are valid for Poisson
    rate_per_trial = np.nan_to_num(rate_per_trial, nan=0.0, posinf=0.0, neginf=0.0)
    rate_per_trial[rate_per_trial < 0] = 0.0

    simulated_trials = np.random.poisson(rate_per_trial[:, None], size=(len(rate_per_trial), n_repeats))

    # Average simulated responses
    simulated_means = simulated_trials.mean(axis=1)

    # # # Convert back to spikes/s
    simulated_means = simulated_means / 0.25

    # Refit model to simulated averages
    try:
        popt_sim, _ = curve_fit(model_wrapper,
                                (contrast_center, contrast_periphery, locations, stim_type_center, stim_type_periphery),
                                simulated_means, p0=initial_guess)
    except RuntimeError:
        continue

    refit_preds = np.array(apply_fitted_model(contrast_center, contrast_periphery, locations,
                                              stim_type_center, stim_type_periphery, *popt_sim))

    # Compute explainable variance (R²)
    ss_resid = np.sum((simulated_means - refit_preds) ** 2)
    ss_total = np.sum((simulated_means - np.mean(simulated_means)) ** 2)
    R2 = 1 - (ss_resid / ss_total)
    explainable_variances.append(R2)
    goodI.append(i)

explainable_variances = np.array(explainable_variances)
# Population-level explainable variance
population_explainable_variance = np.median(explainable_variances)
print("Population-level explainable variance:", population_explainable_variance)


# Scatter plot
plt.figure(figsize=(5, 5))
plt.scatter(explainable_variances, r2ScoresRFWeightFull[goodI], color='blue', alpha=0.7)

minVal = 0.7
maxVal = 1
# Axis limits
plt.xlim(minVal, maxVal)
plt.ylim(minVal, maxVal)

# Unity line
plt.plot([minVal, maxVal], [minVal, maxVal], 'k--')

# Labels and title
plt.xlabel("Explainable variance (R² from simulated data)")
plt.ylabel("Model fit R² (original data)")
plt.title("Explainable variance vs. Original model R²")
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# Compute ratio
ratio = r2ScoresRFWeightFull / explainable_variances

# Plot using matplotlib
plt.figure(figsize=(8, 5))
plt.hist(ratio, bins=20, color="royalblue", edgecolor="black")
plt.axvline(1.0, color='red', linestyle='--', linewidth=2, label='Unity (1.0)')
plt.xlabel("Model R² / Explainable Variance", fontsize=12)
plt.ylabel("Count", fontsize=12)
plt.title("Distribution of Model R² as a Fraction of Explainable Variance", fontsize=14)
plt.legend()
plt.tight_layout()
plt.show()

########################################################################################################################
###################################################### Extra Code ######################################################
########################################################################################################################
# population plots
prefNormalized = np.array(prefNormalized)
nonprefNormalized = np.array(nonprefNormalized)
pnNormalized = np.array(pnNormalized)
npNormalized = np.array(npNormalized)
ppNormalized = np.array(ppNormalized)
nnNormalized = np.array(nnNormalized)
transectNormalized = np.array(transectNormalized)
pnNMIPop = np.array(pnNMIPop)
npNMIPop = np.array(npNMIPop)
ppNMIPop = np.array(ppNMIPop)
nnNMIPop = np.array(nnNMIPop)
offsetDegSepNormPop = np.array(offsetDegSepNormPop)
rfGaborSigmaPop = np.array(rfGaborSigmaPop)

numSteps = 4
numUnits = len(sponNormalized)
offsetDegSepNormPop = offsetDegSepNormPop.reshape((numUnits, numSteps*2+1))
meanSpon = np.mean(sponNormalized)
semSpon = np.std(sponNormalized) / np.sqrt(numUnits)

# main figure
for filler in range(1):
    fig, ((ax2, ax3, ax1, ax4), (ax7, ax8, ax5, ax6)) = plt.subplots(2, 4)
    fig.set_size_inches(16, 8)
    hfont = {'fontname': 'Arial'}

    # bin size = n, exponent for prediction = exp
    n = numUnits
    exp = 1
    fig.suptitle(f'Population responses pref+nonpref neurons, n={numUnits}', **hfont, fontsize=15)

    ## key variables
    x = np.arange(0, numSteps + 1)
    xNorm = np.arange(1, numSteps + 1)
    prefMean = np.mean(prefNormalized, axis=0)
    prefSEM = np.std(prefNormalized, axis=0) / np.sqrt(numUnits)
    nonprefMean = np.mean(nonprefNormalized, axis=0)
    nonprefSEM = np.std(nonprefNormalized, axis=0) / np.sqrt(numUnits)
    pnMean = np.mean(pnNormalized, axis=0)
    pnSEM = np.std(pnNormalized, axis=0) / np.sqrt(numUnits)
    npMean = np.mean(npNormalized, axis=0)
    npSEM = np.std(npNormalized, axis=0) / np.sqrt(numUnits)
    transX = np.concatenate((np.arange(-numSteps, 0),
                             np.arange(0, numSteps + 1)),
                            axis=0)
    transectMean = np.mean(transectNormalized, axis=0)
    transectSEM = np.std(transectNormalized, axis=0) / np.sqrt(numUnits)
    params = gaussFit(transX, transectMean)
    xFull = np.linspace(transX[0], transX[-1], 1000)
    respFull = gauss(xFull, *params)

    # transect response
    ax1.plot(xFull, respFull, color='black', linestyle='--', label='Gauss Fit')
    ax1.plot(transX, transectMean, color='black', label='Pref')
    ax1.errorbar(transX, transectMean, yerr=transectSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax1.plot(x, nonprefMean, color='grey', label='NO')
    ax1.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    # for i in transectNormalized:
    #     ax1.plot(transX, i, color='purple', alpha=0.2)
    ax1.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
    ax1.fill_between(x=transX, y1=meanSpon+semSpon, y2=meanSpon-semSpon,
                     color='blue', alpha=0.2)
    ax1.set_ylabel('Normalized Response', fontsize=15, **hfont)

    ax1.set_xlabel('Receptive Field Offset', fontsize=15, **hfont)
    ax1.set_ylim([0, 1.5])
    ax1.set_title('transect')


    ## normalization with distance PN/NP
    # # RF weight (from fitted Gaussian)
    # halfTransX = np.arange(0, 5)
    # respHalfGauss = gauss(halfTransX, *params)
    # # drivenRespHalfGauss = respHalfGauss - params[0]
    # drivenRespHalfGauss = respHalfGauss
    # weights = drivenRespHalfGauss / drivenRespHalfGauss[0]

    # weight from raw responses
    weights = prefMean / prefMean[0]

    # rf weighted average
    # RF weight is by max response of pref stimulus in transect
    # pCentPred = ((prefMean[0] + nonprefMean[1:])**1) / (
    #             1 + prefMean[1:] / np.max(prefMean))
    # npCentPred = ((nonprefMean[0] + prefMean[1:])**1) / (
    #             1 + prefMean[1:] / np.max(prefMean))

    # RF weight is by max of fitted Gaussian
    # pCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / (
    #             ((prefMean[0] + prefMean[1:])**1) / np.max(respFull))
    # npCentPred = ((nonprefMean[0] + prefMean[1:])**exp) / (
    #              ((prefMean[0] + prefMean[1:])**1) / np.max(respFull))

    # pCentPred = (prefMean[0] / np.max(respFull)) / ((prefMean[0] + prefMean[1:]) / np.max(respFull)) + (
    #             (nonprefMean[1:] * (prefMean[1:] / np.max(respFull))) / ((prefMean[0] + prefMean[1:]) / np.max(respFull)))
    #
    # npCentPred = (nonprefMean[0] / np.max(respFull)) / ((prefMean[0] + prefMean[1:]) / np.max(respFull)) + (
    #              (prefMean[1:] * (prefMean[1:] / np.max(respFull))) / ((prefMean[0] + prefMean[1:]) / np.max(respFull)))

    # # RF weight is applied to numerator as well (USE THIS)
    # pCentPred = ((prefMean[0] / np.max(respFull)) + (nonprefMean[1:] * (prefMean[1:]/np.max(respFull)))) / (
    #             (prefMean[0] + prefMean[1:]) / np.max(respFull))
    #
    # npCentPred = ((nonprefMean[0] / np.max(respFull)) + (prefMean[1:] * (prefMean[1:]/np.max(respFull)))) / (
    #              (prefMean[0] + prefMean[1:]) / np.max(respFull))

    # # # RF weight is applied to numerator as well (USE THIS) : Raw Rate
    # pCentPred = ((prefMean[0] * (prefMean[0] / np.max(respFull))) + (nonprefMean[1:] * (prefMean[1:]/np.max(respFull)))) / (
    #             (prefMean[0] + prefMean[1:]) / np.max(transectMean))
    #
    # npCentPred = ((nonprefMean[0] * (prefMean[0] / np.max(respFull))) + (prefMean[1:] * (prefMean[1:]/np.max(respFull)))) / (
    #              (prefMean[0] + prefMean[1:]) / np.max(respFull))

    # # RF weight is applied to numerator as well : Driven Rate
    # pCentPred = (((prefMean[0] - meanSpon) * ((prefMean[0] - meanSpon) / (np.max(respFull) - meanSpon))) + ((nonprefMean[1:] - meanSpon) * ((prefMean[1:] - meanSpon)/(np.max(respFull)-meanSpon)))) / (
    #             ((prefMean[0] - meanSpon) + (prefMean[1:] - meanSpon)) / (np.max(respFull) - meanSpon)) + meanSpon
    #
    # npCentPred = (((nonprefMean[0] - meanSpon) * ((prefMean[0] - meanSpon) / (np.max(respFull) - meanSpon))) + ((prefMean[1:] - meanSpon) * ((prefMean[1:] - meanSpon) / (np.max(respFull)-meanSpon)))) / (
    #              ((prefMean[0] - meanSpon) + (prefMean[1:] - meanSpon)) / (np.max(respFull)) - meanSpon) + meanSpon

    # # RF weight is applied to numerator as well (OLD WAY)
    # pCentPred = (prefMean[0] + (nonprefMean[1:]) / np.max(respFull)) / (
    #             (prefMean[0] + prefMean[1:]) / np.max(respFull))
    # npCentPred = (nonprefMean[0] + (prefMean[1:]) / np.max(respFull)) / (
    #              (prefMean[0] + prefMean[1:]) / np.max(respFull))

    ## Best fit but ^2 the responses so not great for capturing cross orientation suppression (CLOSEST) ####
    # pCentPred = (((prefMean[0] * (prefMean[0] / np.max(transectMean))) + (nonprefMean[1:] * (nonprefMean[1:]/np.max(transectMean)))) / (
    #             prefMean[0] / np.max(transectMean) + nonprefMean[1:] / np.max(transectMean))) + (meanSpon/2)
    #
    # npCentPred = (((nonprefMean[0] * (nonprefMean[0] / np.max(transectMean))) + (prefMean[1:] * (prefMean[1:]/np.max(transectMean)))) / (
    #               nonprefMean[0] / np.max(transectMean) + prefMean[1:] / np.max(transectMean))) + (meanSpon/2)

    # # RF weight is applied to numerator as well (USE THIS) : Raw Rate
    # pCentPred = ((prefMean[0] * weights[0]) + (nonprefMean[1:] * weights[1:]) + meanSpon) / (
    #             weights[0] + weights[1:])
    #
    # npCentPred = ((nonprefMean[0] * weights[0]) + (prefMean[1:] * weights[1:]) + meanSpon) / (
    #              weights[0] + weights[1:])

    # # # RF weight is applied to numerator as well (USE THIS) : Raw Rate (maybe good)
    # pCentPred = ((prefMean[0] * weights[0]) + (nonprefMean[1:] * weights[1:]) + (meanSpon * (weights[0] + weights[1:] + 0.07))) / (
    #             weights[0] + weights[1:] + 0.07)
    #
    # npCentPred = ((nonprefMean[0] * weights[0]) + (prefMean[1:] * weights[1:]) + (meanSpon * (weights[0] + weights[1:] + 0.07))) / (
    #              weights[0] + weights[1:] + 0.07)

    pCentPred = (((prefMean[0]) + (nonprefMean[1:])) / (
                weights[0] + weights[1:]))

    npCentPred = (((nonprefMean[0]) + (prefMean[1:])) / (
                 weights[0] + weights[1:]))

    # pCentPred = (((prefMean[0]) + (nonprefMean[0] * (weights[1:]**1))) / (
    #             weights[0] + weights[1:]) + 0.07)
    #
    # npCentPred = (((nonprefMean[0]) + (prefMean[0] * (weights[1:]**1))) / (
    #              weights[0] + weights[1:]) + 0.07)

    # simple average
    # pCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / 2
    # npCentPred = ((nonprefMean[0] + prefMean[1:])**exp) / 2

    ax2.plot(x, prefMean, color='black', label='PO')
    ax2.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax2.plot(x, nonprefMean, color='grey', label='NO')
    ax2.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    ax2.plot(xNorm, pnMean, color='green', label='P0 N1')
    ax2.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
                 color='green', markersize=7)
    ax2.plot(xNorm, npMean, color='red', label='N0 P1')
    ax2.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='red',
                 color='red', markersize=7)
    ax2.plot(xNorm, pCentPred, color='green', linestyle='--')
    ax2.scatter(xNorm, pCentPred, facecolors='none', edgecolors='g', s=50)
    ax2.plot(xNorm, npCentPred, color='red', linestyle='--')
    ax2.scatter(xNorm, npCentPred, facecolors='none', edgecolors='r', s=50)
    ax2.set_xlabel('Receptive Field Offset', fontsize=15)
    ax2.set_ylabel('Normalized Response Rate', fontsize=15)
    ax2.set_ylim([0, 1.5])
    ax2.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
    ax2.fill_between(x=x, y1=meanSpon+semSpon, y2=meanSpon-semSpon,
                     color='blue', alpha=0.2)
    ax2.set_title('Normalization vs Distance PN/NP')

    # normalization with distance PP/NN
    x = np.arange(0, numSteps + 1)
    xNorm = np.arange(1, numSteps + 1)
    prefMean = np.mean(prefNormalized, axis=0)
    prefSEM = np.std(prefNormalized, axis=0) / np.sqrt(numUnits)
    nonprefMean = np.mean(nonprefNormalized, axis=0)
    nonprefSEM = np.std(nonprefNormalized, axis=0) / np.sqrt(numUnits)
    ppMean = np.mean(ppNormalized, axis=0)
    pnSEM = np.std(ppNormalized, axis=0) / np.sqrt(numUnits)
    nnMean = np.mean(nnNormalized, axis=0)
    nnSEM = np.std(nnNormalized, axis=0) / np.sqrt(numUnits)

    # RF weighted Avg (best fit but not great for capturing cross orientation suppression because
    # squaring the responses
    # ppCentPred = (((prefMean[0] * (prefMean[0] / np.max(transectMean))) + (prefMean[1:] * (prefMean[1:] / np.max(transectMean)))) / (
    #               (prefMean[0] + prefMean[1:]) / np.max(transectMean))) + (meanSpon/2)
    # npnpCentPred = (((nonprefMean[0] * (nonprefMean[0] / np.max(transectMean))) + (nonprefMean[1:] * (nonprefMean[1:]/np.max(transectMean)))) / (
    #                 (nonprefMean[0] + nonprefMean[1:]) / np.max(transectMean))) + (meanSpon/2)

    # RF weighted Avg
    # ppCentPred = ((((prefMean[0] - meanSpon) * weights[0]) + ((prefMean[1:] - meanSpon) * weights[1:])) / (
    #             weights[0] + weights[1:])) + meanSpon
    #
    # npnpCentPred = ((((nonprefMean[0] - meanSpon) * weights[0]) + ((nonprefMean[1:] - meanSpon) * weights[1:])) / (
    #              weights[0] + weights[1:])) + meanSpon

    # ppCentPred = ((prefMean[0] * weights[0]) + (prefMean[1:] * weights[1:]) + meanSpon) / (
    #             weights[0] + weights[1:])
    #
    # npnpCentPred = ((nonprefMean[0] * weights[0]) + (nonprefMean[1:] * weights[1:]) + meanSpon) / (
    #              weights[0] + weights[1:])

    ppCentPred = ((prefMean[0] * weights[0]) + (prefMean[1:] * weights[1:]) + (meanSpon * (weights[0] + weights[1:] + 0.07))) / (
                weights[0] + weights[1:] + 0.07)

    npnpCentPred = ((nonprefMean[0] * weights[0]) + (nonprefMean[1:] * weights[1:]) + (meanSpon * (weights[0] + weights[1:] + 0.07))) / (
                 weights[0] + weights[1:] + 0.07)


    ax3.plot(x, prefMean, color='black', label='PO')
    ax3.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
                 color='black', markersize=2)
    ax3.plot(x, nonprefMean, color='grey', label='NO')
    ax3.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=2)
    ax3.plot(xNorm, ppMean, color='black', label='P0 P1', linestyle='--')
    ax3.errorbar(xNorm, ppMean, yerr=pnSEM, fmt='o', ecolor='black',
                 color='black', markersize=2)
    ax3.plot(xNorm, nnMean, color='grey', label='N0 N1', linestyle='--')
    ax3.errorbar(xNorm, nnMean, yerr=npSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=2)
    ax3.plot(xNorm, ppCentPred, color='black', linestyle=':', alpha=0.8)
    ax3.plot(xNorm, npnpCentPred, color='grey', linestyle=':', alpha=0.8)
    ax3.set_xlabel('stimulus offset positions')
    ax3.set_ylabel('normalized response')
    ax3.set_ylim([0, 1.5])
    ax3.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
    ax3.fill_between(x=x, y1=meanSpon+semSpon, y2=meanSpon-semSpon,
                     color='blue', alpha=0.2)
    ax3.set_title('Normalization vs Distance PP/NN')

    # # NMI pop plot
    # pnNMIPopMean = np.mean(pnNMIPop, axis=0)
    # pnNMISEM = np.std(pnNMIPop, axis=0) / np.sqrt(numUnits)
    # npNMIPopMean = np.mean(npNMIPop, axis=0)
    # npNMISEM = np.std(npNMIPop, axis=0) / np.sqrt(numUnits)
    # ppNMIPopMean = np.mean(ppNMIPop, axis=0)
    # ppNMISEM = np.std(ppNMIPop, axis=0) / np.sqrt(numUnits)
    # nnNMIPopMean = np.mean(nnNMIPop, axis=0)
    # nnNMISEM = np.std(nnNMIPop, axis=0) / np.sqrt(numUnits)
    #
    # ax4.plot(xNorm, pnNMIPopMean, color='green', label='PN')
    # ax4.errorbar(xNorm, pnNMIPopMean, yerr=pnNMISEM, fmt='o', ecolor='green',
    #              color='green', markersize=2)
    # ax4.scatter(xNorm, pnNMIPopMean, color='green')
    # ax4.plot(xNorm, npNMIPopMean, color='red', label='NP')
    # ax4.errorbar(xNorm, npNMIPopMean, yerr=npNMISEM, fmt='o', ecolor='red',
    #              color='red', markersize=2)
    # ax4.scatter(xNorm, npNMIPopMean, color='red')
    # ax4.plot(xNorm, ppNMIPopMean, color='black', label='PP', linestyle='--')
    # ax4.errorbar(xNorm, ppNMIPopMean, yerr=ppNMISEM, fmt='o', ecolor='black',
    #              color='black', markersize=2)
    # ax4.scatter(xNorm, ppNMIPopMean, color='black')
    # ax4.plot(xNorm, nnNMIPopMean, color='grey', label='NN', linestyle='--')
    # ax4.errorbar(xNorm, nnNMIPopMean, yerr=nnNMISEM, fmt='o', ecolor='grey',
    #              color='grey', markersize=2)
    # ax4.scatter(xNorm, nnNMIPopMean, color='grey')
    # ax4.set_xlabel('stimulus offset position')
    # ax4.set_ylabel('NMI')
    # ax4.set_ylim(bottom=0)
    # ax4.set_title('NMI, Ni and Maunsell, 2017')

    # NMI vs normalized separation
    normSep = offsetDegSepNormPop[:, 5:].reshape(numUnits*numSteps)
    pn = pnNMIPop.reshape(numUnits*numSteps)
    nonp = npNMIPop.reshape(numUnits*numSteps)
    pp = ppNMIPop.reshape(numUnits*numSteps)
    nn = nnNMIPop.reshape(numUnits*numSteps)

    sortIndex = np.argsort(normSep)
    sortedNormSep = normSep[sortIndex]
    sortedPN = pn[sortIndex]
    sortedNP = nonp[sortIndex]
    sortedPP = pp[sortIndex]
    sortedNN = nn[sortIndex]

    # manual bins - equally populated bins
    equalBinsNormSep = [sortedNormSep[i:i + n] for i in range(0, len(sortedNormSep), n)]
    equalBinsPN = [sortedPN[i:i + n] for i in range(0, len(sortedPN), n)]
    equalBinsNP = [sortedNP[i:i + n] for i in range(0, len(sortedNP), n)]
    equalBinsPP = [sortedPP[i:i + n] for i in range(0, len(sortedPP), n)]
    equalBinsNN = [sortedNN[i:i + n] for i in range(0, len(sortedNN), n)]
    binMeanNormSep = np.array([np.mean(i) for i in equalBinsNormSep])
    binMeanPN = np.array([np.mean(i) for i in equalBinsPN])
    binSEMPN = np.array([np.std(i) for i in equalBinsPN]) / np.sqrt(n)
    binMeanNP = np.array([np.mean(i) for i in equalBinsNP])
    binSEMNP = np.array([np.std(i) for i in equalBinsNP]) / np.sqrt(n)
    binMeanPP = np.array([np.mean(i) for i in equalBinsPP])
    binSEMPP = np.array([np.std(i) for i in equalBinsPP]) / np.sqrt(n)
    binMeanNN = np.array([np.mean(i) for i in equalBinsNN])
    binSEMNN = np.array([np.std(i) for i in equalBinsNN]) / np.sqrt(n)

    ax4.plot(binMeanNormSep, binMeanPN, color='green')
    ax4.errorbar(binMeanNormSep, binMeanPN, yerr=binSEMPN, color='green', ecolor='green',
                 markersize=2, fmt='o')
    ax4.plot(binMeanNormSep, binMeanNP, color='red')
    ax4.errorbar(binMeanNormSep, binMeanNP, yerr=binSEMNP, color='red', ecolor='red',
                 markersize=2, fmt='o')
    ax4.plot(binMeanNormSep, binMeanPP, color='black')
    ax4.errorbar(binMeanNormSep, binMeanPP, yerr=binSEMPP, color='black', ecolor='black',
                 markersize=2, fmt='o')
    ax4.plot(binMeanNormSep, binMeanNN, color='grey')
    ax4.errorbar(binMeanNormSep, binMeanNN, yerr=binSEMNN, color='grey', ecolor='grey',
                 markersize=2, fmt='o')
    ax4.set_ylim(bottom=0)
    ax4.set_xlabel('Gabor Sigma Separation normalized by Eccentricity')
    ax4.set_ylabel('NMI')
    ax4.set_title('NMI, Maunsell and Ni 2017')

    # # raw vs predicted norm PN, NP, PP, NN

    # weighted average
    pCentPred = (np.reshape(prefNormalized[:, 0], (numUnits, 1)) + nonprefNormalized[:, 1:]) / (
                (np.reshape(prefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:]) / (
                 (np.reshape(np.max(transectNormalized, axis=1), (numUnits, 1)))))
    meanPCentPred = (prefMean[0] + nonprefMean[1:]) / (
                    ((prefMean[0] + prefMean[1:])**1) / np.max(transectMean))
    nCentPred = ((np.reshape(nonprefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:])**exp) / (
                (np.reshape(prefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:]) / (
                 (np.reshape(np.max(transectNormalized, axis=1), (numUnits, 1)))))
    meanNPCentPred = ((nonprefMean[0] + prefMean[1:])**exp) / (
                     ((prefMean[0] + prefMean[1:])**1) / np.max(transectMean))
    ppPred = ((np.reshape(prefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:])**exp) / (
                (np.reshape(prefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:]) / (
                 (np.reshape(np.max(transectNormalized, axis=1), (numUnits, 1)))))
    meanPPCentPred = ((prefMean[0] + prefMean[1:])**exp) / (
                     ((prefMean[0] + prefMean[1:])**1) / np.max(transectMean))
    nnPred = ((np.reshape(nonprefNormalized[:, 0], (numUnits, 1)) + nonprefNormalized[:, 1:])**exp) / (
                (np.reshape(prefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:]) / (
                 (np.reshape(np.max(transectNormalized, axis=1), (numUnits, 1)))))
    meanNNCentPred = ((nonprefMean[0] + nonprefMean[1:])**exp) / (
                     ((prefMean[0] + prefMean[1:])**1) / np.max(transectMean))

    # # # simple average prediction
    # pCentPred = ((np.reshape(prefNormalized[:, 0], (numUnits, 1)) + nonprefNormalized[:, 1:])**exp) / 2
    # meanPCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / 2
    # nCentPred = ((np.reshape(nonprefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:])**exp) / 2
    # meanNPCentPred = ((nonprefMean[0] + prefMean[1:])**exp) / 2

    pCentPredReshape = np.reshape(pCentPred, numUnits*4)
    pnNormalizedReshape = np.reshape(pnNormalized, numUnits*4)
    nCentPredReshape = np.reshape(nCentPred, numUnits*4)
    npNormalizedReshape = np.reshape(npNormalized, numUnits*4)
    ax5.scatter(pnNormalizedReshape, pCentPredReshape, color='green',
                label='p center', alpha=0.4)
    ax5.scatter(npNormalizedReshape, nCentPredReshape, color='red',
                label='n center', alpha=0.4)
    ppPredReshape = np.reshape(ppPred, numUnits*4)
    ppNormalizedReshape = np.reshape(ppNormalized, numUnits*4)
    nnPredReshape = np.reshape(nnPred, numUnits*4)
    nnNormalizedReshape = np.reshape(nnNormalized, numUnits*4)
    ax5.scatter(ppNormalizedReshape, ppPredReshape, color='black',
                label='pp', alpha=1)
    ax5.scatter(nnNormalizedReshape, nnPredReshape, color='grey',
                label='nn', alpha=0.4)
    ax5.scatter(pnMean, meanPCentPred, color='green', edgecolor='gold',
                s=80)
    ax5.scatter(npMean, meanNPCentPred, color='red', edgecolor='gold',
                s=80)
    ax5.scatter(ppMean, meanPPCentPred, color='blue', edgecolor='gold',
                s=80)
    ax5.scatter(nnMean, meanNNCentPred, color='blue', edgecolor='gold',
                s=80)
    ax5.set_ylim([0, 1.2])
    ax5.set_xlim([0, 1.2])
    ax5.set_xlabel('Actual Response', fontsize=15, **hfont)
    ax5.set_ylabel('Predicted Response', fontsize=15, **hfont)
    ax5.set_title('actual vs predicated normalization response')
    line = lines.Line2D([0, 1], [0, 1], color='black')
    transform = ax5.transAxes
    line.set_transform(transform)
    ax5.add_line(line)

    # histogram of pred vs real ratio
    predResp = np.concatenate((pCentPredReshape, nCentPredReshape,
                               ppPredReshape, nnPredReshape), axis=0)
    realResp = np.concatenate((pnNormalizedReshape, npNormalizedReshape,
                               ppNormalizedReshape, nnNormalizedReshape), axis=0)
    realPredRatio = realResp / predResp
    ax6.hist(realPredRatio, bins=20)
    ax6.axvline(x=np.median(realPredRatio), color='black',
                label=f'median = {np.median(realPredRatio):.2f}')
    ax6.legend()
    ax6.set_title('Real vs Pred response ratio')
    ax6.set_xlabel('Real Response / Predicted Response')
    ax6.set_ylabel('frequency')

    # offset positions converted to normalized deg separation, plots of responses
    # as a function of distance

    # transect response
    transectNormSep = offsetDegSepNormPop.reshape(numUnits*(numSteps*2+1))
    transectResp = transectNormalized.reshape(numUnits*(numSteps*2+1))
    sortTransectIndex = np.argsort(transectNormSep)
    sortedTransectSep = transectNormSep[sortTransectIndex]
    sortedTransectResp = transectResp[sortTransectIndex]
    equalBinsTransectSep = [sortedTransectSep[i:i + n] for i in range(0, len(sortedTransectSep), n)]
    equalBinsTransectResp = [sortedTransectResp[i:i + n] for i in range(0, len(sortedTransectResp), n)]
    binMeanTransectResp = np.array([np.mean(i) for i in equalBinsTransectResp])
    binMeanTransectSep = np.array([np.mean(i) for i in equalBinsTransectSep])

    transParams = gaussFit(binMeanTransectSep, binMeanTransectResp)
    transXFull = np.linspace(binMeanTransectSep[0], binMeanTransectSep[-1], 1000)
    transFitResp = gauss(transXFull, *transParams)

    # pref/nonpref response
    normSep = offsetDegSepNormPop[:, 4:].reshape(numUnits*(numSteps+1))
    prefResp = prefNormalized.reshape(numUnits*(numSteps+1))
    nonprefResp = nonprefNormalized.reshape(numUnits*(numSteps+1))
    sortIndex = np.argsort(normSep)
    sortedNormSep = normSep[sortIndex]
    sortedPrefResp = prefResp[sortIndex]
    sortedNonprefResp = nonprefResp[sortIndex]

    # manual bins - equally populated bins
    equalBinsNormSep = [sortedNormSep[i:i + n] for i in range(0, len(sortedNormSep), n)]
    equalBinsPrefResp = [sortedPrefResp[i:i + n] for i in range(0, len(sortedPrefResp), n)]
    equalBinsNonprefResp = [sortedNonprefResp[i:i + n] for i in range(0, len(sortedNonprefResp), n)]
    binMeanNormSep = np.array([np.mean(i) for i in equalBinsNormSep])
    binMeanPrefResp = np.array([np.mean(i) for i in equalBinsPrefResp])
    binSEMPrefResp = np.array([np.std(i) for i in equalBinsPrefResp]) / np.sqrt(n)
    binMeanNonprefResp = np.array([np.mean(i) for i in equalBinsNonprefResp])
    binSEMNonpreResp = np.array([np.std(i) for i in equalBinsNonprefResp]) / np.sqrt(n)

    # bin weights
    binWeights = binMeanPrefResp / binMeanPrefResp[0]

    # pn/np/pp/nn response
    sep = offsetDegSepNormPop[:, 5:].reshape(numUnits*numSteps)
    pnResp = pnNormalized.reshape(numUnits*numSteps)
    npResp = npNormalized.reshape(numUnits*numSteps)
    ppResp = ppNormalized.reshape(numUnits*numSteps)
    nnResp = nnNormalized.reshape(numUnits*numSteps)
    sortIndex = np.argsort(sep)
    sortedSep = sep[sortIndex]
    sortedPNResp = pnResp[sortIndex]
    sortedNPResp = npResp[sortIndex]
    sortedPPResp = ppResp[sortIndex]
    sortedNNResp = nnResp[sortIndex]
    # manual bins - equally populated bins
    equalBinsSep = [sortedSep[i:i + n] for i in range(0, len(sortedSep), n)]
    equalBinsPNResp = [sortedPNResp[i:i + n] for i in range(0, len(sortedPNResp), n)]
    equalBinsNPResp = [sortedNPResp[i:i + n] for i in range(0, len(sortedNPResp), n)]
    equalBinsPPResp = [sortedPPResp[i:i + n] for i in range(0, len(sortedPPResp), n)]
    equalBinsNNResp = [sortedNNResp[i:i + n] for i in range(0, len(sortedNNResp), n)]
    binMeanSep = np.array([np.mean(i) for i in equalBinsSep])
    binMeanPNResp = np.array([np.mean(i) for i in equalBinsPNResp])
    binSEMPNResp = np.array([np.std(i) for i in equalBinsPNResp]) / np.sqrt(n)
    binMeanNPResp = np.array([np.mean(i) for i in equalBinsNPResp])
    binSEMNPResp = np.array([np.std(i) for i in equalBinsNPResp]) / np.sqrt(n)
    binMeanPPResp = np.array([np.mean(i) for i in equalBinsPPResp])
    binSEMPPResp = np.array([np.std(i) for i in equalBinsPPResp]) / np.sqrt(n)
    binMeanNNResp = np.array([np.mean(i) for i in equalBinsNNResp])
    binSEMNNResp = np.array([np.std(i) for i in equalBinsNNResp]) / np.sqrt(n)

    # # pn/np pred (RF Weight)
    pnPred = (binMeanPrefResp[0] + binMeanNonprefResp[1:]) / (
              binWeights[0] + binWeights[1:])
    npPred = (binMeanNonprefResp[0] + binMeanPrefResp[1:]) / (
              binWeights[0] + binWeights[1:])

    # pnPred = (binMeanPrefResp[0] + (binMeanNonprefResp[1:] * (binMeanPrefResp[1:] / np.max(transFitResp)))) / (
    #     (binMeanPrefResp[0] + binMeanPrefResp[1:]) / np.max(transFitResp))
    # npPred = (binMeanNonprefResp[0] + (binMeanPrefResp[1:] * (binMeanPrefResp[1:] / np.max(transFitResp)))) / (
    #     (binMeanPrefResp[0] + binMeanPrefResp[1:]) / np.max(transFitResp))

    # pnPred = ((gauss(0, *transParams) + binMeanNonprefResp[1:])**exp) / (
    #         (gauss(0, *transParams) + gauss(binMeanSep, *transParams)) / np.max(transFitResp))
    # npPred = ((binMeanNonprefResp[0] + gauss(binMeanSep, *transParams))**exp) / (
    #         (gauss(0, *transParams) + gauss(binMeanSep, *transParams)) / np.max(transFitResp))

    # pn/np pred (simple average)
    # pnPred = ((binMeanPrefResp[0] + binMeanNonprefResp[1:]) / 2)
    # npPred = ((binMeanNonprefResp[0] + binMeanPrefResp[1:]) / 2)

    ax7.plot(binMeanNormSep, binMeanPrefResp, color='black')
    ax7.errorbar(binMeanNormSep, binMeanPrefResp, yerr=binSEMPrefResp,
                 fmt='o', ecolor='black', color='black', markersize=2)
    ax7.plot(binMeanNormSep, binMeanNonprefResp, color='grey')
    ax7.errorbar(binMeanNormSep, binMeanNonprefResp, yerr=binSEMNonpreResp,
                 fmt='o', ecolor='grey', color='grey', markersize=2)
    ax7.plot(binMeanSep, binMeanPNResp, color='green')
    ax7.errorbar(binMeanSep, binMeanPNResp, yerr=binSEMPNResp,
                 fmt='o', ecolor='green', color='green', markersize=2)
    ax7.plot(binMeanSep, binMeanNPResp, color='red')
    ax7.errorbar(binMeanSep, binMeanNPResp, yerr=binSEMNPResp,
                 fmt='o', ecolor='red', color='red', markersize=2)
    ax7.plot(binMeanSep, pnPred, linestyle='--', color='green')
    ax7.plot(binMeanSep, npPred, linestyle='--', color='red')
    ax7.set_xlabel('Gabor Sigma Separation normalized by Eccentricity')
    ax7.set_ylabel('Normalized Binned Resp')
    ax7.set_ylim(bottom=0)
    ax7.set_title(f'Equal sized bins ({n}) of normalized separation vs Response',
                  fontsize=5)
    ax7.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
    ax7.fill_between(x=np.linspace(0, 1, 10), y1=meanSpon+semSpon,
                     y2=meanSpon-semSpon, color='blue', alpha=0.2)

    # pp/nn pred
    ppPred = ((binMeanPrefResp[0] + binMeanPrefResp[1:])**exp) / (
        (binMeanPrefResp[0] + binMeanPrefResp[1:]) / np.max(binMeanTransectResp))
    nnPred = ((binMeanNonprefResp[0] + binMeanNonprefResp[1:])**exp) / (
        (binMeanPrefResp[0] + binMeanPrefResp[1:]) / np.max(binMeanTransectResp))

    ax8.plot(binMeanNormSep, binMeanPrefResp, color='black')
    ax8.errorbar(binMeanNormSep, binMeanPrefResp, yerr=binSEMPrefResp,
                 fmt='o', ecolor='black', color='black', markersize=2)
    ax8.plot(binMeanNormSep, binMeanNonprefResp, color='grey')
    ax8.errorbar(binMeanNormSep, binMeanNonprefResp, yerr=binSEMNonpreResp,
                 fmt='o', ecolor='grey', color='grey', markersize=2)
    ax8.plot(binMeanSep, binMeanPPResp, color='black', linestyle='--')
    ax8.errorbar(binMeanSep, binMeanPPResp, yerr=binSEMPPResp,
                 fmt='o', ecolor='black', color='black', markersize=2)
    ax8.plot(binMeanSep, binMeanNNResp, color='grey', linestyle='--')
    ax8.errorbar(binMeanSep, binMeanNNResp, yerr=binSEMNNResp,
                 fmt='o', ecolor='grey', color='grey', markersize=2)
    ax8.plot(binMeanSep, ppPred, linestyle=':', color='black', alpha=0.7)
    ax8.plot(binMeanSep, nnPred, linestyle=':', color='grey', alpha=0.7)
    ax8.set_xlabel('Gabor Sigma Separation normalized by Eccentricity')
    ax8.set_ylabel('Normalized Binned Resp')
    ax8.set_ylim(bottom=0)
    ax8.set_title(f'Equal sized bins ({n}) of normalized separation vs Response',
                  fontsize=5)
    ax8.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
    ax8.fill_between(x=np.linspace(0, 1, 10), y1=meanSpon+semSpon,
                     y2=meanSpon-semSpon, color='blue', alpha=0.2)


    plt.tight_layout()
    plt.show()


# sub-figures
# Ni and Maunsell 2017 NMI for just the first offset position
data1 = pnNMIPop[:, 0]
data2 = npNMIPop[:, 0]
data3 = ppNMIPop[:, 0]
data4 = nnNMIPop[:, 0]

plt.hist(data1, bins=20, alpha=1, label='pnNMIPop', color='blue')
plt.hist(data2, bins=20, alpha=1, label='npNMIPop', color='blue')
plt.hist(data3, bins=20, alpha=1, label='ppNMIPop', color='blue')
plt.hist(data4, bins=20, alpha=1, label='nnNMIPop', color='blue')
plt.title('Distribution of Data Points at first offset')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
plt.show()


# plot adaptation profiles
meanPrefResp = np.nanmedian(adaptationMat[:numUnits, 0, :], axis=0)
semPrefResp = stats.sem(adaptationMat[:numUnits, 0, :], axis=0, nan_policy='omit')

trials = np.arange(adaptationMat[:numUnits, 0, :].shape[1])
plt.errorbar(trials, meanPrefResp, yerr=semPrefResp, fmt='-o', capsize=5, capthick=2, label='Mean Response with SEM')
plt.xlabel('Block Number')
plt.ylabel('Response')
plt.title('Average Neuron Response with SEM')
plt.legend()

plt.show()


# fit different normalization models to the data
# RF Weight, Heeger, Ni, Bram models

genericNormR2Scores = []
weightedNormR2Scores = []
heegerNormR2Scores = []
weightedNormFewerParamsR2Scores = []
emsNormR2Scores = []
genNormPos1 = []
genNormPos2 = []
genNormPos3 = []
genNormPos4 = []
weightedNormPos1 = []
weightedNormPos2 = []
weightedNormPos3 = []
weightedNormPos4 = []
heegerNormPos1 = []
heegerNormPos2 = []
heegerNormPos3 = []
heegerNormPos4 = []
emsNormPos1 = []
emsNormPos2 = []
emsNormPos3 = []
emsNormPos4 = []
weight3Params = []
weight5Params = []

for i in range(len(prefNormalized)):
    for j in range(1, 5):
        responses, fixedVals = fixedValsForMTND(prefNormalized[i], nonprefNormalized[i],
                                                pnNormalized[i], npNormalized[i],
                                                ppNormalized[i], nnNormalized[i],
                                                sponNormalized[i], j)
        # bram norm model
        pOpt1, pCov1, = curve_fit(genericNorm, fixedVals.T, responses, maxfev=10000000)  # ,
                                # bounds=[[0, 0, 0, 0, 0], [1, 1, 1, 1, 0.10]])
        y_pred1 = genericNorm(fixedVals.T, *pOpt1)
        r2Gen = r2_score(responses.squeeze(), y_pred1)
        genericNormR2Scores.append(r2Gen)

        # heeger norm model
        pOpt2, pCov2, = curve_fit(heegerNorm, fixedVals.T, responses, maxfev=10000000,
                                  bounds=[[0, 0, 0], [np.inf, np.inf, np.inf]])
        y_pred2 = heegerNorm(fixedVals.T, *pOpt2)
        r2Heeger = r2_score(responses.squeeze(), y_pred2)
        heegerNormR2Scores.append(r2Heeger)

        # RF weight model (5 params)
        pOpt3, pCov3, = curve_fit(weightedNorm, fixedVals.T, responses, maxfev=10000000,
                                  bounds=[[0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf, np.inf]])
        y_pred3 = weightedNorm(fixedVals.T, *pOpt3)
        r2Weight = r2_score(responses.squeeze(), y_pred3)
        weightedNormR2Scores.append(r2Weight)
        residuals3 = responses - y_pred3

        # EMS norm model (5 params)
        pOpt5, pCov5, = curve_fit(emsNorm, fixedVals.T, responses, maxfev=10000000,
                                  bounds=[[-np.inf, -np.inf, -10, -10, -np.inf], [np.inf, np.inf, 10, 10, np.inf]])
        y_pred5 = emsNorm(fixedVals.T, *pOpt5)
        r2EMS = r2_score(responses.squeeze(), y_pred5)
        emsNormR2Scores.append(r2EMS)


        # RF weight model fewer params (3 params)
        responses, fixedVals = fixedValsForRFWeightFewerParams(prefNormalized[i], nonprefNormalized[i],
                                                               pnNormalized[i], npNormalized[i],
                                                               ppNormalized[i], nnNormalized[i],
                                                               sponNormalized[i], j)

        pOpt4, pCov4, = curve_fit(weightedNormFewerParams, fixedVals.T, responses, maxfev=10000000,
                                bounds=[[0, 0, 0], [1, 1, 0.10]])
        y_pred4 = weightedNormFewerParams(fixedVals.T, *pOpt4)
        r2WeightFewerParams = r2_score(responses.squeeze(), y_pred4)
        weightedNormFewerParamsR2Scores.append(r2WeightFewerParams)
        residuals4 = responses - y_pred4

        if j == 1:
            genNormPos1.append(r2Gen)
            heegerNormPos1.append(r2Heeger)
            weightedNormPos1.append(r2Weight)
            emsNormPos1.append(r2EMS)
        elif j == 2:
            genNormPos2.append(r2Gen)
            heegerNormPos2.append(r2Heeger)
            weightedNormPos2.append(r2Weight)
            emsNormPos2.append(r2EMS)
        elif j == 3:
            genNormPos3.append(r2Gen)
            heegerNormPos3.append(r2Heeger)
            weightedNormPos3.append(r2Weight)
            emsNormPos3.append(r2EMS)
        else:
            genNormPos4.append(r2Gen)
            heegerNormPos4.append(r2Heeger)
            weightedNormPos4.append(r2Weight)
            emsNormPos4.append(r2EMS)

        # adjusted R2 comparing weighted model (3 parameters) vs weighted model
        # (5 parameters)
        # Calculate RSS
        rss3 = np.sum(residuals3 ** 2)
        rss4 = np.sum(residuals4 ** 2)

        # Calculate TSS
        tss = np.sum((responses - np.mean(responses)) ** 2)

        # Calculate R-squared and Adjusted R-squared
        r2_3 = 1 - (rss3 / tss)
        r2_4 = 1 - (rss4 / tss)

        n = len(responses)
        p3 = len(pOpt3)
        p4 = len(pOpt4)

        adjR2_3 = 1 - ((1 - r2_3) * (n - 1) / (n - p3 - 1))
        adjR2_4 = 1 - ((1 - r2_4) * (n - 1) / (n - p4 - 1))

        weight3Params.append(adjR2_4)
        weight5Params.append(adjR2_3)


genericMedian = np.median(genericNormR2Scores)
weightedMedian = np.median(weightedNormR2Scores)
heegerMedian = np.median(heegerNormR2Scores)
weightedFewerParamsMedian = np.median(weightedNormFewerParamsR2Scores)
weight3ParamsMedian = np.median(weight3Params)
weight5ParamsMedian = np.median(weight5Params)
emsNormMedian = np.median(emsNormR2Scores)

# figure 0 comparing heeger model to rf weighted model
for i in range(1):
    # Perform the Mann-Whitney U test (Wilcoxon rank-sum test)
    statistic, p_value = mannwhitneyu(heegerNormR2Scores, weightedNormR2Scores, alternative='two-sided')

    # figure
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # Scatter plot
    axes[0, 0].scatter(heegerNormR2Scores, weightedNormR2Scores, alpha=0.7)
    axes[0, 0].set_title('Scatter Plot of R2 Scores')
    axes[0, 0].set_xlabel('heegerNormR2Scores')
    axes[0, 0].set_ylabel('weightedNormR2Scores')
    axes[0, 0].text(0.05, 0.95, f'Median generic: {heegerMedian:.2f}\nMedian weighted: {weightedMedian:.2f}',
                    transform=axes[0, 0].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))
    min_val = min(min(heegerNormR2Scores), min(weightedNormR2Scores))
    max_val = max(max(heegerNormR2Scores), max(weightedNormR2Scores))
    axes[0, 0].plot([min_val, max_val], [min_val, max_val], linestyle='--', color='grey')

    # Histogram of genericNormR2Scores
    axes[0, 1].hist(heegerNormR2Scores, bins=30, alpha=0.7, color='blue')
    axes[0, 1].set_title('Histogram of heegerNormR2Scores')
    axes[0, 1].set_xlabel('R2 Score')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].text(0.95, 0.95, f'Median: {heegerMedian:.2f}',
                    transform=axes[0, 1].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

    # Histogram of weightedNormR2Scores
    axes[1, 0].hist(weightedNormR2Scores, bins=30, alpha=0.7, color='red')
    axes[1, 0].set_title('Histogram of weightedNormR2Scores')
    axes[1, 0].set_xlabel('R2 Score')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].text(0.95, 0.95, f'Median: {weightedMedian:.2f}',
                    transform=axes[1, 0].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

    # Combined histogram
    axes[1, 1].hist(heegerNormR2Scores, bins=30, alpha=0.5, color='blue')
    axes[1, 1].hist(weightedNormR2Scores, bins=30, alpha=0.5, color='red')
    axes[1, 1].set_title('Combined Histogram of R2 Scores')
    axes[1, 1].set_xlabel('R2 Score')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].axvline(heegerMedian, color='blue', linestyle='dashed', linewidth=1)
    axes[1, 1].axvline(weightedMedian, color='red', linestyle='dashed', linewidth=1)
    axes[1, 1].text(0.05, 0.95, f'Mann-Whitney U Test\nStatistic: {statistic:.2f}\nP-value: {p_value:.2e}',
                    transform=axes[1, 1].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

    # Add a star between the vertical lines if distributions are significantly different
    if p_value < 0.05:
        y_max = max(np.histogram(heegerNormR2Scores, bins=30)[0].max(),
                    np.histogram(weightedNormR2Scores, bins=30)[0].max())
        y_star = y_max * 1.1  # Position the star a little above the highest bin
        x_star = (heegerMedian + weightedMedian) / 2  # Position the star between the medians
        axes[1, 1].plot(x_star, y_star, 'k*', markersize=15)
        axes[1, 1].text(x_star, y_star, ' * ', fontsize=12, ha='center', va='bottom')

    plt.tight_layout()
    plt.show()


# figure 1 comparing bram model to rf weighted model
for i in range(1):

    # Perform the Mann-Whitney U test (Wilcoxon rank-sum test)
    statistic, p_value = mannwhitneyu(genericNormR2Scores, weightedNormR2Scores, alternative='two-sided')

    # figure
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # Scatter plot
    axes[0, 0].scatter(genericNormR2Scores, weightedNormR2Scores, alpha=0.7)
    axes[0, 0].set_title('Scatter Plot of R2 Scores')
    axes[0, 0].set_xlabel('genericNormR2Scores')
    axes[0, 0].set_ylabel('weightedNormR2Scores')
    axes[0, 0].text(0.05, 0.95, f'Median generic: {genericMedian:.2f}\nMedian weighted: {weightedMedian:.2f}',
                    transform=axes[0, 0].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))
    min_val = min(min(genericNormR2Scores), min(weightedNormR2Scores))
    max_val = max(max(genericNormR2Scores), max(weightedNormR2Scores))
    axes[0, 0].plot([min_val, max_val], [min_val, max_val], linestyle='--', color='grey')

    # Histogram of genericNormR2Scores
    axes[0, 1].hist(genericNormR2Scores, bins=30, alpha=0.7, color='blue')
    axes[0, 1].set_title('Histogram of genericNormR2Scores')
    axes[0, 1].set_xlabel('R2 Score')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].text(0.95, 0.95, f'Median: {genericMedian:.2f}',
                    transform=axes[0, 1].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

    # Histogram of weightedNormR2Scores
    axes[1, 0].hist(weightedNormR2Scores, bins=30, alpha=0.7, color='red')
    axes[1, 0].set_title('Histogram of weightedNormR2Scores')
    axes[1, 0].set_xlabel('R2 Score')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].text(0.95, 0.95, f'Median: {weightedMedian:.2f}',
                    transform=axes[1, 0].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

    # Combined histogram
    axes[1, 1].hist(genericNormR2Scores, bins=30, alpha=0.5, color='blue')
    axes[1, 1].hist(weightedNormR2Scores, bins=30, alpha=0.5, color='red')
    axes[1, 1].set_title('Combined Histogram of R2 Scores')
    axes[1, 1].set_xlabel('R2 Score')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].axvline(genericMedian, color='blue', linestyle='dashed', linewidth=1)
    axes[1, 1].axvline(weightedMedian, color='red', linestyle='dashed', linewidth=1)
    axes[1, 1].text(0.05, 0.95, f'Mann-Whitney U Test\nStatistic: {statistic:.2f}\nP-value: {p_value:.2e}',
                    transform=axes[1, 1].transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

    # Add a star between the vertical lines if distributions are significantly different
    if p_value < 0.05:
        y_max = max(np.histogram(genericNormR2Scores, bins=30)[0].max(),
                    np.histogram(weightedNormR2Scores, bins=30)[0].max())
        y_star = y_max * 1.1  # Position the star a little above the highest bin
        x_star = (genericMedian + weightedMedian) / 2  # Position the star between the medians
        axes[1, 1].plot(x_star, y_star, 'k*', markersize=15)
        axes[1, 1].text(x_star, y_star, ' * ', fontsize=12, ha='center', va='bottom')

    plt.tight_layout()
    plt.show()


# figure 2 violin/box plots of all 4 models and their performance at different offset points
for i in range(1):
    genNormData = [genNormPos1, genNormPos2, genNormPos3, genNormPos4]
    heegerNormData = [heegerNormPos1, heegerNormPos2, heegerNormPos3, heegerNormPos4]
    weightedNormData = [weightedNormPos1, weightedNormPos2, weightedNormPos3, weightedNormPos4]
    emsNormData = [emsNormPos1, emsNormPos2, emsNormPos3, emsNormPos4]

    fig, axes = plt.subplots(1, 4, figsize=(15, 4))

    # Plot the violin plots for genNorm arrays
    sns.boxplot(data=genNormData, ax=axes[0])
    axes[0].set_title('Box Plot for bramNorm Arrays')
    axes[0].set_xlabel('Position')
    axes[0].set_ylabel('Values')
    axes[0].set_xticks([0, 1, 2, 3])
    axes[0].set_xticklabels(['Pos1', 'Pos2', 'Pos3', 'Pos4'])
    axes[0].set_ylim(0, 1.2)

    # Calculate and annotate median values
    medians = [np.median(data) for data in genNormData]
    maxY = max(max(data) for data in genNormData)
    yText = maxY + 0.1
    for j, median in enumerate(medians):
        axes[0].text(j, yText, f'{median:.2f}', ha='center', va='bottom')

    # Plot the violin plots for heegerNorm arrays
    sns.boxplot(data=heegerNormData, ax=axes[1])
    axes[1].set_title('Box Plot for heegerNorm Arrays')
    axes[1].set_xlabel('Position')
    axes[1].set_ylabel('Values')
    axes[1].set_xticks([0, 1, 2, 3])
    axes[1].set_xticklabels(['Pos1', 'Pos2', 'Pos3', 'Pos4'])
    axes[1].set_ylim(0, 1.2)

    # Calculate and annotate median values
    medians = [np.median(data) for data in heegerNormData]
    maxY = max(max(data) for data in heegerNormData)
    yText = maxY + 0.1
    for j, median in enumerate(medians):
        axes[1].text(j, yText, f'{median:.2f}', ha='center', va='bottom')

    # Plot the violin plots for emsNorm arrays
    sns.boxplot(data=emsNormData, ax=axes[2])
    axes[2].set_title('Box Plot for emsNorm Arrays')
    axes[2].set_xlabel('Position')
    axes[2].set_ylabel('Values')
    axes[2].set_xticks([0, 1, 2, 3])
    axes[2].set_xticklabels(['Pos1', 'Pos2', 'Pos3', 'Pos4'])
    axes[2].set_ylim(0, 1.2)

    # Calculate and annotate median values
    medians = [np.median(data) for data in emsNormData]
    maxY = max(max(data) for data in emsNormData)
    yText = maxY + 0.1
    for j, median in enumerate(medians):
        axes[2].text(j, yText, f'{median:.2f}', ha='center', va='bottom')

    # Plot the violin plots for weightedNorm arrays
    sns.boxplot(data=weightedNormData, ax=axes[3])
    axes[3].set_title('Box Plot for weightedNorm Arrays')
    axes[3].set_xlabel('Position')
    axes[3].set_ylabel('Values')
    axes[3].set_xticks([0, 1, 2, 3])
    axes[3].set_xticklabels(['Pos1', 'Pos2', 'Pos3', 'Pos4'])
    axes[3].set_ylim(0, 1.2)

    # Calculate and annotate median values
    medians = [np.median(data) for data in weightedNormData]
    maxY = max(max(data) for data in weightedNormData)
    yText = maxY + 0.1
    for j, median in enumerate(medians):
        axes[3].text(j, yText, f'{median:.2f}', ha='center', va='bottom')

    plt.tight_layout()
    plt.show()
