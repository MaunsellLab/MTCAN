"""
MTND main analysis script. This script will extract spike counts during
valid trials from units active in a session and put them in the appropriate
matrix. From there, this script plots how normalization changes with distance.
We also plot the preferred response across the diameter of the receptive field

Chery - June 2023

line 1465 : population plots
"""
import matplotlib.pyplot as plt

# Imports
from usefulFns import *

# subpar sessions: 230606, 230706, 230712

p0 = time.time()


############# MEETZ ##########

# subpar sessions: 230606, 230706, 230712
# fileList = ['230719']

# # master file list (MEETZ)
# fileList = ['230607', '230608', '230612', '230613', '230615',
#             '230620', '230626', '230627', '230628', '230630',
#             '230707', '230710', '230711', '230718', '230719',
#             '230720']
#
# # # master unit list
# # unitList = ['230607_147', '230607_146', '230607_144', '230607_140', '230607_126',
# #             '230607_125', '230608_72', '230608_137', '230608_143', '230608_145',
# #             '230608_146', '230608_149', '230612_213', '230613_187', '230615_83',
# #             '230615_164', '230615_166', '230620_58', '230626_131', '230627_147',
# #             '230627_141', '230627_106', '230627_71', '230627_70', '230627_14',
# #             '230628_143', '230630_50', '230707_140', '230710_137', '230711_45',
# #             '230711_50', '230718_178', '230719_147', '230719_156', '230720_149']
#
# # master unit list with null resp > baseline
# unitList = ['230607_147', '230607_146', '230607_144', '230607_126','230607_125',
#             '230608_143', '230608_145', '230608_149', '230612_213', '230613_187',
#             '230615_83', '230615_166', '230620_58', '230626_131', '230627_147',
#             '230627_141', '230627_106', '230627_71', '230627_70', '230627_14',
#             '230628_143', '230630_50', '230707_140', '230710_137', '230711_45',
#             '230711_50', '230718_178', '230719_147', '230719_156', '230720_149']


# # pref + null (MEETZ)
# fileList = ['230607', '230608', '230612', '230613', '230615']
#
# unitList = ['230607_147', '230607_146', '230607_144', '230607_140', '230607_126',
#             '230607_125', '230608_72', '230608_137', '230608_143', '230608_145',
#             '230608_146', '230608_149', '230612_213', '230613_187', '230615_83',
#             '230615_164', '230615_166']
#
# # pref + null units with null above baseline
# unitList = ['230607_147', '230607_146', '230607_144', '230607_126','230607_125',
#             '230608_143', '230608_145', '230608_149', '230612_213', '230613_187',
#             '230615_83', '230615_166']

# pref + nonpref (MEETZ)
# fileList = ['230620', '230626', '230627', '230628', '230630',
#             '230707', '230710', '230711', '230718', '230719',
#             '230720']
#
# unitList = ['230620_58', '230626_131', '230627_147', '230627_141', '230627_106',
#             '230627_71', '230627_70', '230627_14', '230628_143', '230630_50',
#             '230707_140', '230710_137', '230711_45', '230711_50', '230718_178',
#             '230719_147', '230719_156', '230720_149']


####### AKSHAN #######

# pref + nonpref (AKSHAN) (potential good units/sessions)
# fileList = ['240529', '240530', '240603', '240605', '240606', '240607' ,
#             '240610', '240611', '240613']
# unitList = ['240529_31', '240530_55', '240603_167', '240605_147', '240606_176', '240607_137', '240607_155',
#             '240607_164', '240607_170', '240607_179', '240610_169']
'''
# 240607, check units some actually prefer the 'non-pref'
# 240611, lots of work to clean up data
'''
# cherry-picked good sessions/units for population plots
# fileList = ['240530', '240603', '240606', '240607', '240611', '240613', '240628']
# unitList = ['240530_55', '240603_167', '240606_176', '240607_137', '240607_170', '240607_179',
#             '240611_136', '240613_138', '240628_19']

# list if I consider poorly centered stimuli
# fileList = ['240530', '240603', '240606', '240607', '240611', '240613', '240703']
# unitList = ['240530_55', '240603_167', '240606_176', '240607_137', '240607_170', '240607_179',
#             '240611_62', '240611_67', '240611_80', '240611_131', '240611_132', '240611_135',
#             '240611_136', '240611_140', '240611_145', '240611_151', '240613_138', '240703_91']

fileList = ['240709']
unitList = ['240709_7']


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
adaptationMat = np.zeros((100, 4, 100))
adaptationMat[:] = np.nan
adaptC = 0
allBlocksDone = []

for file in fileList:
    # Load relevant file here with pyMat reader
    monkeyName = 'Akshan'  # 'Meetz'
    seshDate = file
    withinSesh = 1
    if withinSesh == 1:
        fileName = f'{monkeyName}_{seshDate}_MTND_Spikes.mat'
    else:
        fileName = f'{monkeyName}_{seshDate}_MTND_2_Spikes.mat'
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

    if withinSesh == 1:
        if not os.path.exists('Normalization'):
            os.makedirs('Normalization')
        os.chdir('Normalization/')
    else:
        if not os.path.exists('Normalization Unit 2'):
            os.makedirs('Normalization Unit 2')
        os.chdir('Normalization Unit 2/')

    # list of indices of correctTrials (non-instruct, valid trialCertify)
    corrTrials = correctTrialsMTX(allTrials)

    # generate list of unique active units, and their channel
    units = activeUnits('spikeData', allTrials)
    unitCluster = allTrials[corrTrials[0]]['spikeTempInfo']['cgs']
    unitsChannel = unitsInfo(units, corrTrials, allTrials)

    # # adds sessions unit to master list of units across sessions
    # for unit in units:
    #     unitIdentifier.append(f'{fileIterator}_{unit}')

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
    allBlocksDone.append(blocksDone)
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

                    for unitCount, unit in enumerate(units):
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

    # direction tuning arrays for each unit
    allDirTuning = np.load('../Direction Tuning/unitsDirTuningMat.npy')
    allDirTuningSEM = np.load('../Direction Tuning/unitsDirTuningSEM.npy')

    # RF location arrays for each unit
    allRFLoc = np.load('../RFLoc Tuning/unitsRFLocMat.npy')
    eleLabel = np.load('../RFLoc Tuning/eleLabels.npy')
    aziLabel = np.load('../RFLoc Tuning/aziLabels.npy')

    # # correlations with distance
    # corrUnits = []
    # for unit in units:
    #     if f'{seshDate}_{unit}' in unitList:
    #         corrUnits.append(unit)
    # if len(corrUnits) > 1:
    #     # offset1Corr = []
    #     # offset2Corr = []
    #     # offset3Corr = []
    #     # offset4Corr = []
    #     combs = [i for i in combinations(corrUnits, 2)]
    #
    #     for comb in combs:
    #         n1, n2 = comb[0], comb[1]
    #         n1Index = np.where(units == n1)[0][0]
    #         n2Index = np.where(units == n2)[0][0]
    #
    #         tempCorr = []
    #         tempSingleCorr = []
    #         tempPairedCorr = []
    #         tempRFWeight = []
    #         tempRFWeightPaired = []
    #         tempRFWeightSingle = []
    #         # spontaneous correlation
    #         n1SpikeMat = spikeCountMat[n1Index, :blocksDone, 0]
    #         n2SpikeMat = spikeCountMat[n2Index, :blocksDone, 0]
    #         pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)
    #         sponCorrPop.append(pairStimCorr)
    #
    #         # n1 center response
    #         n1CenterP = meanSpikeReshaped[n1Index, 2, 0]
    #         n1CenterNP = meanSpikeReshaped[n1Index, 1, 0]
    #
    #         # n2 center response
    #         n2CenterP = meanSpikeReshaped[n2Index, 2, 0]
    #         n2CenterNP = meanSpikeReshaped[n2Index, 1, 0]
    #
    #         for i in range(1, 5):
    #             stimIndices = stimCountIndex[1:, 2*i-1:2*i+1].reshape(4)
    #             x = np.linspace(-4, 4, 9)
    #             # n1 pref transect
    #             transReverse = transectMeanSpike[n1Index, ::-1]
    #             prefTransect = np.concatenate((transReverse,
    #                                           [meanSpikeReshaped[n1Index, 2, 0]],
    #                                           meanSpikeReshaped[n1Index, 0, 2::2]),
    #                                           axis=0)
    #             params = gaussFit(x, prefTransect)
    #             xFull = np.linspace(-4, 4, 1000)
    #             respFull = gauss(xFull, *params)
    #             rfWeightN1Cent = gauss(0, *params) / np.max(respFull)
    #             rfWeightN1 = gauss(i, *params) / np.max(respFull)
    #             combWeightN1 = (rfWeightN1Cent + rfWeightN1) / 2
    #             # rfWeightN1 = prefTransect[4+i] / np.max(prefTransect)
    #
    #             # n2 pref transect
    #             transReverse = transectMeanSpike[n2Index, ::-1]
    #             prefTransect = np.concatenate((transReverse,
    #                                           [meanSpikeReshaped[n2Index, 2, 0]],
    #                                           meanSpikeReshaped[n2Index, 0, 2::2]),
    #                                           axis=0)
    #             params = gaussFit(x, prefTransect)
    #             respFull = gauss(xFull, *params)
    #             rfWeightN2Cent = gauss(0, *params) / np.max(respFull)
    #             rfWeightN2 = gauss(i, *params) / np.max(respFull)
    #             combWeightN2 = (rfWeightN2Cent + rfWeightN2) / 2
    #             # rfWeightN2 = prefTransect[4+i] / np.max(prefTransect)
    #
    #             # rfWeight = np.sqrt(combWeightN1 * combWeightN2)
    #             rfWeight = np.power(rfWeightN1Cent*rfWeightN1*rfWeightN2Cent*rfWeightN2, (1/4))
    #
    #             # n1 offset response
    #             n1OffsetP = meanSpikeReshaped[n1Index, 0, i*2]
    #             n1OffsetNP = meanSpikeReshaped[n1Index, 0, i*2-1]
    #
    #             # n2 center response
    #             n2OffsetP = meanSpikeReshaped[n2Index, 0, i*2]
    #             n2OffsetNP = meanSpikeReshaped[n2Index, 0, i*2-1]
    #
    #             # paired stimulus correlations for offset
    #             for count, j in enumerate(stimIndices):
    #                 # correlation for that Gabor pair b/w 2 units excluding trials where
    #                 # spike counts exceeded 3 SD from mean
    #                 n1SpikeMat = spikeCountMat[n1Index, :blocksDone, j]
    #                 n2SpikeMat = spikeCountMat[n2Index, :blocksDone, j]
    #                 pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)
    #                 # totCorrPop.append(pairStimCorr)
    #                 # rfWeightPop.append(rfWeight)
    #                 tempCorr.append(pairStimCorr)
    #                 tempPairedCorr.append(pairStimCorr)
    #                 tempRFWeight.append(rfWeight)
    #                 tempRFWeightPaired.append(rfWeight)
    #
    #                 # n1/n2 condition resp normalized to max resp
    #                 n1NormalizedResp = meanSpike[n1Index, j] / np.max(meanSpike[n1Index])
    #                 n2NormalizedResp = meanSpike[n2Index, j] / np.max(meanSpike[n2Index])
    #                 normalizedRespWeightPop.append((n1NormalizedResp + n2NormalizedResp)/2)
    #
    #                 if i == 1:
    #                     # offset1Corr.append(pairStimCorr)
    #                     off1CorrPop.append(pairStimCorr)
    #                 elif i == 2:
    #                     # offset2Corr.append(pairStimCorr)
    #                     off2CorrPop.append(pairStimCorr)
    #                 elif i == 3:
    #                     # offset3Corr.append(pairStimCorr)
    #                     off3CorrPop.append(pairStimCorr)
    #                 elif i == 4:
    #                     # offset4Corr.append(pairStimCorr)
    #                     off4CorrPop.append(pairStimCorr)
    #                 if count == 0:
    #                     n1NMI = (n1CenterNP + n1OffsetNP - meanSpikeReshaped[n1Index, 1, i*2-1]) / (
    #                              n1CenterNP + n1OffsetNP + meanSpikeReshaped[n1Index, 1, i*2-1])
    #                     n2NMI = (n2CenterNP + n2OffsetNP - meanSpikeReshaped[n2Index, 1, i*2-1]) / (
    #                             n2CenterNP + n2OffsetNP + meanSpikeReshaped[n2Index, 1, i*2-1])
    #                     pairNMI = (n1NMI + n2NMI) / 2
    #                     pairNMIPop.append(pairNMI)
    #                 elif count == 1:
    #                     n1NMI = (n1CenterNP + n1OffsetP - meanSpikeReshaped[n1Index, 1, i*2]) / (
    #                              n1CenterNP + n1OffsetP + meanSpikeReshaped[n1Index, 1, i*2])
    #                     n2NMI = (n2CenterNP + n2OffsetP - meanSpikeReshaped[n2Index, 1, i*2]) / (
    #                             n2CenterNP + n2OffsetP + meanSpikeReshaped[n2Index, 1, i*2])
    #                     pairNMI = (n1NMI + n2NMI) / 2
    #                     pairNMIPop.append(pairNMI)
    #                 elif count == 2:
    #                     n1NMI = (n1CenterP + n1OffsetNP - meanSpikeReshaped[n1Index, 2, i*2-1]) / (
    #                              n1CenterP + n1OffsetNP + meanSpikeReshaped[n1Index, 2, i*2-1])
    #                     n2NMI = (n2CenterP + n2OffsetNP - meanSpikeReshaped[n2Index, 2, i*2-1]) / (
    #                             n2CenterP + n2OffsetNP + meanSpikeReshaped[n2Index, 2, i*2-1])
    #                     pairNMI = (n1NMI + n2NMI) / 2
    #                     pairNMIPop.append(pairNMI)
    #                 elif count == 3:
    #                     n1NMI = (n1CenterP + n1OffsetP - meanSpikeReshaped[n1Index, 2, i*2]) / (
    #                              n1CenterP + n1OffsetP + meanSpikeReshaped[n1Index, 2, i*2])
    #                     n2NMI = (n2CenterP + n2OffsetP - meanSpikeReshaped[n2Index, 2, i*2]) / (
    #                             n2CenterP + n2OffsetP + meanSpikeReshaped[n2Index, 2, i*2])
    #                     pairNMI = (n1NMI + n2NMI) / 2
    #                     pairNMIPop.append(pairNMI)
    #
    #             # single stimulus correlations for offset
    #             stimIndices = stimCountIndex[0, 2*i-1:2*i+1]
    #             for j in stimIndices:
    #                 # correlation for that Gabor pair b/w 2 units excluding trials where
    #                 # spike counts exceeded 3 SD from mean
    #                 n1SpikeMat = spikeCountMat[n1Index, :blocksDone, j]
    #                 n2SpikeMat = spikeCountMat[n2Index, :blocksDone, j]
    #                 pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)
    #                 # totCorrPop.append(pairStimCorr)
    #                 # rfWeightPop.append(rfWeight)
    #                 tempCorr.append(pairStimCorr)
    #                 tempSingleCorr.append(pairStimCorr)
    #                 tempRFWeight.append(np.sqrt(rfWeightN1 * rfWeightN2))
    #                 tempRFWeightSingle.append(np.sqrt(rfWeightN1 * rfWeightN2))
    #
    #                 # n1/n2 condition resp normalized to max resp
    #                 n1NormalizedResp = meanSpike[n1Index, j] / np.max(meanSpike[n1Index])
    #                 n2NormalizedResp = meanSpike[n2Index, j] / np.max(meanSpike[n2Index])
    #                 normalizedRespWeightPop.append((n1NormalizedResp + n2NormalizedResp)/2)
    #
    #             # single stimulus correlations for center
    #             stimIndices = stimCountIndex[1:, 0]
    #             for j in stimIndices:
    #                 # correlation for that Gabor pair b/w 2 units excluding trials where
    #                 # spike counts exceeded 3 SD from mean
    #                 n1SpikeMat = spikeCountMat[n1Index, :blocksDone, j]
    #                 n2SpikeMat = spikeCountMat[n2Index, :blocksDone, j]
    #                 pairStimCorr, pairDCov, pairDSD = pairCorrExclude3SD(n1SpikeMat, n2SpikeMat)
    #                 # totCorrPop.append(pairStimCorr)
    #                 # rfWeightPop.append(rfWeight)
    #                 tempCorr.append(pairStimCorr)
    #                 tempSingleCorr.append(pairStimCorr)
    #                 tempRFWeight.append(np.sqrt(rfWeightN1Cent * rfWeightN2Cent))
    #                 tempRFWeightSingle.append(np.sqrt(rfWeightN1Cent * rfWeightN2Cent))
    #
    #                 # n1/n2 condition resp normalized to max resp
    #                 n1NormalizedResp = meanSpike[n1Index, j] / np.max(meanSpike[n1Index])
    #                 n2NormalizedResp = meanSpike[n2Index, j] / np.max(meanSpike[n2Index])
    #                 normalizedRespWeightPop.append((n1NormalizedResp + n2NormalizedResp)/2)
    #
    #         tempCorr = np.array(tempCorr)
    #         tempPairedCorr = np.array(tempPairedCorr)
    #         tempSingleCorr = np.array(tempSingleCorr)
    #         normVal = max(abs(sponCorrPop[-1]), np.nanmax(np.abs(tempCorr)))
    #         tempCorr = tempCorr / normVal
    #         tempPairedCorr = tempPairedCorr / normVal
    #         tempSingleCorr = tempSingleCorr / normVal
    #         sponCorrPop[-1] = sponCorrPop[-1] / normVal
    #
    #         for x in range(len(tempCorr)):
    #             totCorrPop.append(tempCorr[x])
    #             rfWeightPop.append(tempRFWeight[x])
    #
    #         for x in range(len(tempPairedCorr)):
    #             pairedCorrPop.append(tempPairedCorr[x])
    #             rfWeightPairedPop.append(tempRFWeightPaired[x])
    #
    #         for x in range(len(tempSingleCorr)):
    #             singleCorrPop.append(tempSingleCorr[x])
    #             rfWeightSinglePop.append(tempRFWeightSingle[x])
    #
    #
    #
    #     # off1CorrMean = np.mean(offset1Corr)
    #     # off2CorrMean = np.mean(offset2Corr)
    #     # off3CorrMean = np.mean(offset3Corr)
    #     # off4CorrMean = np.mean(offset4Corr)
    #     # off1CorrPop.append(off1CorrMean)
    #     # off2CorrPop.append(off2CorrMean)
    #     # off3CorrPop.append(off3CorrMean)
    #     # off4CorrPop.append(off4CorrMean)

    # figure (line plots)
    for count in range(len(units)):
        fig = plt.figure()
        fig.suptitle(f'unit {units[count]}')
        fig.set_size_inches(15, 7)
        gs0 = gridspec.GridSpec(2, 4)
        gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0, 0])
        gs01 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0, 1])
        gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1, 0])
        gs03 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0, 2])
        gs04 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1, 2])
        gs05 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1, 1])
        gs06 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1, 3])
        gs07 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0, 3])

        # polar
        prefDirIndex = np.where(np.arange(0, 390, 360/len(allDirTuning[0])) == prefDir)[0]
        ax = fig.add_subplot(gs00[0, 0], polar=True)
        theta = np.radians(np.arange(0, 390, 360/len(allDirTuning[0])))
        r = (np.append(allDirTuning[count], allDirTuning[count][0]))
        err = (np.append(allDirTuningSEM[count], allDirTuningSEM[count][0]))
        ax.plot(theta, r, markersize=2, color='black')
        ax.errorbar(theta, r, yerr=err, fmt='o', ecolor='black',
                    color='black', markersize=2)
        # ax.errorbar(theta[prefDirIndex], r[prefDirIndex], yerr=err[prefDirIndex],
        #             fmt='o', ecolor='gold', markersize=4, color='gold')
        ax.set_theta_zero_location("W")
        ax.set_title('Direction Tuning Polar Plot', fontsize=8)
        ax.set_ylim([0, np.max(r) * 1.1])
        ax.scatter(np.radians(prefDir), np.max(r)*1.07, color='gold', s=70)
        ax.scatter(np.radians(nonPrefDir), np.max(r)*1.07, color='grey', s=70)
        extraTicks = [np.radians(prefDir), np.radians(nonPrefDir)]
        ax.set_xticks(np.radians(np.linspace(0, 330, 12)))
        ax.set_xticks(list(plt.xticks()[0]) + extraTicks)
        ax.set_rticks(np.linspace(0, np.max(r)*1.1, 3))

        # RF Heatmap
        ax1 = fig.add_subplot(gs03[0, 0])

        dx = abs(aziLabel[0] - aziLabel[1]) / 2
        dy = abs(eleLabel[0] - eleLabel[1]) / 2
        extent = [aziLabel[0] - dx, aziLabel[-1] + dx,
                  eleLabel[0] - dy, eleLabel[-1] + dy]
        ax1.imshow(np.flip(allRFLoc[count], axis=0), extent=extent, vmin=0, cmap='gist_heat',
                   aspect='auto')
        ax1.set_xlabel('azimuth (˚)', fontsize=8)
        ax1.set_ylabel('elevation (˚)', fontsize=8)
        ax1.set_title('Heatmap of unit RF location', fontsize=9)
        for offCount, i in enumerate(offsetAziEle):
            if offCount == 0:
                col = 'blue'
            elif offCount > numSteps:
                col = 'green'
            else:
                col = 'gold'
            ax1.scatter(i[0], i[1], s=50, marker='o', color=col)

        # pref response across diameter of RF
        ax3 = fig.add_subplot(gs01[0, 0])
        x = np.concatenate((np.arange(-numSteps, 0),
                            np.arange(0, numSteps+1)),
                           axis=0)
        transReverse = transectMeanSpike[count, ::-1]
        transSEMReverse = transectSEM[count, ::-1]
        prefTransect = np.concatenate((transReverse,
                                       [meanSpikeReshaped[count, 2, 0]],
                                       meanSpikeReshaped[count, 0, 2::2]),
                                      axis=0)
        prefTransectSEM = np.concatenate((transSEMReverse,
                                          [SEMReshaped[count, 2, 0]],
                                          SEMReshaped[count, 0, 2::2]),
                                         axis=0)
        spon = meanSpikeReshaped[count, 0, 0]
        ax3.plot(x, prefTransect, color='black')
        ax3.errorbar(x, prefTransect, yerr=prefTransectSEM, fmt='o', ecolor='black',
                     color='black', markersize=2)
        ax3.set_xlabel('stimulus offset positions')
        ax3.set_ylabel('Response (spikes/sec)')
        ax3.set_ylim(bottom=0)
        ax3.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)

        # normalization vs distance (PN, NP)
        ax2 = fig.add_subplot(gs02[0, 0])
        x = np.arange(0, numSteps+1)
        xNorm = np.arange(1, numSteps+1)
        nullOnly = np.concatenate(([meanSpikeReshaped[count, 1, 0]],
                                  meanSpikeReshaped[count, 0, 1::2]), axis=0)
        nullSEM = np.concatenate(([SEMReshaped[count, 1, 0]],
                                  SEMReshaped[count, 0, 1::2]), axis=0)
        prefOnly = np.concatenate(([meanSpikeReshaped[count, 2, 0]],
                                  meanSpikeReshaped[count, 0, 2::2]), axis=0)
        prefSEM = np.concatenate(([SEMReshaped[count, 2, 0]],
                                  SEMReshaped[count, 0, 2::2]), axis=0)
        # RF weighted average
        pCentPred = (prefOnly[0] + (nullOnly[1:] / np.max(prefTransect))) / (
                    (prefOnly[0] + prefOnly[1:]) / np.max(prefTransect))
        npCentPred = (nullOnly[0] + (prefOnly[1:] / np.max(prefTransect))) / (
                     (prefOnly[0] + prefOnly[1:]) / np.max(prefTransect))
        # # simple average
        # pCentPred = (prefOnly[0] + nullOnly[1:]) / 2
        # npCentPred = (nullOnly[0] + prefOnly[1:]) / 2
        NP = meanSpikeReshaped[count, 1, 2::2]
        npSEM = SEMReshaped[count, 1, 2::2]
        PN = meanSpikeReshaped[count, 2, 1::2]
        pnSEM = SEMReshaped[count, 2, 1::2]
        ax2.plot(x, prefOnly, color='black', label='P0')
        ax2.errorbar(x, prefOnly, yerr=prefSEM, fmt='o', ecolor='black',
                     color='black', markersize=2)
        ax2.plot(x, nullOnly, color='grey', label='N0')
        ax2.errorbar(x, nullOnly, yerr=nullSEM, fmt='o', ecolor='grey',
                     color='grey', markersize=2)
        ax2.plot(xNorm, NP, color='red', label='N0 P1')
        ax2.errorbar(xNorm, NP, yerr=npSEM, fmt='o', ecolor='red',
                     color='red', markersize=2)
        ax2.plot(xNorm, PN, color='green', label='P0 N1')
        ax2.errorbar(xNorm, PN, yerr=pnSEM, fmt='o', ecolor='green',
                     color='green', markersize=2)
        ax2.plot(xNorm, pCentPred, color='green', linestyle='--', alpha=0.8)
        ax2.plot(xNorm, npCentPred, color='red', linestyle='--', alpha=0.8)
        ax2.set_xlabel('stimulus offset positions')
        ax2.set_ylabel('Response (spikes/sec)')
        ax2.set_ylim(bottom=0, top=np.max(meanSpikeReshaped[count])*1.1)
        ax2.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
        ax2.legend(fontsize=4)

        # normalization vs distance (PP, NN)
        ax5 = fig.add_subplot(gs05[0, 0])
        x = np.arange(0, numSteps+1)
        xNorm = np.arange(1, numSteps+1)
        spon = meanSpikeReshaped[count, 0, 0]
        nullOnly = np.concatenate(([meanSpikeReshaped[count, 1, 0]],
                                  meanSpikeReshaped[count, 0, 1::2]), axis=0)
        nullSEM = np.concatenate(([SEMReshaped[count, 1, 0]],
                                  SEMReshaped[count, 0, 1::2]), axis=0)
        prefOnly = np.concatenate(([meanSpikeReshaped[count, 2, 0]],
                                  meanSpikeReshaped[count, 0, 2::2]), axis=0)
        prefSEM = np.concatenate(([SEMReshaped[count, 2, 0]],
                                  SEMReshaped[count, 0, 2::2]), axis=0)
        PP = meanSpikeReshaped[count, 2, 2::2]
        ppSEM = SEMReshaped[count, 2, 2::2]
        NN = meanSpikeReshaped[count, 1, 1::2]
        nnSEM = SEMReshaped[count, 1, 1::2]
        ax5.plot(x, prefOnly, color='black', label='P0')
        ax5.errorbar(x, prefOnly, yerr=prefSEM, fmt='o', ecolor='black',
                     color='black', markersize=2)
        ax5.plot(x, nullOnly, color='grey', label='N0')
        ax5.errorbar(x, nullOnly, yerr=nullSEM, fmt='o', ecolor='grey',
                     color='grey', markersize=2)
        ax5.plot(xNorm, PP, color='black', linestyle='--', label='PO P1', alpha=0.8)
        ax5.errorbar(xNorm, PP, yerr=ppSEM, fmt='o', ecolor='black',
                     color='black', markersize=2)
        ax5.plot(xNorm, NN, color='grey', linestyle='--', label='N0 N1', alpha=0.7)
        ax5.errorbar(xNorm, NN, yerr=nnSEM, fmt='o', ecolor='grey',
                     color='grey', markersize=2)
        ax5.set_xlabel('stimulus offset positions')
        ax5.set_ylabel('Response (spikes/sec)')
        ax5.set_ylim(bottom=0, top=np.max(meanSpikeReshaped[count])*1.1)
        ax5.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
        ax5.legend(fontsize=4)

        # NMI vs distance (Ruff and Cohen, 2016)
        ax4 = fig.add_subplot(gs04[0, 0])
        x = np.arange(1, numSteps + 1)
        pCentNMI = []
        nCentNMI = []
        ppNMI = []
        nnNMI = []
        for i in np.arange(1, numSteps+1):
            # PN
            centerP = meanSpikeReshaped[count, 2, 0]
            offsetN = meanSpikeReshaped[count, 0, i*2-1]
            PN = meanSpikeReshaped[count, 2, i*2-1]
            pnNMI = (centerP + offsetN) / PN
            pCentNMI.append(pnNMI)

            # NP
            centerN = meanSpikeReshaped[count, 1, 0]
            offsetP = meanSpikeReshaped[count, 0, i*2]
            NP = meanSpikeReshaped[count, 1, i*2]
            npNMI = (centerN + offsetP) / NP
            nCentNMI.append(npNMI)

            # PP
            PP = meanSpikeReshaped[count, 2, i*2]
            tempNMI = (centerP + offsetP) / PP
            ppNMI.append(tempNMI)

            # NN
            NN = meanSpikeReshaped[count, 1, i*2-1]
            tempNMI = (centerN + offsetN) / NN
            nnNMI.append(tempNMI)

        ax4.axhline(y=1, linestyle='--', color='grey', alpha=0.6)
        ax4.axhline(y=2, linestyle='--', color='grey', alpha=0.6)
        ax4.plot(x, pCentNMI, color='green', label='PN')
        ax4.scatter(x, pCentNMI, color='green')
        ax4.plot(x, nCentNMI, color='red', label='NP')
        ax4.scatter(x, nCentNMI, color='red')
        ax4.plot(x, ppNMI, color='black', label='PP', linestyle='--')
        ax4.scatter(x, ppNMI, color='black')
        ax4.plot(x, nnNMI, color='grey', label='NN', linestyle='--')
        ax4.scatter(x, nnNMI, color='grey')
        ax4.set_xlabel('stimulus offset position')
        ax4.set_ylabel('NMI')
        ax4.set_ylim(bottom=0)
        ax4.set_title('NMI, Ruff and Cohen, 2016')
        ax4.legend()

        # NMI vs distance (Ni and Maunsell, 2017)
        ax5 = fig.add_subplot(gs06[0, 0])
        x = np.arange(1, numSteps + 1)
        pCentNMI = []
        pCentNMIBootstrapMean = []
        pCentNMIBootstrap2P = []
        pCentNMIBootstrap97P = []
        nCentNMI = []
        nCentNMIBootstrapMean = []
        nCentNMIBootstrap2P = []
        nCentNMIBootstrap97P = []
        ppNMI = []
        ppNMIBootstrapMean = []
        ppNMIBootstrap2P = []
        ppNMIBootstrap97P = []
        nnNMI = []
        nnNMIBootstrapMean = []
        nnNMIBootstrap2p = []
        nnNMIBootstrap97p = []

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

            if f'{seshDate}_{units[count]}' in unitList:
                totNMIPop.append(pnNMI)
                # tempWeight = meanSpikeReshaped[count, 0, i*2] / meanSpikeReshaped[count, 2, 0]
                # tempWeight = prefTransect[4+i] / np.max(prefTransect)
                # tempWeight = (gauss(i, *params)/np.max(respFull)) + (gauss(0, *params)/np.max(respFull))
                tempWeight = np.sqrt((gauss(i, *params) / np.max(respFull)) *
                                     (gauss(0, *params) / np.max(respFull)))
                totRFWeightPop.append(tempWeight)

            # bootstrap PN
            pIndex = stimCountIndex[2, 0]
            offIndex = stimCountIndex[0, i*2-1]
            PNIndex = stimCountIndex[2, i*2-1]
            pSpikes = spikeCountMat[count, :blocksDone, pIndex]
            bootstrapP = np.mean(np.random.choice(pSpikes, size=(10000, blocksDone),
                                                  replace=True), axis=1)
            offSpikes = spikeCountMat[count, :blocksDone, offIndex]
            bootstrapOff = np.mean(np.random.choice(offSpikes, size=(10000, blocksDone),
                                                    replace=True), axis=1)
            pnSpikes = spikeCountMat[count, :blocksDone, PNIndex]
            bootstrapPN = np.mean(np.random.choice(pnSpikes, size=(10000, blocksDone),
                                                   replace=True), axis=1)
            bootstrapNMI = ((bootstrapP + bootstrapOff) - bootstrapPN) / (
                            (bootstrapP + bootstrapOff) + bootstrapPN)
            pCentNMIBootstrapMean.append(np.mean(bootstrapNMI))
            pCentNMIBootstrap2P.append(np.percentile(bootstrapNMI, 2.5))
            pCentNMIBootstrap97P.append(np.percentile(bootstrapNMI, 97.5))

            # NP
            centerN = meanSpikeReshaped[count, 1, 0]
            offsetP = meanSpikeReshaped[count, 0, i*2]
            NP = meanSpikeReshaped[count, 1, i*2]
            npNMI = ((centerN + offsetP) - NP) / ((centerN + offsetP) + NP)
            nCentNMI.append(npNMI)

            if f'{seshDate}_{units[count]}' in unitList:
                totNMIPop.append(npNMI)
                # tempWeight = prefTransect[4+i] / np.max(prefTransect)
                tempWeight = np.sqrt((gauss(i, *params) / np.max(respFull)) *
                                     (gauss(0, *params) / np.max(respFull)))
                totRFWeightPop.append(tempWeight)

            # bootstrap NP
            cIndex = stimCountIndex[1, 0]
            offIndex = stimCountIndex[0, i*2]
            pairIndex = stimCountIndex[1, i*2]
            cSpikes = spikeCountMat[count, :blocksDone, cIndex]
            bootstrapC = np.mean(np.random.choice(cSpikes, size=(10000, blocksDone),
                                                  replace=True), axis=1)
            offSpikes = spikeCountMat[count, :blocksDone, offIndex]
            bootstrapOff = np.mean(np.random.choice(offSpikes, size=(10000, blocksDone),
                                                    replace=True), axis=1)
            pairSpikes = spikeCountMat[count, :blocksDone, pairIndex]
            bootstrapPair = np.mean(np.random.choice(pairSpikes, size=(10000, blocksDone),
                                                   replace=True), axis=1)
            bootstrapNMI = ((bootstrapC + bootstrapOff) - bootstrapPair) / (
                            (bootstrapC + bootstrapOff) + bootstrapPair)
            nCentNMIBootstrapMean.append(np.mean(bootstrapNMI))
            nCentNMIBootstrap2P.append(np.percentile(bootstrapNMI, 2.5))
            nCentNMIBootstrap97P.append(np.percentile(bootstrapNMI, 97.5))

            # PP
            PP = meanSpikeReshaped[count, 2, i*2]
            tempNMI = ((centerP + offsetP) - PP) / ((centerP + offsetP) + PP)
            ppNMI.append(tempNMI)

            if f'{seshDate}_{units[count]}' in unitList:
                totNMIPop.append(tempNMI)
                # tempWeight = prefTransect[4+i] / np.max(prefTransect)
                tempWeight = np.sqrt((gauss(i, *params) / np.max(respFull)) *
                                     (gauss(0, *params) / np.max(respFull)))
                totRFWeightPop.append(tempWeight)

            # NN
            NN = meanSpikeReshaped[count, 1, i*2-1]
            tempNMI = ((centerN + offsetN) - NN) / ((centerN + offsetN) + NN)
            nnNMI.append(tempNMI)

            if f'{seshDate}_{units[count]}' in unitList:
                totNMIPop.append(tempNMI)
                # tempWeight = prefTransect[4+i] / np.max(prefTransect)
                tempWeight = np.sqrt((gauss(i, *params) / np.max(respFull)) *
                                     (gauss(0, *params) / np.max(respFull)))
                totRFWeightPop.append(tempWeight)

        ax5.plot(x, pCentNMI, color='green', label='PN')
        ax5.scatter(x, pCentNMI, color='green')
        ax5.plot(x, nCentNMI, color='red', label='NP')
        ax5.scatter(x, nCentNMI, color='red')
        ax5.plot(x, ppNMI, color='black', linestyle='--', label='PP')
        ax5.scatter(x, ppNMI, color='black')
        ax5.plot(x, nnNMI, color='grey', linestyle='--', label='NN')
        ax5.scatter(x, nnNMI, color='grey')
        ax5.set_xlabel('stimulus offset position')
        ax5.set_ylabel('NMI')
        ax5.legend()
        ax5.set_title('NMI, Ni and Maunsell, 2017')

        ax6 = fig.add_subplot(gs07[0, 0])
        ax6.plot(x, pCentNMIBootstrapMean, color='green', label='PN')
        ax6.scatter(x, pCentNMIBootstrapMean, color='green')
        ax6.fill_between(x, pCentNMIBootstrap2P, pCentNMIBootstrap97P,
                         color='green', alpha=0.4)
        ax6.plot(x, nCentNMIBootstrapMean, color='red', label='NP')
        ax6.scatter(x, nCentNMIBootstrapMean, color='red')
        ax6.fill_between(x, nCentNMIBootstrap2P, nCentNMIBootstrap97P,
                         color='red', alpha=0.4)
        ax6.set_xlabel('stimulus offset position')
        ax6.set_ylabel('Bootstrapped NMI w. 95% CI')
        ax6.legend()

        # # predicted PN weighted
        # x = np.concatenate((np.arange(-numSteps, 0),
        #                     np.arange(0, numSteps + 1)),
        #                    axis=0)
        # transReverse = transectMeanSpike[count, ::-1]
        # prefTransect = np.concatenate((transReverse,
        #                                [meanSpikeReshaped[count, 2, 0]],
        #                                meanSpikeReshaped[count, 0, 2::2]),
        #                               axis=0)
        #
        # prefNonPrefRatio = meanSpikeReshaped[count, 1, 0] / meanSpikeReshaped[count, 2, 0]
        # params = gaussFit(x, prefTransect - spon)
        # nonPrefParams = [i if j != 1 else i * prefNonPrefRatio for j, i in enumerate(params)]
        # xFull = np.linspace(0, x[-1], 1000)
        # prefFull = gauss(xFull, *params)
        # nonPrefFull = gauss(xFull, *nonPrefParams)
        #
        # # pCentPred = (params[1] + nonPrefFull) / (1 + prefFull / params[1])
        # # npCentPred = (nonPrefParams[1] + prefFull) / (1 + prefFull / params[1])
        # pCentPred = (gauss(0, *params) + nonPrefFull) / (1 + prefFull / params[1])
        # npCentPred = (gauss(0, *nonPrefParams) + prefFull) / (1 + prefFull / params[1])
        #
        # # plt.plot(xFull, prefFull, color='black', linestyle='--')
        # # plt.plot(x, prefTransect, color='black')
        # # plt.plot(xFull, nonPrefFull, color='grey', linestyle='--')
        # ax2.plot(xFull, pCentPred, color='green', linestyle='--')
        # ax2.plot(xFull, npCentPred, color='red', linestyle='--')
        #
        # xFull = np.linspace(x[0], x[-1], 1000)
        # prefFull = gauss(xFull, *params)
        # ax3.plot(xFull, prefFull, color='black', linestyle='--')

        plt.tight_layout()
        plt.savefig(f'{seshDate}_{units[count]}.pdf')
        plt.close('all')

        if f'{seshDate}_{units[count]}' in unitList:
            normVal = np.max(meanSpike[count]) * 1000 / trueStimDurMS  # np.max(prefTransect)
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

            p0Indx = stimCountIndex[2, 0]
            adaptationMat[adaptC, 0, :blocksDone] = (spikeCountMat[count, :blocksDone, p0Indx] /
                                                     np.max(spikeCountMat[count, :blocksDone, p0Indx]))
            n0Indx = stimCountIndex[1, 0]
            adaptationMat[adaptC, 1, :blocksDone] = (spikeCountMat[count, :blocksDone, n0Indx] /
                                                     np.max(spikeCountMat[count, :blocksDone, n0Indx]))
            p1Indx = stimCountIndex[0, 2]
            adaptationMat[adaptC, 2, :blocksDone] = (spikeCountMat[count, :blocksDone, p1Indx] /
                                                     np.max(spikeCountMat[count, :blocksDone, p1Indx]))
            pT1Indx = transectCountIndex[0]
            adaptationMat[adaptC, 3, :blocksDone] = (spikeCountMat[count, :blocksDone, pT1Indx] /
                                                     np.max(spikeCountMat[count, :blocksDone, pT1Indx]))
            adaptC += 1

    # figure (hists)
    for count in range(len(units)):
        fig = plt.figure()
        fig.suptitle(f'unit {units[count]}')
        fig.set_size_inches(18, 10)
        gs0 = gridspec.GridSpec(4, 1)
        gs00 = gridspec.GridSpecFromSubplotSpec(1, numSteps+1, subplot_spec=gs0[0, 0])
        gs01 = gridspec.GridSpecFromSubplotSpec(1, numSteps+1, subplot_spec=gs0[1, 0])
        gs02 = gridspec.GridSpecFromSubplotSpec(1, numSteps+1, subplot_spec=gs0[2, 0])
        gs03 = gridspec.GridSpecFromSubplotSpec(1, numSteps+1, subplot_spec=gs0[3, 0])

        # initialize some variables
        yMax = 0
        smoothValue = 5
        spon = meanSpikeReshaped[count, 0, 0]
        axList = []

        # P center, N center
        ax = fig.add_subplot(gs00[0, 0])
        axList.append(ax)
        pIndx = stimCountIndex[2, 0]
        nIndx = stimCountIndex[1, 0]
        pHist = spikeHists[count, pIndx] * 1000 / blocksDone
        nHist = spikeHists[count, nIndx] * 1000 / blocksDone
        gaussSmoothP = gaussian_filter1d(pHist, smoothValue)
        gaussSmoothN = gaussian_filter1d(nHist, smoothValue)
        if np.max([gaussSmoothP, gaussSmoothN]) > yMax:
            yMax = np.max([gaussSmoothP, gaussSmoothN])
        ax.plot(gaussSmoothP, label='P0', color='black')
        ax.plot(gaussSmoothN, label='N0', color='grey')
        ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                       2 * histPrePostMS + trueStimDurMS])
        ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                            trueStimDurMS + histPrePostMS], fontsize=5)
        ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                   color='grey', alpha=0.2)
        ax.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
        ax.set_xlabel('Time (ms)', fontsize=7)
        ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
        ax.set_title('P Center N Center')
        ax.legend(fontsize=5)

        # P Center, N Offsets, and PN
        for i in range(1, numSteps+1):
            ax = fig.add_subplot(gs00[0, i])
            axList.append(ax)
            nOffIndx = stimCountIndex[0, i*2-1]
            pnIndx = stimCountIndex[2, i*2-1]
            nOffHist = spikeHists[count, nOffIndx] * 1000 / blocksDone
            pnHist = spikeHists[count, pnIndx] * 1000 / blocksDone
            gaussSmoothNOff = gaussian_filter1d(nOffHist, smoothValue)
            gaussSmoothPN = gaussian_filter1d(pnHist, smoothValue)
            if np.max([gaussSmoothNOff, gaussSmoothPN]) > yMax:
                yMax = np.max([gaussSmoothNOff, gaussSmoothPN])
            ax.plot(gaussSmoothP, label='P0', color='black')
            ax.plot(gaussSmoothNOff, label='N Offset', color='grey')
            ax.plot(gaussSmoothPN, label='PN', color='green')
            ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                           2 * histPrePostMS + trueStimDurMS])
            ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                                trueStimDurMS + histPrePostMS], fontsize=5)
            ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                       color='grey', alpha=0.2)
            ax.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
            ax.set_xlabel('Time (ms)', fontsize=7)
            ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
            ax.set_title(f'P Center, N Offset {i}')
            ax.legend(fontsize=5)

        # P center, N center
        ax = fig.add_subplot(gs01[0, 0])
        axList.append(ax)
        pIndx = stimCountIndex[2, 0]
        nIndx = stimCountIndex[1, 0]
        pHist = spikeHists[count, pIndx] * 1000 / blocksDone
        nHist = spikeHists[count, nIndx] * 1000 / blocksDone
        gaussSmoothP = gaussian_filter1d(pHist, smoothValue)
        gaussSmoothN = gaussian_filter1d(nHist, smoothValue)
        if np.max([gaussSmoothP, gaussSmoothN]) > yMax:
            yMax = np.max([gaussSmoothP, gaussSmoothN])
        ax.plot(gaussSmoothP, label='P0', color='black')
        ax.plot(gaussSmoothN, label='N0', color='grey')
        ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                       2 * histPrePostMS + trueStimDurMS])
        ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                            trueStimDurMS + histPrePostMS], fontsize=5)
        ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                   color='grey', alpha=0.2)
        ax.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
        ax.set_xlabel('Time (ms)', fontsize=7)
        ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
        ax.set_title('P Center N Center')
        ax.legend(fontsize=5)

        # N Center, P Offset, and NP
        for i in range(1, numSteps+1):
            ax = fig.add_subplot(gs01[0, i])
            axList.append(ax)
            pOffIndx = stimCountIndex[0, i*2]
            npIndx = stimCountIndex[1, i*2]
            pOffHist = spikeHists[count, pOffIndx] * 1000 / blocksDone
            npHist = spikeHists[count, npIndx] * 1000 / blocksDone
            gaussSmoothPOff = gaussian_filter1d(pOffHist, smoothValue)
            gaussSmoothNP = gaussian_filter1d(npHist, smoothValue)
            if np.max([gaussSmoothPOff, gaussSmoothNP]) > yMax:
                yMax = np.max([gaussSmoothPOff, gaussSmoothNP])
            ax.plot(gaussSmoothN, label='N0', color='grey')
            ax.plot(gaussSmoothPOff, label='P Offset', color='black')
            ax.plot(gaussSmoothNP, label='PN', color='red')
            ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                           2 * histPrePostMS + trueStimDurMS])
            ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                                trueStimDurMS + histPrePostMS], fontsize=5)
            ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                       color='grey', alpha=0.2)
            ax.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
            ax.set_xlabel('Time (ms)', fontsize=7)
            ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
            ax.set_title(f'N Center, P Offset {i}')
            ax.legend(fontsize=5)

        # P center, N center
        ax = fig.add_subplot(gs02[0, 0])
        axList.append(ax)
        pIndx = stimCountIndex[2, 0]
        nIndx = stimCountIndex[1, 0]
        pHist = spikeHists[count, pIndx] * 1000 / blocksDone
        nHist = spikeHists[count, nIndx] * 1000 / blocksDone
        gaussSmoothP = gaussian_filter1d(pHist, smoothValue)
        gaussSmoothN = gaussian_filter1d(nHist, smoothValue)
        if np.max([gaussSmoothP, gaussSmoothN]) > yMax:
            yMax = np.max([gaussSmoothP, gaussSmoothN])
        ax.plot(gaussSmoothP, label='P0', color='black')
        ax.plot(gaussSmoothN, label='N0', color='grey')
        ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                       2 * histPrePostMS + trueStimDurMS])
        ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                            trueStimDurMS + histPrePostMS], fontsize=5)
        ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                   color='grey', alpha=0.2)
        ax.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
        ax.set_xlabel('Time (ms)', fontsize=7)
        ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
        ax.set_title('P Center N Center')
        ax.legend(fontsize=5)

        # P Center, P Offsets, and PP
        for i in range(1, numSteps+1):
            ax = fig.add_subplot(gs02[0, i])
            axList.append(ax)
            pOffIndx = stimCountIndex[0, i*2]
            ppIndx = stimCountIndex[2, i*2]
            pOffHist = spikeHists[count, pOffIndx] * 1000 / blocksDone
            ppHist = spikeHists[count, ppIndx] * 1000 / blocksDone
            gaussSmoothPOff = gaussian_filter1d(pOffHist, smoothValue)
            gaussSmoothPP = gaussian_filter1d(ppHist, smoothValue)
            if np.max([gaussSmoothPOff, gaussSmoothPP]) > yMax:
                yMax = np.max([gaussSmoothPOff, gaussSmoothPP])
            ax.plot(gaussSmoothP, label='P0', color='black')
            ax.plot(gaussSmoothPOff, label='P Offset', color='black', alpha=0.5)
            ax.plot(gaussSmoothPP, label='PP', color='black', linestyle='--')
            ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                           2 * histPrePostMS + trueStimDurMS])
            ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                                trueStimDurMS + histPrePostMS], fontsize=5)
            ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                       color='grey', alpha=0.2)
            ax.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
            ax.set_xlabel('Time (ms)', fontsize=7)
            ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
            ax.set_title(f'P Center, P Offset {i}')
            ax.legend(fontsize=5)

        # P center, N center
        ax = fig.add_subplot(gs03[0, 0])
        axList.append(ax)
        pIndx = stimCountIndex[2, 0]
        nIndx = stimCountIndex[1, 0]
        pHist = spikeHists[count, pIndx] * 1000 / blocksDone
        nHist = spikeHists[count, nIndx] * 1000 / blocksDone
        gaussSmoothP = gaussian_filter1d(pHist, smoothValue)
        gaussSmoothN = gaussian_filter1d(nHist, smoothValue)
        if np.max([gaussSmoothP, gaussSmoothN]) > yMax:
            yMax = np.max([gaussSmoothP, gaussSmoothN])
        ax.plot(gaussSmoothP, label='P0', color='black')
        ax.plot(gaussSmoothN, label='N0', color='grey')
        ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                       2 * histPrePostMS + trueStimDurMS])
        ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                            trueStimDurMS + histPrePostMS], fontsize=5)
        ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                   color='grey', alpha=0.2)
        ax.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
        ax.set_xlabel('Time (ms)', fontsize=7)
        ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
        ax.set_title('P Center N Center')
        ax.legend(fontsize=5)

        # N Center, N Offsets, and NN
        yMaxNull = 0
        for i in range(1, numSteps+1):
            ax = fig.add_subplot(gs03[0, i])
            axList.append(ax)
            nOffIndx = stimCountIndex[0, i*2-1]
            nnIndx = stimCountIndex[1, i*2-1]
            nOffHist = spikeHists[count, nOffIndx] * 1000 / blocksDone
            nnHist = spikeHists[count, nnIndx] * 1000 / blocksDone
            gaussSmoothNOff = gaussian_filter1d(nOffHist, smoothValue)
            gaussSmoothNN = gaussian_filter1d(nnHist, smoothValue)
            if np.max([gaussSmoothNOff, gaussSmoothNN]) > yMax:
                yMax = np.max([gaussSmoothNOff, gaussSmoothNN])
            if np.max([gaussSmoothNOff, gaussSmoothNN]) > yMaxNull:
                yMaxNull = np.max([gaussSmoothNOff, gaussSmoothNN])
            ax.plot(gaussSmoothN, label='N0', color='grey')
            ax.plot(gaussSmoothNOff, label='N Offset', color='grey', alpha=0.5)
            ax.plot(gaussSmoothNN, label='NN', color='grey', linestyle='--')
            ax.set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                           2 * histPrePostMS + trueStimDurMS])
            ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                                trueStimDurMS + histPrePostMS], fontsize=5)
            ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                       color='grey', alpha=0.2)
            ax.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)
            ax.set_xlabel('Time (ms)', fontsize=7)
            ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
            ax.set_title(f'N Center, N Offset {i}')
            ax.legend(fontsize=5)

        for ax in axList:
            ax.set_ylim([0, yMax * 1.1])

        for ax in axList[-4:]:
            ax.set_ylim([0, yMaxNull * 1.1])

        plt.tight_layout()
        plt.savefig(f'{seshDate}_{units[count]} histogram.pdf')
        plt.close('all')

    # use this if I need to run another data file immediately afterwards
    os.chdir('../../../Python Code')

print(time.time()-p0)


# condition resp normalized by max resp vs correlation
fig, ax = plt.subplots()
fig.suptitle('normalized correlations vs condition response normalized by max response')
a = np.array(normalizedRespWeightPop)
b = np.array(totCorrPop)
indx = np.where((~np.isnan(a)) & (~np.isnan(b)))
c = stats.linregress(a[indx],b[indx], alternative='two-sided')
ax.scatter(normalizedRespWeightPop, totCorrPop)
ax.set_xlabel('pair normalized resp')
ax.set_ylabel('pair correlation')
ax.set_title(f'corr: {c[2]:.2f}, slope = {c[0]:.2f}, pval of slope= {c[3]:.5f}')
plt.show()

# pair NMI with offset position
pairNMIPop = np.array(pairNMIPop)
pairedCorrPop = np.array(pairedCorrPop)
a = pairNMIPop.reshape((len(sponCorrPop), int(len(pairNMIPop)/len(sponCorrPop))))
a1 = a[:, :4].flatten()
a2 = a[:, 4:8].flatten()
a3 = a[:, 8:12].flatten()
a4 = a[:, 12:16].flatten()
# a1 = a[:, :2].flatten()
# a2 = a[:, 2:4].flatten()
# a3 = a[:, 4:6].flatten()
# a4 = a[:, 6:8].flatten()

a1Mean = np.nanmean(a1)
a1SEM = np.nanstd(a1) / np.sqrt(len(a1))
a2Mean = np.nanmean(a2)
a2SEM = np.nanstd(a2) / np.sqrt(len(a2))
a3Mean = np.nanmean(a3)
a3SEM = np.nanstd(a3) / np.sqrt(len(a3))
a4Mean = np.nanmean(a4)
a4SEM = np.nanstd(a4) / np.sqrt(len(a4))

# figure
fig, (ax, ax1, ax2) = plt.subplots(1, 3)
fig.set_size_inches(14, 4)
fig.suptitle('NMI vs offset position and correlation, for paired conditions (PN, NP), neurons (P+NP)')

ax.plot([1, 2, 3, 4], [a1Mean, a2Mean, a3Mean, a4Mean])
ax.errorbar([1, 2, 3, 4], [a1Mean, a2Mean, a3Mean, a4Mean],
            yerr=[a1SEM, a2SEM, a3SEM, a4SEM], ecolor='black', color='black',
            markersize=4)
ax.set_xlabel('offset position')
ax.set_ylabel('pair NMI')
ax.set_title('pair NMI vs offset position')

# correlations vs pair NMI
indx = np.where((~np.isnan(pairNMIPop)) & (~np.isnan(pairedCorrPop)))[0]
res = stats.linregress(pairNMIPop[indx], pairedCorrPop[indx],
                       alternative='two-sided')
ax1.scatter(pairNMIPop, pairedCorrPop)
ax1.set_xlabel('pair average NMI')
ax1.set_ylabel('pair correlation')
ax1.set_title(f'pair NMI vs Rsc, corr = {res[2]:.2f}, slope = {res[0]:.2f}, '
              f'pval = {res[3]:.4f}')

# correlations changes with offset
off1PopMean = np.nanmean(off1CorrPop)
off1SEM = np.nanstd(off1CorrPop) / np.sqrt(len(off1CorrPop))
off2PopMean = np.nanmean(off2CorrPop)
off2SEM = np.nanstd(off2CorrPop) / np.sqrt(len(off2CorrPop))
off3PopMean = np.nanmean(off3CorrPop)
off3SEM = np.nanstd(off3CorrPop) / np.sqrt(len(off3CorrPop))
off4PopMean = np.nanmean(off4CorrPop)
off4SEM = np.nanstd(off4CorrPop) / np.sqrt(len(off4CorrPop))

corrArrMean = [off1PopMean, off2PopMean, off3PopMean, off4PopMean]
corrArr = [off1CorrPop, off2CorrPop, off3CorrPop, off4CorrPop]

SEMArr = [off1SEM, off2SEM, off3SEM, off4SEM]
x = np.arange(1, 5, 1)

ax2.scatter(x, corrArrMean)
ax2.plot(x, corrArrMean)
# ax.violinplot(corrArr, x, widths=0.5, showmeans=True, showextrema=True, showmedians=True,
#               bw_method=0.5)
ax2.errorbar(x, corrArrMean, yerr=SEMArr, fmt='o', ecolor='black', color='black',
            markersize=2)
ax2.set_ylabel('Spike Count Correlation (pop mean)')
ax2.set_xlabel('Offset Position')
ax2.set_title('Offset Position vs Rsc blocks (all blocks)')

plt.tight_layout()
plt.show()


# rf weight vs correlations (binned)
n = int(len(totCorrPop) / 4)
rfWeightPop = np.array(rfWeightPop)
rfWeightPairedPop = np.array(rfWeightPairedPop)
rfWeightSinglePop = np.array(rfWeightSinglePop)
totCorrPop = np.array(totCorrPop)
pairedCorrPop = np.array(pairedCorrPop)
singleCorrPop = np.array(singleCorrPop)
sponCorrPop = np.array(sponCorrPop)
sortIndex = np.argsort(rfWeightPop)
sortedRFWeight = rfWeightPop[sortIndex]
sortedPopCorr = totCorrPop[sortIndex]

# manual bins - equally populated bins (total pop)
equalBinsRFWeight = [sortedRFWeight[i:i + n] for i in range(0, len(sortedRFWeight), n)]
equalBinsCorr = [sortedPopCorr[i:i + n] for i in range(0, len(sortedPopCorr), n)]
binMeanRFWeight = np.array([np.nanmean(i) for i in equalBinsRFWeight])
binSEMRFWeight = np.array([np.nanstd(i) for i in equalBinsRFWeight]) / np.sqrt(n)
binMeanCorr = np.array([np.nanmean(i) for i in equalBinsCorr])
binSEMCorr = np.array([np.nanstd(i) for i in equalBinsCorr]) / np.sqrt(n)

# manual bins - equally populated bins (paired pop)
n = int(len(pairedCorrPop) / 4)
sortIndex = np.argsort(rfWeightPairedPop)
sortedRFWeightPaired = rfWeightPairedPop[sortIndex]
sortedPopCorrPaired = pairedCorrPop[sortIndex]
binsRFWeightPaired = [sortedRFWeightPaired[i:i + n] for i in range(0, len(sortedRFWeightPaired), n)]
binsCorrPaired = [sortedPopCorrPaired[i:i + n] for i in range(0, len(sortedPopCorrPaired), n)]
binMeanRFWeightPaired = np.array([np.nanmean(i) for i in binsRFWeightPaired])
binMeanCorrPaired = np.array([np.nanmean(i) for i in binsCorrPaired])
binSEMCorrPaired = np.array([np.nanstd(i) for i in binsCorrPaired]) / np.sqrt(n)

# manual bins - equally populated bins (single pop)
n = int(len(singleCorrPop) / 4)
sortIndex = np.argsort(rfWeightSinglePop)
sortedRFWeightSingle = rfWeightSinglePop[sortIndex]
sortedPopCorrSingle = singleCorrPop[sortIndex]
binsRFWeightSingle = [sortedRFWeightSingle[i:i + n] for i in range(0, len(sortedRFWeightSingle), n)]
binsCorrSingle = [sortedPopCorrSingle[i:i + n] for i in range(0, len(sortedPopCorrSingle), n)]
binMeanRFWeightSingle = np.array([np.nanmean(i) for i in binsRFWeightSingle])
binMeanCorrSingle = np.array([np.nanmean(i) for i in binsCorrSingle])
binSEMCorrSingle = np.array([np.nanstd(i) for i in binsCorrSingle]) / np.sqrt(n)

# ANOVA
m = np.where((~np.isnan(equalBinsCorr[0])) & (~np.isnan(equalBinsCorr[1]))
             & (~np.isnan(equalBinsCorr[2])) & (~np.isnan(equalBinsCorr[3])))[0]
fTest = f_oneway(equalBinsCorr[0][m], equalBinsCorr[1][m],
                 equalBinsCorr[2][m], equalBinsCorr[3][m])
fStat, fPValue = fTest[0], fTest[1]

# # paired corr ANOVA
# m = np.where((~np.isnan(binsCorrPaired[0])) & (~np.isnan(binsCorrPaired[1]))
#              & (~np.isnan(binsCorrPaired[2])) & (~np.isnan(binsCorrPaired[3])))[0]
# fTest = f_oneway(binsCorrPaired[0][m], binsCorrPaired[1][m],
#                  binsCorrPaired[2][m], binsCorrPaired[3][m])
# fStat, fPValue = fTest[0], fTest[1]
#
# # single corr ANOVA
# m = np.where((~np.isnan(binsCorrSingle[0])) & (~np.isnan(binsCorrSingle[1]))
#              & (~np.isnan(binsCorrSingle[2])) & (~np.isnan(binsCorrSingle[3])))[0]
# fTest = f_oneway(binsCorrSingle[0][m], binsCorrSingle[1][m],
#                  binsCorrSingle[2][m], binsCorrSingle[3][m])
# fStat, fPValue = fTest[0], fTest[1]

# figure
fig, (ax, ax1, ax2) = plt.subplots(1, 3)
fig.set_size_inches(15, 4)
fig.suptitle('correlations vs rf weight for all stim (single+paired) and neurons (P+NP)')
ax.plot(binMeanRFWeight, binMeanCorr, color='black')
ax.errorbar(binMeanRFWeight, binMeanCorr, yerr=binSEMCorr, fmt='o', ecolor='black',
            color='black', markersize=2, label='total stimuli')
ax.plot(binMeanRFWeightPaired, binMeanCorrPaired, color='red')
ax.errorbar(binMeanRFWeightPaired, binMeanCorrPaired, yerr=binSEMCorrPaired,
            fmt='o', ecolor='red', color='red', markersize=2, label='paired stimuli')
ax.plot(binMeanRFWeightSingle, binMeanCorrSingle, color='green')
ax.errorbar(binMeanRFWeightSingle, binMeanCorrSingle, yerr=binSEMCorrSingle,
            fmt='o', ecolor='green', color='green', markersize=2, label='Single stimuli')
ax.set_ylabel('Spike Count Correlation (bin mean)')
ax.set_xlabel('binned RF weights, pair average')
ax.axhline(y=np.nanmean(sponCorrPop), linestyle='--', color='blue', label='spon')
# polyfit
lineParams = np.polyfit(binMeanRFWeight, binMeanCorr, 1)
trendLine = np.poly1d(lineParams)
xLine = np.linspace(min(binMeanRFWeight), max(binMeanRFWeight), 100)
ax.plot(xLine, trendLine(xLine), color='black', linestyle='--', label="fitted line")
if fPValue < 0.05:
    ax.set_title(f'ANOVA p-value < 0.05, pval={fPValue:.5f}')
else:
    ax.set_title('ANOVA p-value > 0.05')
ax.invert_xaxis()
ax.legend()

# scatter of RF weight vs correlation
indx = np.where((~np.isnan(rfWeightPop)) & (~np.isnan(totCorrPop)))[0]
linRegressResults = stats.linregress(rfWeightPop[indx],
                                     totCorrPop[indx],
                                     alternative='two-sided')
ax1.scatter(rfWeightPop, totCorrPop)
ax1.set_xlabel('rf weight, pair average')
ax1.set_ylabel('pair spike count correlation')
nonNanIndex = np.where((~np.isnan(totCorrPop)) & (~np.isnan(rfWeightPop)))[0]
corr = stats.pearsonr(rfWeightPop[nonNanIndex], totCorrPop[nonNanIndex])[0]
# polyfit
indx = np.where((~np.isnan(rfWeightPop)) & (~np.isnan(totCorrPop)))[0]
scatterLineParams = np.polyfit(rfWeightPop[indx], totCorrPop[indx], 1)
trendLine = np.poly1d(scatterLineParams)
xLine = np.linspace(min(rfWeightPop), max(rfWeightPop), 100)
ax1.plot(xLine, trendLine(xLine), color='black', linestyle='--', label="fitted line")
ax1.set_title(f'pair rf weight vs rsc, corr = {corr:.2f}, slope = {scatterLineParams[0]:.2f}, '
              f'p_value = {linRegressResults[3]:.4f}')
ax1.invert_xaxis()

# shuffle arrays and compare slope of shuffled arrays to actual slope
shuffSlopes = []
for i in np.arange(5000):
    a = np.copy(equalBinsCorr)
    np.random.shuffle(a.flat)
    binMeans = np.nanmean(a, axis=1)
    shuffParams = np.polyfit(binMeanRFWeight, binMeans, 1)
    shuffSlopes.append(shuffParams[0])
ax2.hist(shuffSlopes, bins=20)
ax2.axvline(lineParams[0], linestyle='--', label='Real Slope')
ax2.set_xlabel('shuffled slopes')
ax2.set_ylabel('count')
ax2.set_title(f'real slope: {lineParams[0]:.2f}, mean shuffled slope: {np.mean(shuffSlopes):.2f}')

plt.tight_layout()
plt.show()

# plot of rf weight with NMI
n = int(len(totRFWeightPop) / 4)
totRFWeightPop = np.array(totRFWeightPop)
totNMIPop = np.array(totNMIPop)
sortIndex = np.argsort(totRFWeightPop)
sortedRFWeight = totRFWeightPop[sortIndex]
sortedNMIPop = totNMIPop[sortIndex]

# manual bins - equally populated bins
equalBinsRFWeight = [sortedRFWeight[i:i + n] for i in range(0, len(sortedRFWeight), n)]
equalBinsNMIPop = [sortedNMIPop[i:i + n] for i in range(0, len(sortedNMIPop), n)]
binMeanRFWeight = np.array([np.nanmean(i) for i in equalBinsRFWeight])
binMeanNMIPop = np.array([np.nanmean(i) for i in equalBinsNMIPop])
binSEMNMIPop = np.array([np.std(i) for i in equalBinsNMIPop]) / np.sqrt(n)

# figure
fig, (ax, ax1) = plt.subplots(1, 2)
fig.suptitle('individual neurons pref+nonpref sessions')
fig.set_size_inches(10, 4)
ax.plot(binMeanRFWeight, binMeanNMIPop, color='black')
ax.errorbar(binMeanRFWeight, binMeanNMIPop, yerr=binSEMNMIPop, fmt='o', ecolor='black',
            color='black', markersize=2)
ax.set_ylabel('NMI (bin median)')
ax.set_xlabel('binned RF weights')
ax.set_title('binned RF weight vs NMI')
ax.invert_xaxis()

indx = np.where((~np.isnan(totRFWeightPop)) & (~np.isnan(totNMIPop)))
c = stats.linregress(totRFWeightPop[indx], totNMIPop[indx], alternative='two-sided')
ax1.scatter(totRFWeightPop, totNMIPop)
ax1.set_xlabel('RF weights')
ax1.set_ylabel('NMI')
ax1.set_title(f'RF weight vs NMI, corr = {c[2]:.2f}, slope = {c[0]:.2f}, '
              f'slope pval= {c[3]:.5f}')
ax1.invert_xaxis()

plt.show()


# # Population plots
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

# figure
for filler in range(1):
    fig, ((ax2, ax3, ax1, ax4), (ax7, ax8, ax5, ax6)) = plt.subplots(2, 4)
    fig.set_size_inches(16, 8)

    # bin size = n, exponent for prediction = exp
    n = numUnits
    exp = 1
    fig.suptitle(f'population responses pref+nonpref neurons, exp={exp}')
    hfont = {'fontname': 'Arial'}

    # transect response
    transX = np.concatenate((np.arange(-numSteps, 0),
                             np.arange(0, numSteps + 1)),
                            axis=0)
    transectMean = np.mean(transectNormalized, axis=0)
    transectSEM = np.std(transectNormalized, axis=0) / np.sqrt(numUnits)
    params = gaussFit(transX, transectMean)
    xFull = np.linspace(transX[0], transX[-1], 1000)
    respFull = gauss(xFull, *params)
    ax1.plot(xFull, respFull, color='black', linestyle='--', label='Gauss Fit')
    ax1.plot(transX, transectMean, color='black', label='Pref')
    ax1.errorbar(transX, transectMean, yerr=transectSEM, fmt='o', ecolor='black',
                 color='black', markersize=2)
    for i in transectNormalized:
        ax1.plot(transX, i, color='purple', alpha=0.2)
    ax1.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
    ax1.fill_between(x=transX, y1=meanSpon+semSpon, y2=meanSpon-semSpon,
                     color='blue', alpha=0.2)
    ax1.set_ylabel('Normalized Response', fontsize=15, **hfont)

    ax1.set_xlabel('Stimulus Offset Positions', fontsize=15, **hfont)
    ax1.set_ylim([0, 1])
    ax1.set_title('transect')

    # normalization with distance PN/NP
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
    # params = gaussFit(x, prefMean)
    # xFull = np.linspace(x[0], x[-1], 1000)
    # respFull = gauss(xFull, *params)

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


    # RF weight is applied to numerator as well (USE THIS)
    # pCentPred = ((prefMean[0] / np.max(respFull)) + (nonprefMean[1:] * (prefMean[1:]/np.max(respFull)))) / (
    #             (prefMean[0] + prefMean[1:]) / np.max(respFull))
    #
    # npCentPred = ((nonprefMean[0] / np.max(respFull)) + (prefMean[1:] * (prefMean[1:]/np.max(respFull)))) / (
    #              (prefMean[0] + prefMean[1:]) / np.max(respFull))

    # RF weight is by max of fitted Gaussian (OLD WAY)
    pCentPred = (prefMean[0] + (nonprefMean[1:]) / np.max(respFull)) / (
                (prefMean[0] + prefMean[1:]) / np.max(respFull))
    npCentPred = (nonprefMean[0] + (prefMean[1:]) / np.max(respFull)) / (
                 (prefMean[0] + prefMean[1:]) / np.max(respFull))

    # # simple average
    # pCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / 2
    # npCentPred = ((nonprefMean[0] + prefMean[1:])**exp) / 2

    ax2.plot(x, prefMean, color='black', label='PO')
    ax2.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
                 color='black', markersize=2)
    ax2.plot(x, nonprefMean, color='grey', label='NO')
    ax2.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=2)
    ax2.plot(xNorm, pnMean, color='green', label='P0 N1')
    ax2.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
                 color='green', markersize=2)
    ax2.plot(xNorm, npMean, color='red', label='N0 P1')
    ax2.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='red',
                 color='red', markersize=2)
    ax2.plot(xNorm, pCentPred, color='green', linestyle='--')
    ax2.plot(xNorm, npCentPred, color='red', linestyle='--')
    ax2.set_xlabel('stimulus offset positions')
    ax2.set_ylabel('normalized response')
    ax2.set_ylim([0, 1])
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
    ppCentPred = ((prefMean[0] + prefMean[1:])**exp) / (
                  (prefMean[0] + prefMean[1:]) / np.max(transectMean))
    npnpCentPred = ((nonprefMean[0] + nonprefMean[1:])**exp) / (
                    (prefMean[0] + prefMean[1:]) / np.max(transectMean))
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
    ax3.set_ylim([0, 1])
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

    # raw vs predicted norm PN, NP, PP, NN
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

    # # simple average prediction
    # pCentPred = ((np.reshape(prefNormalized[:, 0], (numUnits, 1)) + nonprefNormalized[:, 1:])**exp) / 2
    # meanPCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / 2
    # nCentPred = ((np.reshape(nonprefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:])**exp) / 2
    # meanNPCentPred = ((nonprefMean[0] + prefMean[1:])**exp) / 2

    pCentPredReshape = np.reshape(pCentPred, numUnits*4)
    pnNormalizedReshape = np.reshape(pnNormalized, numUnits*4)
    nCentPredReshape = np.reshape(nCentPred, numUnits*4)
    npNormalizedReshape = np.reshape(npNormalized, numUnits*4)
    ppPredReshape = np.reshape(ppPred, numUnits*4)
    ppNormalizedReshape = np.reshape(ppNormalized, numUnits*4)
    nnPredReshape = np.reshape(nnPred, numUnits*4)
    nnNormalizedReshape = np.reshape(nnNormalized, numUnits*4)
    ax5.scatter(pnNormalizedReshape, pCentPredReshape, color='green',
                label='p center', alpha=0.4)
    ax5.scatter(npNormalizedReshape, nCentPredReshape, color='red',
                label='n center', alpha=0.4)
    # ax5.scatter(ppNormalizedReshape, ppPredReshape, color='black',
    #             label='pp', alpha=1)
    # ax5.scatter(nnNormalizedReshape, nnPredReshape, color='grey',
    #             label='nn', alpha=0.4)
    ax5.scatter(pnMean, meanPCentPred, color='green', edgecolor='gold',
                s=80)
    ax5.scatter(npMean, meanNPCentPred, color='red', edgecolor='gold',
                s=80)
    # ax5.scatter(ppMean, meanPPCentPred, color='blue', edgecolor='gold',
    #             s=80)
    # ax5.scatter(nnMean, meanNNCentPred, color='blue', edgecolor='gold',
    #             s=80)
    ax5.set_ylim([0, 1.5])
    ax5.set_xlim([0, 1.5])
    ax5.set_xlabel('Real Response', fontsize=15, **hfont)
    ax5.set_ylabel('Predicted Response', fontsize=15, **hfont)
    ax5.set_title('real vs predicated normalization response')
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

    # plt.plot(binMeanTransectSep, binMeanTransectResp)
    # plt.plot(transXFull, transFitResp)
    # plt.show()

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
    pnPred = (binMeanPrefResp[0] + (binMeanNonprefResp[1:] * (binMeanPrefResp[1:] / np.max(transFitResp)))) / (
        (binMeanPrefResp[0] + binMeanPrefResp[1:]) / np.max(transFitResp))
    npPred = (binMeanNonprefResp[0] + (binMeanPrefResp[1:] * (binMeanPrefResp[1:] / np.max(transFitResp)))) / (
        (binMeanPrefResp[0] + binMeanPrefResp[1:]) / np.max(transFitResp))

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


    # # rf weight vs correlations (binned)
    # n = 68
    # rfWeightPop = np.array(rfWeightPop)
    # totCorrPop = np.array(totCorrPop)
    # sortIndex = np.argsort(rfWeightPop)
    # sortedRFWeight = rfWeightPop[sortIndex]
    # sortedPopCorr = totCorrPop[sortIndex]
    # # manual bins - equally populated bins
    # equalBinsRFWeight = [sortedRFWeight[i:i + n] for i in range(0, len(sortedRFWeight), n)]
    # equalBinsCorr = [sortedPopCorr[i:i + n] for i in range(0, len(sortedPopCorr), n)]
    # binMeanRFWeight = np.array([np.mean(i) for i in equalBinsRFWeight])
    # binSEMRFWeight = np.array([np.std(i) for i in equalBinsRFWeight]) / np.sqrt(n)
    # binMeanCorr = np.array([np.mean(i) for i in equalBinsCorr])
    # binSEMCorr = np.array([np.std(i) for i in equalBinsCorr]) / np.sqrt(n)
    #
    # ax8.plot(binMeanRFWeight, binMeanCorr, color='black')
    # ax8.errorbar(binMeanRFWeight, binMeanCorr, yerr=binSEMCorr, fmt='o', ecolor='black',
    #             color='black', markersize=2)
    # ax8.set_ylabel('Spike Count Correlation (bin mean)')
    # ax8.set_xlabel('binned RF weights, pair average')
    # ax8.set_ylim(bottom=0)
    # ax8.invert_xaxis()

    plt.tight_layout()
    plt.show()


# distribution of sigmas scaled by eccentricity for RF gauss, gabor, offset sep
# gauss fit for RF tuning
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
fig.set_size_inches(14, 4)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# distribution of RF gauss fit sigma / eccentricity
rfGaussSigma = []
for i in range(len(offsetDegSepNormPop)):
    params = gaussFit(offsetDegSepNormPop[i], transectNormalized[i])
    rfGaussSigma.append(params[3])
ax1.hist(rfGaussSigma, bins=6)
ax1.scatter(rfGaussSigma, np.random.uniform(0.5, 1.5, len(rfGaussSigma)),
            color='black')
ax1.set_title('RF Gaussian sigmas scaled w. Eccentricity', fontsize=8)
ax1.set_xlabel('RF Gaussian sigma / Eccentricity')
ax1.set_ylabel('count')
ax1.text(0.05, 0.95, f'median: {np.median(rfGaussSigma):.2f}',
         transform=ax1.transAxes, fontsize=10, verticalalignment='top',
         bbox=props)

# distribution of offset pos 1 deg separation / eccentricity
ax2.hist(offsetDegSepNormPop[:, 5], bins=6)
ax2.scatter(offsetDegSepNormPop[:, 5], np.random.uniform(0.5, 1.5, len(rfGaussSigma)),
            color='black')
ax2.set_title('Offset 1 separation scaled w. Eccentricity', fontsize=8)
ax2.set_xlabel('Offset 1 Deg Separation / Eccentricity')
ax2.set_ylabel('count')
ax2.text(0.05, 0.95, f'median: {np.median(offsetDegSepNormPop[:,5]):.2f}',
         transform=ax2.transAxes, fontsize=10, verticalalignment='top',
         bbox=props)

# distribution of offset pos 2 deg separation / eccentricity
ax3.hist(offsetDegSepNormPop[:, 6], bins=6)
ax3.scatter(offsetDegSepNormPop[:, 6], np.random.uniform(0.5, 1.5, len(rfGaussSigma)),
            color='black')
ax3.set_title('Offset 2 separation scaled w. Eccentricity', fontsize=8)
ax3.set_xlabel('Offset 2 Deg Separation / Eccentricity')
ax3.set_ylabel('count')
ax3.text(0.05, 0.95, f'median: {np.median(offsetDegSepNormPop[:,6]):.2f}',
         transform=ax3.transAxes, fontsize=10, verticalalignment='top',
         bbox=props)

# distribution of rfGabor sigmas/eccentricity
ax4.hist(rfGaborSigmaPop, bins=6)
ax4.scatter(rfGaborSigmaPop, np.random.uniform(0.5, 1.5, len(rfGaborSigmaPop)),
            color='black')
ax4.set_title('rfGabor Sigma scaled w. Eccentricity', fontsize=8)
ax4.set_xlabel('RF Gabor sigma / Eccentricity')
ax4.set_ylabel('count')
ax4.text(0.05, 0.95, f'median: {np.median(rfGaborSigmaPop):.2f}',
         transform=ax4.transAxes, fontsize=10, verticalalignment='top',
         bbox=props)

plt.tight_layout()
plt.show()

# adaptation mat
adaptationMean = np.nanmean(adaptationMat[:numUnits, :, :50], axis=0)
plt.plot(adaptationMean[0]); plt.show()




# retreat poster figures
fig, ax2 = plt.subplots()
hfont = {'fontname':'Arial'}
# normalization with distance PN/NP
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
params = gaussFit(x, prefMean)
xFull = np.linspace(x[0], x[-1], 1000)
respFull = gauss(xFull, *params)

# rf weighted average
pCentPred = ((prefMean[0] + nonprefMean[1:])**1) / (
            1 + prefMean[1:] / np.max(prefMean))
npCentPred = ((nonprefMean[0] + prefMean[1:])**1) / (
            1 + prefMean[1:] / np.max(prefMean))
pCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / (
            ((prefMean[0] + prefMean[1:])**1) / np.max(respFull))
npCentPred = ((nonprefMean[0] + prefMean[1:])**exp) / (
             ((prefMean[0] + prefMean[1:])**1) / np.max(respFull))

# # simple average
# pCentPred = ((prefMean[0] + nonprefMean[1:]) ** exp) / 2
# npCentPred = ((nonprefMean[0] + prefMean[1:]) ** exp) / 2

ax2.plot(x, prefMean, color='black', label='Preferred Stimulus')
ax2.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
             color='black', markersize=2)
ax2.plot(x, nonprefMean, color='grey', label='Nonpreferred Stimulus')
ax2.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
             color='grey', markersize=2)

# ax2.plot(xNorm, pnMean, color='green', label='Preferred Center, Non-preferred Offset')
# ax2.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
#              color='green', markersize=2)
# ax2.plot(xNorm, nonprefMean[1:], color='grey', label='Nonpreferred Stimulus Offset')
# ax2.errorbar(xNorm, nonprefMean[1:], yerr=nonprefSEM[1:], fmt='o', ecolor='grey',
#              color='grey', markersize=2)
# ax2.plot(xNorm, [prefMean[0]] * 4, color='black', label='Preferred Stimulus Center')
# ax2.errorbar(xNorm, [prefMean[0]] * 4, yerr=[prefSEM[0]] * 4, fmt='o', ecolor='black',
#              color='black', markersize=2)

# ax2.plot(xNorm, npMean, color='red', label='Nonpreferred Center, Preferred Offset')
# ax2.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='red',
#              color='red', markersize=2)
# ax2.plot(xNorm, [nonprefMean[0]] * 4, color='grey', label='Nonpreferred Stimulus Center')
# ax2.errorbar(xNorm, [nonprefMean[0]] * 4, yerr=[nonprefSEM[0]] * 4, fmt='o', ecolor='grey',
#              color='grey', markersize=2)
#
# ax2.plot(xNorm, prefMean[1:], color='black', label='Preferred Stimulus Offset')
# ax2.errorbar(xNorm, prefMean[1:], yerr=prefSEM[1:], fmt='o', ecolor='black',
#              color='black', markersize=2)

# ax2.plot(xNorm, pCentPred, color='green', linestyle='--')
# ax2.plot(xNorm, npCentPred, color='red', linestyle='--')
ax2.set_xlabel('Stimulus Offset Positions', fontsize=15, **hfont)
ax2.set_ylabel('Normalized Response', fontsize=15, **hfont)
ax2.set_ylim([0, 1])
ax2.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
ax2.fill_between(x=x, y1=meanSpon + semSpon, y2=meanSpon - semSpon,
                 color='blue', alpha=0.2)
ax2.set_xticks(np.arange(0, 5, step=1))
plt.show()


"""
# fit gaussian to pref transect
count = 7

x = np.concatenate((np.arange(-numSteps, 0),
                    np.arange(0, numSteps + 1)),
                   axis=0)
transReverse = transectMeanSpike[count, ::-1]
prefTransect = np.concatenate((transReverse,
                               [meanSpikeReshaped[count, 2, 0]],
                               meanSpikeReshaped[count, 0, 2::2]),
                              axis=0)

prefNonPrefRatio = meanSpikeReshaped[count, 1, 0] / meanSpikeReshaped[count, 2, 0]
params = gaussFit(x, prefTransect - spon)
nonPrefParams = [i if j != 1 else i*prefNonPrefRatio for j, i in enumerate(params)]
xFull = np.linspace(x[0], x[-1], 1000)
prefFull = gauss(xFull, *params)
nonPrefFull = gauss(xFull, *nonPrefParams)

# pCentPred = (params[1] + nonPrefFull) / (1 + prefFull / params[1])
# npCentPred = (nonPrefParams[1] + prefFull) / (1 + prefFull / params[1])

pCentPred = (gauss(0, *params) + nonPrefFull) / (1 + prefFull / params[1])
npCentPred = (gauss(0, *nonPrefParams) + prefFull) / (1 + prefFull / params[1])

# pCentPred = (meanSpikeReshaped[count, 2, 0] + nonPrefFull) / (1 + prefFull / params[1])
# npCentPred = (meanSpikeReshaped[count, 1, 0] + prefFull) / (1 + prefFull / params[1])

plt.plot(x, prefTransect, color='black')
plt.plot(xFull, prefFull, color='black', linestyle='--')
plt.plot(xFull, nonPrefFull, color='grey', linestyle='--')
plt.plot(xFull, pCentPred, color='green', linestyle='--')
plt.plot(xFull, npCentPred, color='red', linestyle='--')

plt.show()

# plot individual gaussian fits to RF using offset positions converted to separation
plt.figure()
for i in range(len(offsetDegSepNormPop)):
    params = gaussFit(offsetDegSepNormPop[i], transectNormalized[i])
    xFull = np.linspace(offsetDegSepNormPop[i][0], offsetDegSepNormPop[i][-1], 1000)
    respFull = gauss(xFull, *params)
    respPred = gauss(offsetDegSepNormPop[i], *params)
    plt.plot(xFull, respFull)
    r2 = r2_score(transectNormalized[i], respPred)
    print(r2)
plt.show()

"""

fig, ax2 = plt.subplots()
fig.set_size_inches(10, 8)

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


# # rf weighted average
# params = gaussFit(x, prefMean)
# xFull = np.linspace(x[0], x[-1], 1000)
# respFull = gauss(xFull, *params)
# # pCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / (
# #             ((prefMean[0] + prefMean[1:])**1) / np.max(respFull))
pCentPred = ((prefMean[0] + nonprefMean[1:])**1) / (
            1 + prefMean[1:] / np.max(prefMean))

# # simple average
# pCentPred = ((prefMean[0] + nonprefMean[1:]) ** exp) / 2

ax2.plot(x, prefMean, color='black', label='PO')
ax2.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
             color='black', markersize=2)
ax2.plot(x, nonprefMean, color='grey', label='NO')
ax2.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
             color='grey', markersize=2)
ax2.plot(xNorm, pnMean, color='green', label='P0 N1')
ax2.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
             color='green', markersize=2)

# ax2.plot(xNorm, pCentPred, color='green', linestyle='--')
# ax2.scatter(xNorm, pCentPred, color='green')

ax2.set_xlabel('Stimulus Offset Position', fontsize=20)
ax2.set_ylabel('Normalized Response', fontsize=20)
ax2.set_ylim([0, 1])

ax2.set_xticks([0, 1, 2, 3, 4])
ax2.set_yticks([0, 0.25, 0.50, 0.75, 1])
ax2.set_xticklabels([0, 1, 2, 3, 4], fontsize=15)
ax2.set_yticklabels([0, 0.25, 0.50, 0.75, 1], fontsize=15)

plt.show()
