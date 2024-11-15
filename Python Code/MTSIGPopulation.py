"""
MTSigma population main analysis script. This script will extract spike counts during
valid trials from units active in a session and put them in the appropriate
matrix. Then we will add this to a population matrix to measure how c50 changes
with rMax at a population level.

Chery - July 2024

Notes: 240910 and 240911 are the same unit (just labelled different files to help with analysis)
       240912 and 240913 are the same units (240913 has better stimulus optimization to drive more units in an optimized way)
"""

###### NEED TO ADD UNITS from KS4 sessions #######

# Imports
from usefulFns import *
from normalizationFunctions import *
import numpy as np
import psignifit as ps
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
import warnings

t0 = time.time()

fileList = ['Akshan_240610', 'Akshan_240701', 'Akshan_240903', 'Akshan_240904',
            'Akshan_240906', 'Akshan_240911', 'Akshan_240913', 'Akshan_240916',
            'Akshan_240917', 'Akshan_240922', 'Akshan_240924', 'Akshan_241029',
            'Akshan_241104', 'Akshan_241112']

# # unit list with 50x difference b/w CI (for inclusion criteria) 241104, 240913 resorted
unitList = ['240610_92', '240701_5', '240903_3', '240904_8', '240906_25',
            '240911_17', '240913_0', '240913_4', '240913_13', '240913_44',
            '240916_7', '240917_44', '240922_6', '240924_21', '240924_35',
            '241104_3', '241104_4', '241104_23', '241112_18', '241112_19']


# # # unit list with 50x difference b/w CI (for inclusion criteria);
# unitList = ['240610_92', '240701_5', '240903_3', '240904_8', '240906_25',
#             '240911_17', '240913_47', '240916_7', '240917_44', '240922_6',
#             '240924_21', '240924_35', '241104_41', '241104_43', '241112_18',
#             '241112_19']

# # unit list with 200x difference b/w CI (for inclusion criteria)
# unitList = ['240610_92', '240701_5', '240903_3', '240904_8', '240906_25',
#             '240911_17', '240913_47', '240916_7', '240916_13', '240917_44',
#             '240922_6', '240924_21', '240924_35', '241104_41', '241104_43']


# # non-curated file List Akshan (C50 at center and peri for pref unmatched)
# fileList = ['Akshan_240603', 'Akshan_240606', 'Akshan_240610', 'Akshan_240701',
#             'Akshan_240903', 'Akshan_240904', 'Akshan_240905', 'Akshan_240906',
#             'Akshan_240909', 'Akshan_240911', 'Akshan_240913', 'Akshan_240916',
#             'Akshan_240917', 'Akshan_240922', 'Akshan_240924', 'Akshan_241029',
#             'Akshan_241104']
#
# # # kilosort 4 list
# unitList = ['240603_118', '240606_46', '240606_99', '240610_169', '240701_1',
#             '240903_3', '240903_8_2', '240905_4', '240906_25', '240909_10',
#             '240911_17', '240913_1', '240913_27', '240913_36', '240913_47',
#             '240913_48', '240916_7', '240917_44', '240922_6', '240924_22',
#             '240924_35', '240924_21', '241029_25', '241104_25', '241104_35']

### OLD WORKING VERSION
# # # kilosort 2 list
# unitList = ['240603_167', '240606_176', '240610_169', '240701_1', '240903_3',
#             '240903_8_2', '240905_4', '240906_25', '240909_10', '240911_17',
#             '240913_1', '240913_27', '240913_36', '240913_47', '240913_48',
#             '240916_7', '240917_44', '240922_6', '240924_22', '240924_35',
#             '240924_21', '241029_25', '241104_25', '241104_35']

# excluded one unit
# unitList = ['240603_167', '240606_176', '240610_169', '240701_1', '240903_3',
#             '240903_8_2', '240905_4', '240906_25', '240909_10', '240911_17',
#             '240913_27', '240913_36', '240913_47', '240913_48',
#             '240916_7', '240917_44', '240922_6', '240924_22', '240924_35',
#             '240924_21', '241029_25']


# # curated file List Akshan (C50 at center and peri for pref unmatched)
# fileList = ['Akshan_240603', 'Akshan_240606', 'Akshan_240610', 'Akshan_240903_2',
#             'Akshan_240905', 'Akshan_240906', 'Akshan_240913']
# unitList = ['240603_167', '240606_176', '240610_169', '240903_8_2', '240905_4',
#             '240906_25', '240913_27', '240913_48', '240913_36']

# fileList = ['Akshan_240701']
# unitList = ['240701_1']

rMaxPerChange = []
c50PerChange = []
prefCenterPopResp = []
prefPeriPopResp = []
nonprefCenterPopResp = []
nonprefPeriPopResp = []
centerC50Pop = []
periC50Pop = []
centerRMaxNormPop = []
periRMaxNormPop = []
nonprefCenterC50Pop = []
nonprefPeriC50Pop = []
nonprefCenterRMaxNormPop = []
nonprefPeriRMaxNormPop = []
badUnits = []
popBlocksDone = []
popUnitCluster = []

for file in fileList:

    # Load relevant file here with pyMat reader
    if len(file.split('_')) == 2:
        monkeyName, seshDate = file.split('_')
        fileName = f'{monkeyName}_{seshDate}_MTSIG_Spikes.mat'
        dayFileNumber = 1
    else:
        monkeyName, seshDate, dayFileNumber = file.split('_')
        dayFileNumber = int(dayFileNumber)
        fileName = f'{monkeyName}_{seshDate}_MTSIG_unit2_Spikes.mat'
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
        if len(unit.split('_')) == 2 and dayFileNumber == 1:
            date, unitID = unit.split('_')
            if date == seshDate:
                sessionGoodUnits.append(int(unitID))
        if len(unit.split('_')) == 3 and dayFileNumber == 2:
            date, unitID, filler = unit.split('_')
            if date == seshDate:
                sessionGoodUnits.append(int(unitID))
    sessionGoodUnits = np.array(sessionGoodUnits)

    # change stimDesc to be list of dictionaries
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        nStim = len(currTrial['stimDesc']['data']['locType'])
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
            if stim['locType'] == 0:
                frameDiff = stim['stimOffFrame'] - stim['stimOnFrame']
                stimDurFrame.append(frameDiff)

    if len(set(stimDurFrame)) != 1:
        print('stimulus frame duration not consistent for mapping stimuli')
    else:
        trueStimDurMS = np.int32(np.around(1000 / frameRateHz * stimDurFrame[0]))

    # initialize lists/arrays/dataframes for counting spikeCounts and for analysis
    if 'blockStatus' in allTrials[corrTrials[-1]]:
        blocksDone = allTrials[corrTrials[-1]]['blockStatus']['data']['blocksDone'][0]
    else:
        blocksDone = allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone'][0]
    numContrasts = header['blockStatus']['data']['numContrasts'][0]
    contrasts = header['blockStatus']['data']['contrasts'][0] * 100
    contrasts = np.insert(contrasts, 0, 0)
    spikeCountMat = np.zeros((len(units), blocksDone+1, numContrasts*4))
    stimCount = np.zeros((2, 2, numContrasts), dtype=int)
    sponRate = np.zeros((len(units), numContrasts*4*(blocksDone+1)))
    sponIndex = 0
    stimCountIndex = np.arange(numContrasts*4)
    stimCountIndex = stimCountIndex.reshape(2, 2, numContrasts)
    onLatency = 50 / 1000  # time in MS for counting window latency after stim on
    offLatency = 50 / 1000  # time in MS for counting window latency after stim off
    histPrePostMS = 100  # 100ms window pre/post stimulus on/off
    spikeHists = np.zeros((len(units), numContrasts*4, trueStimDurMS + (2*histPrePostMS+1)))

    # insert spikes from valid stimulus presentations into spike count matrices
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        if 'spikeData' in currTrial:
            stimDesc = currTrial['stimDesc']['data']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            for stim in stimDesc:
                if stim['stimType'] == 1 and stim['locType'] != 2:
                    stimOnTimeS = ((1000 / frameRateHz * stim['stimOnFrame'])
                                   / 1000) + stim1TimeS
                    stimOffTimeS = ((1000 / frameRateHz * stim['stimOffFrame'])
                                    / 1000) + stim1TimeS
                    stimLocation = stim['locType']
                    stimDirection = stim['dirType']
                    stimContrast = stim['contrastIndex']
                    stCount = int(stimCount[stimLocation, stimDirection, stimContrast])
                    stimCount[stimLocation, stimDirection, stimContrast] += 1
                    stimIndex = stimCountIndex[stimLocation, stimDirection, stimContrast]

                    for unitCount, unit in enumerate(units):
                        if unit in currTrial['spikeData']['unit']:
                            unitIndex = np.where(currTrial['spikeData']['unit'] == unit)[0]
                            unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                            stimSpikes = np.where((unitTimeStamps >= (stimOnTimeS + onLatency)) &
                                                  (unitTimeStamps <= (stimOffTimeS + offLatency)))[0]
                            spikeCountMat[unitCount][stCount][stimIndex] \
                                = len(stimSpikes)
                            sponRate[unitCount][sponIndex] = len(np.where((unitTimeStamps >= (stimOnTimeS - (100/1000))) &
                                                                          (unitTimeStamps <= stimOnTimeS))[0])

                            # PSTHs
                            stimOnPreSNEV = stimOnTimeS - (histPrePostMS / 1000)
                            stimOffPostSNEV = stimOffTimeS + (histPrePostMS / 1000)
                            histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                                             & (unitTimeStamps <= stimOffPostSNEV)
                                                             )] - stimOnPreSNEV
                            histStimSpikes = np.int32(histStimSpikes * 1000)
                            spikeHists[unitCount, stimIndex, histStimSpikes] += 1

                    sponIndex += 1

    # mean, SEM, and reshaping of spikeCount matrices (if pref and non pref are different)
    meanSpike = np.mean(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSD = np.std(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSEM = spikeCountSD/np.sqrt(blocksDone)
    meanSpikeReshaped = np.zeros((len(units), 2, 2, numContrasts))
    SEMReshaped = np.zeros((len(units), 2, 2, numContrasts))
    for count, i in enumerate(meanSpikeReshaped):
        meanSpikeReshaped[count] = (meanSpike[count].reshape(2, 2, numContrasts) *
                                    1000/trueStimDurMS)
        SEMReshaped[count] = (spikeCountSEM[count].reshape(2, 2, numContrasts) *
                              1000/trueStimDurMS)


    # some sessions pref/nonpref and cent/peri need to be flipped as stimuli were optimized
    # for other neurons from that session
    if seshDate == '241112':
        indx = np.where(units == 18)[0][0]
        meanSpikeReshaped[indx] = meanSpikeReshaped[indx][::-1, ::-1, :]
    if seshDate == '241107':
        indx = np.where(units == 6)[0][0]
        meanSpikeReshaped[indx] = meanSpikeReshaped[indx][::-1, :, :]
        indx = np.where(units == 16)[0][0]
        meanSpikeReshaped[indx] = meanSpikeReshaped[indx][::-1, :, :]

    # fit CRF
    for unit in sessionGoodUnits:
        count = np.where(units == unit)[0][0]
        popUnitCluster.append(unitCluster[count])
        baselineResp = np.mean(sponRate[count][:numContrasts * 4 * blocksDone]) * 1000 / trueStimDurMS
        warnings.simplefilter('ignore', OptimizeWarning)

        unitRMax = []
        unitC50 = []
        centPeriRMax = []
        centPeriC50 = []
        nonprefCentPeriRMax = []
        nonprefCentPeriC50 = []

        for i in range(2):
            for j in range(2):
                response = meanSpikeReshaped[count][i][j]
                response = np.insert(response, 0, baselineResp)
                try:
                    initialGuess = [baselineResp, max(response), np.median(contrasts), 2.0]
                    pOpt, pCov = curve_fit(contrastFn, contrasts, response,
                                           bounds=([baselineResp, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
                    # driven rate
                    unitRMax.append(pOpt[1]-pOpt[0])
                    unitC50.append(pOpt[2])
                    if j == 0:
                        centPeriRMax.append(pOpt[1])
                        centPeriC50.append(pOpt[2])
                    if j == 1:
                        nonprefCentPeriRMax.append(pOpt[1])
                        nonprefCentPeriC50.append(pOpt[2])
                except (RuntimeError, ValueError) as e:
                    print('no good fit found')
                    print(seshDate, unit, file)

        prefRMaxChange = (unitRMax[2] - unitRMax[0]) / unitRMax[0] * 100
        prefC50Change = 1 / ((unitC50[2] - unitC50[0]) / unitC50[0]) * 100
        # npRMaxChange = (unitRMax[1] - unitRMax[3]) / unitRMax[1] * 100
        # npC50Change = 1 / ((unitC50[3] - unitC50[1]) / unitC50[1]) * 100

        # create population avg CRF curve for pref resp at center and peri
        # response matrix is baseline subtracted and normalized to center resp
        # at 100% contrast
        prefCenterResp = meanSpikeReshaped[count][0][0]
        prefCenterResp = np.insert(prefCenterResp, 0, baselineResp)
        prefCenterResp = prefCenterResp - baselineResp
        prefPeriResp = meanSpikeReshaped[count][1][0]
        prefPeriResp = np.insert(prefPeriResp, 0, baselineResp)
        prefPeriResp = prefPeriResp - baselineResp
        prefCenterPopResp.append(prefCenterResp/prefCenterResp[-1])
        prefPeriPopResp.append(prefPeriResp/prefCenterResp[-1])

        # create population avg CRF curve for nonpref resp at center and peri
        # response matrix is baseline subtracted and normalized to center resp
        # at 100% contrast
        nonprefCenterResp = meanSpikeReshaped[count][0][1]
        nonprefCenterResp = np.insert(nonprefCenterResp, 0, baselineResp)
        nonprefCenterResp = nonprefCenterResp - baselineResp
        nonprefPeriResp = meanSpikeReshaped[count][1][1]
        nonprefPeriResp = np.insert(nonprefPeriResp, 0, baselineResp)
        nonprefPeriResp = nonprefPeriResp - baselineResp
        nonprefCenterPopResp.append(nonprefCenterResp/prefCenterResp[-1])
        nonprefPeriPopResp.append(nonprefPeriResp/prefCenterResp[-1])

        # add normalized rMax and c50 values for pref direction at center/peri to population lists
        centPeriRMax = centPeriRMax/np.max(centPeriRMax)
        centerC50Pop.append(centPeriC50[0])
        periC50Pop.append(centPeriC50[1])
        centerRMaxNormPop.append(centPeriRMax[0])
        periRMaxNormPop.append(centPeriRMax[1])

        # add normalized rMax and c50 values for nonpref direction at center/peri to population lists
        nonprefCentPeriRMax = nonprefCentPeriRMax/np.max(nonprefCentPeriRMax)
        nonprefCenterC50Pop.append(nonprefCentPeriC50[0])
        nonprefPeriC50Pop.append(nonprefCentPeriC50[1])
        nonprefCenterRMaxNormPop.append(nonprefCentPeriRMax[0])
        nonprefPeriRMaxNormPop.append(nonprefCentPeriRMax[1])

        popBlocksDone.append(blocksDone)

        rMaxPerChange.append(prefRMaxChange)
        c50PerChange.append(prefC50Change)
        # rMaxPerChange.append(npRMaxChange)
        # c50PerChange.append(npC50Change)

        # bootstrap the naka-rusthon function to find CI of c50
        prefCentRawSpikeCounts = spikeCountMat[count][:blocksDone, :6] * 1000/trueStimDurMS
        prefCentRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
                                                    prefCentRawSpikeCounts))
        prefCentResp = prefCentRawSpikeCountsWithBase - baselineResp

        nonprefCentRawSpikeCounts = spikeCountMat[count][:blocksDone, 6:12] * 1000/trueStimDurMS
        nonprefCentRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
                                                      nonprefCentRawSpikeCounts))
        nonprefCentResp = nonprefCentRawSpikeCountsWithBase - baselineResp

        prefPeriRawSpikeCounts = spikeCountMat[count][:blocksDone, 12:18] * 1000/trueStimDurMS
        prefPeriRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
                                                    prefPeriRawSpikeCounts))
        prefPeriResp = prefPeriRawSpikeCountsWithBase - baselineResp

        nonprefPeriRawSpikeCounts = spikeCountMat[count][:blocksDone, 18:] * 1000/trueStimDurMS
        nonprefPeriRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
                                                      nonprefPeriRawSpikeCounts))
        nonprefPeriResp = nonprefPeriRawSpikeCountsWithBase - baselineResp

        n_bootstraps = 1000
        c50bootstrapPrefCent = []
        c50bootstrapNonprefCent = []
        c50bootstrapPrefPeri = []
        c50bootstrapNonprefPeri = []
        # Bootstrap loop
        # for _ in range(n_bootstraps):
        #     # pref center
        #     # Sample rows with replacement to get a bootstrap sample
        #     bootstrap_sample = prefCentResp[
        #         np.random.choice(prefCentResp.shape[0],
        #                          size=prefCentResp.shape[0],
        #                          replace=True)]
        #     meanCount = np.mean(bootstrap_sample, axis=0)
        #     pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
        #                         bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
        #                         maxfev=10000000)
        #     c50bootstrapPrefCent.append(pOpt[1])
        #
        #     # nonpref cent
        #     bootstrap_sample = nonprefCentResp[
        #         np.random.choice(nonprefCentResp.shape[0],
        #                          size=nonprefCentResp.shape[0],
        #                          replace=True)]
        #     meanCount = np.mean(bootstrap_sample, axis=0)
        #     pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
        #                         bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
        #                         maxfev=10000000)
        #     c50bootstrapNonprefCent.append(pOpt[1])
        #
        #     # pref peri
        #     bootstrap_sample = prefPeriResp[
        #         np.random.choice(prefPeriResp.shape[0],
        #                          size=prefPeriResp.shape[0],
        #                          replace=True)]
        #     meanCount = np.mean(bootstrap_sample, axis=0)
        #     pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
        #                         bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
        #                         maxfev=10000000)
        #     c50bootstrapPrefPeri.append(pOpt[1])
        #
        #     # nonpref peri
        #     bootstrap_sample = nonprefPeriResp[
        #         np.random.choice(nonprefPeriResp.shape[0],
        #                          size=nonprefPeriResp.shape[0],
        #                          replace=True)]
        #     meanCount = np.mean(bootstrap_sample, axis=0)
        #     pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
        #                         bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
        #                         maxfev=10000000)
        #     c50bootstrapNonprefPeri.append(pOpt[1])
        #
        # c50bootstrapPrefCent = np.array(c50bootstrapPrefCent)
        # c50bootstrapNonprefCent = np.array(c50bootstrapNonprefCent)
        # c50bootstrapPrefPeri = np.array(c50bootstrapPrefPeri)
        # c50bootstrapNonprefPeri = np.array(c50bootstrapNonprefPeri)
        #
        # # Calculate confidence interval (e.g., 95%)
        # confidence_interval = np.percentile(c50bootstrapPrefCent, [2.5, 97.5])
        #
        # prefCentCI = np.percentile(c50bootstrapPrefCent, [2.5, 97.5])
        # nonprefCentCI = np.percentile(c50bootstrapNonprefCent, [2.5, 97.5])
        # prefPeriCI = np.percentile(c50bootstrapPrefPeri, [2.5, 97.5])
        # nonprefPeriCI = np.percentile(c50bootstrapNonprefPeri, [2.5, 97.5])
        #
        # randVar = 0
        # if prefCentCI[1]/prefCentCI[0] > 50:
        #     randVar = 1
        # if nonprefCentCI[1]/nonprefCentCI[0] > 50:
        #     randVar = 1
        # if prefPeriCI[1]/prefPeriCI[0] > 50:
        #     randVar = 1
        # if nonprefPeriCI[1]/nonprefPeriCI[0] > 50:
        #     randVar = 1
        #
        # if randVar == 1:
        #     date = file.split('_')[-1]
        #     badU = f'{date}_{unit}'
        #     badUnits.append(badU)

    # to open another file in the loop
    os.chdir('../../Python Code')

print(time.time()-t0)

#################################################### Plotting ####################################################
unitList = np.array(unitList)
badUnits = np.array(badUnits)

# population plots
prefCenterPopResp = np.array(prefCenterPopResp)
prefPeriPopResp = np.array(prefPeriPopResp)
nonprefCenterPopResp = np.array(nonprefCenterPopResp)
nonprefPeriPopResp = np.array(nonprefPeriPopResp)

prefCenterMean = np.mean(prefCenterPopResp, axis=0)
prefPeriMean = np.mean(prefPeriPopResp, axis=0)
prefCenterSEM = stats.sem(prefCenterPopResp, axis=0)
prefPeriSEM = stats.sem(prefPeriPopResp, axis=0)

nonprefCenterMean = np.mean(nonprefCenterPopResp, axis=0)
nonprefPeriMean = np.mean(nonprefPeriPopResp, axis=0)
nonprefCenterSEM = stats.sem(nonprefCenterPopResp, axis=0)
nonprefPeriSEM = stats.sem(nonprefPeriPopResp, axis=0)

# nonmatches = ~np.isin(unitList, badUnits)
# goodUnitsIndx = np.where(nonmatches)[0]
#
# prefCenterMean = np.mean(prefCenterPopResp[goodUnitsIndx], axis=0)
# prefPeriMean = np.mean(prefPeriPopResp[goodUnitsIndx], axis=0)
# prefCenterSEM = stats.sem(prefCenterPopResp[goodUnitsIndx], axis=0)
# prefPeriSEM = stats.sem(prefPeriPopResp[goodUnitsIndx], axis=0)
#
# nonprefCenterMean = np.mean(nonprefCenterPopResp[goodUnitsIndx], axis=0)
# nonprefPeriMean = np.mean(nonprefPeriPopResp[goodUnitsIndx], axis=0)
# nonprefCenterSEM = stats.sem(nonprefCenterPopResp[goodUnitsIndx], axis=0)
# nonprefPeriSEM = stats.sem(nonprefPeriPopResp[goodUnitsIndx], axis=0)

hfont = {'fontname': 'Arial'}

# plot pref and non pref center CRF using Naka-Rushton
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')

    x_min = ax1.get_xlim()[0]
    initialGuess = [max(prefCenterMean), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFnNoBaseline, contrasts, prefCenterMean,
                           bounds=([0, 0, 0],
                                   [np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFnNoBaseline(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='black', label=f'c50={pOpt[1]:.2f}, rMax={pOpt[0]:.2f}')

    # Add truncated lines for pref curve
    c50_pref = pOpt[1]
    rMax_pref = pOpt[0]
    half_response_pref = 0.5 * rMax_pref
    ax1.plot([c50_pref, c50_pref], [0, half_response_pref],
             color='black', linestyle=':', label='Pref c50')
    ax1.plot([x_min, c50_pref], [half_response_pref, half_response_pref],
             color='black', linestyle=':')

    # plot nonpref Center
    initialGuess = [max(nonprefCenterMean), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFnNoBaseline, contrasts, nonprefCenterMean,
                           bounds=([0, 0, 0],
                                   [np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFnNoBaseline(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='gray', label=f'c50={pOpt[1]:.2f}, rMax={pOpt[0]:.2f}')

    # Add truncated lines for non-pref curve
    c50_nonpref = pOpt[1]
    rMax_nonpref = pOpt[0]
    half_response_nonpref = 0.5 * rMax_nonpref
    ax1.plot([c50_nonpref, c50_nonpref], [0, half_response_nonpref],
             color='gray', linestyle=':', label='Non-pref c50')
    ax1.plot([x_min, c50_nonpref], [half_response_nonpref, half_response_nonpref],
             color='gray', linestyle=':')

    ax1.errorbar(contrasts, prefCenterMean, yerr=prefCenterSEM, fmt='o', color='black')
    ax1.errorbar(contrasts, nonprefCenterMean, yerr=nonprefCenterSEM, fmt='o', color='gray')
    ax1.set_xscale('symlog', linthresh=0.1)
    ax1.set_xlabel('Contrast (%)')
    ax1.set_ylabel('Normalized Spike Rate')
    ax1.set_title(f'Population Contrast Response Function n={len(prefCenterPopResp)}')
    ax1.legend()
    ax1.set_ylim([0, 1.1])
    sns.despine(offset=5)

    plt.show()

# plot pref and non pref peri CRF using Naka-Rushton
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')

    # plot pref peri
    initialGuess = [max(prefPeriMean), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFnNoBaseline, contrasts, prefPeriMean,
                           bounds=([0, 0, 0],
                                   [np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFnNoBaseline(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='black',
             label=f'c50={pOpt[1]:.2f}, rMax={pOpt[0]:.2f}')

    # Add truncated lines for pref curve
    c50_pref = pOpt[1]
    rMax_pref = pOpt[0]
    half_response_pref = 0.5 * rMax_pref
    ax1.plot([c50_pref, c50_pref], [0, half_response_pref],
             color='black', linestyle=':', label='Pref c50')
    ax1.plot([0, c50_pref], [half_response_pref, half_response_pref],
             color='black', linestyle=':')

    # plot nonpref peri
    initialGuess = [max(nonprefPeriMean), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFnNoBaseline, contrasts, nonprefPeriMean,
                           bounds=([0, 0, 0],
                                   [np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFnNoBaseline(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='gray',
             label=f'c50={pOpt[1]:.2f}, rMax={pOpt[0]:.2f}')

    # Add truncated lines for non-pref curve
    c50_nonpref = pOpt[1]
    rMax_nonpref = pOpt[0]
    half_response_nonpref = 0.5 * rMax_nonpref
    ax1.plot([c50_nonpref, c50_nonpref], [0, half_response_nonpref],
             color='gray', linestyle=':', label='Non-pref c50')
    ax1.plot([0, c50_nonpref], [half_response_nonpref, half_response_nonpref],
             color='gray', linestyle=':')

    ax1.errorbar(contrasts, prefPeriMean, yerr=prefPeriSEM,
                 fmt='o', color='black')
    ax1.errorbar(contrasts, nonprefPeriMean, yerr=nonprefPeriSEM,
                 fmt='o', color='gray')

    ax1.set_xscale('symlog', linthresh=0.1)
    ax1.set_xlabel('Contrast (%)')
    ax1.set_ylabel('Normalized Spike Rate')
    ax1.set_title(f'Population Contrast Response Function n={len(prefCenterPopResp)}')
    ax1.legend()
    ax1.set_ylim([0, 1.1])
    sns.despine(offset=5)

    plt.show()

# plot box plots of center/peri pref/non-pref fit c50 naka-rushton
centerPrefPopC50 = []
centerPrefPopRMax = []
centerNonprefPopC50 = []
centerNonprefPopRMax = []
periPrefPopC50 = []
periPrefPopRMax = []
periNonprefPopC50 = []
for i in range(len(prefCenterPopResp)):
    # pref center
    initialGuess = [max(prefCenterPopResp[i]), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFnNoBaseline, contrasts, prefCenterPopResp[i],
                           bounds=([0, 0, 0],
                                   [np.inf, np.inf, np.inf]),
                           p0=initialGuess)
    centerPrefPopC50.append(pOpt[1])
    centerPrefPopRMax.append(pOpt[0])

    # nonpref center
    initialGuess = [max(nonprefCenterPopResp[i]), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFnNoBaseline, contrasts, nonprefCenterPopResp[i],
                           bounds=([0, 0, 0],
                                   [np.inf, np.inf, np.inf]),
                           p0=initialGuess)
    centerNonprefPopC50.append(pOpt[1])
    centerNonprefPopRMax.append(pOpt[0])

    # pref peri
    initialGuess = [max(prefPeriPopResp[i]), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFnNoBaseline, contrasts, prefPeriPopResp[i],
                           bounds=([0, 0, 0],
                                   [np.inf, np.inf, np.inf]),
                           p0=initialGuess)
    periPrefPopC50.append(pOpt[1])
    periPrefPopRMax.append(pOpt[0])

    # nonpref peri
    initialGuess = [max(nonprefPeriPopResp[i]), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFnNoBaseline, contrasts, nonprefPeriPopResp[i],
                           bounds=([0, 0, 0],
                                   [np.inf, np.inf, np.inf]),
                           p0=initialGuess)
    periNonprefPopC50.append(pOpt[1])

centerPrefPopC50 = np.array(centerPrefPopC50)
centerPrefPopRMax = np.array(centerPrefPopRMax)
centerNonprefPopRMax = np.array(centerNonprefPopRMax)
centerNonprefPopC50 = np.array(centerNonprefPopC50)
periPrefPopRMax = np.array(periPrefPopRMax)
periPrefPopC50 = np.array(periPrefPopC50)

# stats test for box plots
# Defining only the specific pairs you want to compare
data_pairs = [
    (centerPrefPopC50, centerNonprefPopC50, "Center Pref vs Center Non-Pref"),
    (centerPrefPopC50, periPrefPopC50, "Center Pref vs Peri Pref"),
    (centerNonprefPopC50, periNonprefPopC50, "Center Non-Pref vs Peri Non-Pref"),
    (periPrefPopC50, periNonprefPopC50, "Peri Pref vs Peri Non-Pref")
]
# Bonferroni correction for 4 comparisons
alpha = 0.05
alpha_adjusted = alpha / len(data_pairs)

# Conduct Mann-Whitney U tests for each pair
for data1, data2, comparison_name in data_pairs:
    stat, p = mannwhitneyu(data1, data2, alternative='two-sided')
    # stat, p, _, _ = median_test(data1, data2)
    print(f"{comparison_name} - Test statistic: {stat}, p-value: {p}")
    if p < alpha_adjusted:
        print(f"Significant difference after Bonferroni correction (p < {alpha_adjusted})")
    else:
        print(f"No significant difference after Bonferroni correction (p >= {alpha_adjusted})")

for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')

    ax1.boxplot(centerPrefPopC50, positions=[0.8],
                widths=0.35, patch_artist=True,
                boxprops=dict(facecolor='black', color='black'),
                medianprops=dict(color='white'))

    ax1.boxplot(centerNonprefPopC50, positions=[1.2],
                widths=0.35, patch_artist=True,
                boxprops=dict(facecolor='gray', color='gray'),
                medianprops=dict(color='white'))


    ax1.boxplot(periPrefPopC50, positions=[2.8],
                widths=0.35, patch_artist=True,
                boxprops=dict(facecolor='black', color='black'),
                medianprops=dict(color='white'))

    ax1.boxplot(periNonprefPopC50, positions=[3.2],
                widths=0.35, patch_artist=True,
                boxprops=dict(facecolor='gray', color='gray'),
                medianprops=dict(color='white'))

    ax1.set_xticks([1, 3])
    ax1.set_xticklabels(['Center', 'Peripheral'], **hfont, fontsize=20)
    ax1.set_ylabel('Fit C50', **hfont, fontsize=25)
    ax1.set_xlabel('Stimulus Location', **hfont, fontsize=25)
    ax1.set_ylim([0, 100])
    ax1.tick_params(axis='both', width=2, length=8)
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_ticks = ax1.get_yticks()
    y_tick_labels = [
        '0' if tick == 0 else
        '100' if tick == 100 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)

    # # Draw significance lines and stars across groups
    # y_max = max(max(centerPrefPopC50), max(centerNonprefPopC50), max(periPrefPopC50), max(periNonprefPopC50))
    # line_height = y_max + 5  # height of the line for significance
    # star_height = line_height + 0.02  # height of stars
    #
    # # Draw line from "Center" to "Peripheral" group
    # ax1.plot([1, 3], [line_height, line_height], color='black', linewidth=1.5)
    # ax1.text(2, star_height, '***', ha='center', va='bottom', color='black', fontsize=16)

    plt.show()

# plot c50 as a function of rMax
rmaxCentNPtoCentP = centerNonprefPopRMax/centerPrefPopRMax
rmaxPeriPtoCentP = periPrefPopRMax/centerPrefPopRMax
c50CentNPtoCentP = centerNonprefPopC50/centerPrefPopC50
c50PeriPtoCentP = periPrefPopC50/centerPrefPopC50

fig, ax1 = plt.subplots()
ax1.scatter(rmaxCentNPtoCentP, c50CentNPtoCentP, color='orange')
ax1.scatter(rmaxPeriPtoCentP, c50PeriPtoCentP, color='green')
plt.show()

stats.pearsonr(c50CentNPtoCentP, rmaxCentNPtoCentP)
stats.pearsonr(c50PeriPtoCentP, rmaxPeriPtoCentP)

# plot distributions of concatenated center vs peri
a = np.concatenate((centerPrefPopC50, centerNonprefPopC50), axis=0)
b = np.concatenate((periPrefPopC50, periNonprefPopC50), axis=0)

fig, ax = plt.subplots()
ax.hist(b, color="orange", alpha=0.5, edgecolor="black", label="Distribution B")
ax.hist(a, color="blue", alpha=0.5, edgecolor="black", label="Distribution A")
ax.legend()  # Add a legend to label each distribution
plt.show()


# bootstrap CI for population plots
for filler in range(1):
    n_bootstraps = 1000
    c50bootstrapPrefCent = []
    c50bootstrapNonprefCent = []
    c50bootstrapPrefPeri = []
    c50bootstrapNonprefPeri = []

    # Bootstrap loop
    for _ in range(n_bootstraps):
        # pref center
        # Sample rows with replacement to get a bootstrap sample
        bootstrap_sample = prefCenterPopResp[
            np.random.choice(prefCenterPopResp.shape[0],
                             size=prefCenterPopResp.shape[0],
                             replace=True)]
        meanCount = np.mean(bootstrap_sample, axis=0)
        pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
                            bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
                            maxfev=10000000)
        c50bootstrapPrefCent.append(pOpt[1])

        # nonpref cent
        bootstrap_sample = nonprefCenterPopResp[
            np.random.choice(nonprefCenterPopResp.shape[0],
                             size=nonprefCenterPopResp.shape[0],
                             replace=True)]
        meanCount = np.mean(bootstrap_sample, axis=0)
        pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
                            bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
                            maxfev=10000000)
        c50bootstrapNonprefCent.append(pOpt[1])

        # pref peri
        bootstrap_sample = prefPeriPopResp[
            np.random.choice(prefPeriPopResp.shape[0],
                             size=prefPeriPopResp.shape[0],
                             replace=True)]
        meanCount = np.mean(bootstrap_sample, axis=0)
        pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
                            bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
                            maxfev=10000000)
        c50bootstrapPrefPeri.append(pOpt[1])

        # nonpref peri
        bootstrap_sample = nonprefPeriPopResp[
            np.random.choice(nonprefPeriPopResp.shape[0],
                             size=nonprefPeriPopResp.shape[0],
                             replace=True)]
        meanCount = np.mean(bootstrap_sample, axis=0)
        pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
                            bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
                            maxfev=10000000)
        c50bootstrapNonprefPeri.append(pOpt[1])

    c50bootstrapPrefCent = np.array(c50bootstrapPrefCent)
    c50bootstrapNonprefCent = np.array(c50bootstrapNonprefCent)
    c50bootstrapPrefPeri = np.array(c50bootstrapPrefPeri)
    c50bootstrapNonprefPeri = np.array(c50bootstrapNonprefPeri)

    # Calculate confidence interval (e.g., 95%)
    prefCentCI = np.percentile(c50bootstrapPrefCent, [2.5, 97.5])
    nonprefCentCI = np.percentile(c50bootstrapNonprefCent, [2.5, 97.5])
    prefPeriCI = np.percentile(c50bootstrapPrefPeri, [2.5, 97.5])
    nonprefPeriCI = np.percentile(c50bootstrapNonprefPeri, [2.5, 97.5])

    print(prefCentCI)
    print(nonprefCentCI)
    print(prefPeriCI)
    print(nonprefPeriCI)


####### extra figures
# plotting figure super figure (3 subfigures)
for filler in range(1):
    fig, ax = plt.subplots(1, 3, figsize=(14, 5))

    # first subplot
    # plot pref Center
    initialGuess = [prefCenterMean[0], max(prefCenterMean), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFn, contrasts, prefCenterMean,
                           bounds=([prefCenterMean[0], 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax[0].plot(xFit, yFit, color='green', label=f'c50={pOpt[2]:.2f}, rMax={pOpt[1]:.2f}')

    # plot pref peri
    initialGuess = [prefPeriMean[0], max(prefPeriMean), np.median(contrasts), 2.0]
    # pOpt, pCov = curve_fit(contrastFn, contrasts, prefPeriMean, bounds=([prefPeriMean[0], 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
    pOpt, pCov = curve_fit(contrastFn, contrasts, prefPeriMean,
                           bounds=([prefPeriMean[0], 0.8812, 0, 0], [np.inf, 0.8812+0.001, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax[0].plot(xFit, yFit, color='green', alpha=0.5, label=f'c50={pOpt[2]:.2f}, rMax={pOpt[1]:.2f}')

    # plot nonpref Center
    initialGuess = [prefCenterMean[0], max(nonprefCenterMean), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFn, contrasts, nonprefCenterMean,
                           bounds=([nonprefCenterMean[0], 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax[0].plot(xFit, yFit, color='red', label=f'c50={pOpt[2]:.2f}, rMax={pOpt[1]:.2f}')

    # plot nonpref peri
    initialGuess = [nonprefPeriMean[0], max(nonprefPeriMean), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFn, contrasts, nonprefPeriMean,
                           bounds=([nonprefPeriMean[0], 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax[0].plot(xFit, yFit, color='red', alpha=0.5, label=f'c50={pOpt[2]:.2f}, rMax={pOpt[1]:.2f}')

    ax[0].errorbar(contrasts, prefCenterMean, yerr=prefCenterSEM, fmt='o', color='green')
    ax[0].errorbar(contrasts, prefPeriMean, yerr=prefPeriSEM, fmt='o', color='green', alpha=0.5)
    ax[0].errorbar(contrasts, nonprefCenterMean, yerr=nonprefCenterSEM, fmt='o', color='red')
    ax[0].errorbar(contrasts, nonprefPeriMean, yerr=nonprefPeriSEM, fmt='o', color='red', alpha=0.5)

    ax[0].set_xscale('symlog', linthresh=0.1)
    ax[0].set_xlabel('Contrast (%)')
    ax[0].set_ylabel('Normalized Spike Rate')
    ax[0].set_title(f'Population Contrast Response Function n={len(prefCenterPopResp)}')
    ax[0].legend()
    sns.despine(offset=5)

    # second subplot
    # plot individual data points (scatter/boxplots)
    centerC50Pop = np.array(centerC50Pop)
    periC50Pop = np.array(periC50Pop)
    centerRMaxNormPop = np.array(centerRMaxNormPop)
    periRMaxNormPop = np.array(periRMaxNormPop)

    nonprefCenterC50Pop = np.array(nonprefCenterC50Pop)
    nonprefPeriC50Pop = np.array(nonprefPeriC50Pop)
    nonprefCenterRMaxNormPop = np.array(nonprefCenterRMaxNormPop)
    nonprefPeriRMaxNormPop = np.array(nonprefPeriRMaxNormPop)

    # 1. Add boxplots for each category
    a = np.concatenate((centerC50Pop, nonprefCenterC50Pop), axis=0)
    b = np.concatenate((periC50Pop, nonprefPeriC50Pop), axis=0)

    ax[1].boxplot([a, b], positions=[1, 2], widths=0.4, patch_artist=True,
                  boxprops=dict(facecolor='lightblue', color='blue'),
                  medianprops=dict(color='darkblue'))

    # 2. Calculate medians and plot a line connecting them
    center_median = np.median(a)
    peri_median = np.median(b)
    ax[1].plot([1, 2], [center_median, peri_median], color='black', linestyle='-', marker='o',
             label='Median Line', zorder=3)

    # 3. Create a scatter plot of the two categories
    ax[1].scatter(np.ones_like(centerC50Pop), centerC50Pop, color='green',
                  label='Pref Center C50')
    ax[1].scatter(np.ones_like(periC50Pop) * 2, periC50Pop, color='green',
                  alpha=0.5, label='Pref Peri C50')

    ax[1].scatter(np.ones_like(nonprefCenterC50Pop), nonprefCenterC50Pop, color='red',
                  label='Nonpref Center C50')
    ax[1].scatter(np.ones_like(nonprefPeriC50Pop) * 2, nonprefPeriC50Pop, color='red',
                  alpha=0.5, label='Nonpref Peri C50')

    # 4. Add lines connecting the corresponding data points between the two categories
    for i in range(len(centerC50Pop)):
        ax[1].plot([1, 2], [centerC50Pop[i], periC50Pop[i]], color='gray', linestyle='--')
    for i in range(len(nonprefCenterC50Pop)):
        ax[1].plot([1, 2], [nonprefCenterC50Pop[i], nonprefPeriC50Pop[i]], color='gray', linestyle='--')

    # 5. Add plot labels, legend, and customize x-axis
    ax[1].set_xticks([1, 2], ['Center C50', 'Peri C50'])  # Custom x-axis labels
    ax[1].set_ylabel('C50 Values')
    ax[1].set_ylim(0, 50)
    ax[1].legend()

    # plot of ratio peripheral/center C50
    ax[2].scatter(np.ones_like(nonprefCenterC50Pop) +
                np.random.uniform(low=-0.025, high=0.025, size=len(nonprefCenterC50Pop)),
                periC50Pop / centerC50Pop, color='green',
                label='Nonpref Center C50')
    ax[2].scatter(np.ones_like(nonprefCenterC50Pop) * 1.5 +
                np.random.uniform(low=-0.025, high=0.025, size=len(nonprefCenterC50Pop)),
                nonprefPeriC50Pop / nonprefCenterC50Pop, color='red',
                label='Nonpref Center C50')
    ax[2].axhline(y=1)
    ax[2].set_xlim((0.75, 1.75))
    ax[2].set_ylim(bottom=0)
    ax[2].set_ylabel('Ratio Peripheral C50 to Center C50 (pC50/cC50)', fontsize=12)
    ax[2].set_xticks([1, 1.5], ['preferred', 'non-preferred'], fontsize=12)

    plt.tight_layout()
    plt.show()

t_stat, p_value_ttest = stats.ttest_ind(a, b)

# 2. Mann-Whitney U test (non-parametric, does not assume normality)
u_stat, p_value_mannwhitney = stats.mannwhitneyu(a, b)

# 3. Kolmogorov-Smirnov test (tests whether the two distributions differ)
ks_stat, p_value_ks = stats.ks_2samp(a, b)

# Output the results
print(f"Two-sample t-test: p-value = {p_value_ttest}")
print(f"Mann-Whitney U test: p-value = {p_value_mannwhitney}")
print(f"Kolmogorov-Smirnov test: p-value = {p_value_ks}")


# box plots of c50 (naka rushton fit for individual neurons) center vs peri
for filler in range(1):
    # 1. Add boxplots for each category
    a = np.concatenate((centerC50Pop, nonprefCenterC50Pop), axis=0)
    b = np.concatenate((periC50Pop, nonprefPeriC50Pop), axis=0)

    fix, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    ax1.boxplot([a, b], positions=[1, 2], widths=0.4, patch_artist=True,
                  boxprops=dict(facecolor='lightblue', color='blue'),
                  medianprops=dict(color='darkblue'))

    # 2. Calculate medians and plot a line connecting them
    center_median = np.median(a)
    peri_median = np.median(b)
    ax1.plot([1, 2], [center_median, peri_median], color='black', linestyle='-', marker='o',
             label='Median Line', zorder=3)

    # 3. Create a scatter plot of the two categories
    ax1.scatter(np.ones_like(centerC50Pop), centerC50Pop, color='green',
                label='Pref Center C50')
    ax1.scatter(np.ones_like(periC50Pop) * 2, periC50Pop, color='green',
                alpha=0.5, label='Pref Peri C50')

    ax1.scatter(np.ones_like(nonprefCenterC50Pop), nonprefCenterC50Pop, color='red',
                label='Nonpref Center C50')
    ax1.scatter(np.ones_like(nonprefPeriC50Pop) * 2, nonprefPeriC50Pop, color='red',
                alpha=0.5, label='Nonpref Peri C50')

    # 4. Add lines connecting the corresponding data points between the two categories
    for i in range(len(centerC50Pop)):
        ax1.plot([1, 2], [centerC50Pop[i], periC50Pop[i]], color='gray', linestyle='--')
    for i in range(len(nonprefCenterC50Pop)):
        ax1.plot([1, 2], [nonprefCenterC50Pop[i], nonprefPeriC50Pop[i]],
                 color='gray', linestyle='--')


    # 5. Add plot labels, legend, and customize x-axis
    ax1.set_xticks([1, 2], ['Center C50', 'Peri C50'])  # Custom x-axis labels
    ax1.set_ylabel('C50 Values')
    ax1.set_ylim(0, 30)
    ax1.legend()

# scatter of pref/non-pref ratio c50 peri/center
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    # plot of ratio peripheral/center C50
    ax1.scatter(np.ones_like(nonprefCenterC50Pop) +
                np.random.uniform(low=-0.025, high=0.025, size=len(nonprefCenterC50Pop)),
                periC50Pop / centerC50Pop, edgecolor='black', facecolor='none',
                label='Nonpref Center C50')
    ax1.scatter(np.ones_like(nonprefCenterC50Pop) * 1.5 +
                np.random.uniform(low=-0.025, high=0.025, size=len(nonprefCenterC50Pop)),
                nonprefPeriC50Pop / nonprefCenterC50Pop, edgecolor='gray', facecolor='none',
                label='Nonpref Center C50')
    ax1.axhline(y=1)
    ax1.set_xlim((0.75, 1.75))
    ax1.set_ylim(bottom=0)
    ax1.set_ylabel('Peripheral C50/Center C50', **hfont, fontsize=25)
    ax1.set_xticks([1, 1.5], ['Preferred', 'non-preferred'], **hfont, fontsize=20)
    ax1.set_xlabel('Stimulus Identity', **hfont, fontsize=25)

    plt.show()


# a = nonprefCenterC50Pop
# b = nonprefPeriC50Pop
# 1. Two-sample t-test (assuming normality)
t_stat, p_value_ttest = stats.ttest_ind(a, b)

# 2. Mann-Whitney U test (non-parametric, does not assume normality)
u_stat, p_value_mannwhitney = stats.mannwhitneyu(a, b)

# 3. Kolmogorov-Smirnov test (tests whether the two distributions differ)
ks_stat, p_value_ks = stats.ks_2samp(a, b)

# Output the results
print(f"Two-sample t-test: p-value = {p_value_ttest}")
print(f"Mann-Whitney U test: p-value = {p_value_mannwhitney}")
print(f"Kolmogorov-Smirnov test: p-value = {p_value_ks}")


# fit simple normalization equation to get sigma
pOpt, pCov = curve_fit(normMTSIG, contrasts, prefCenterMean,
                       bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
xFit = np.logspace(-1, 2, 100)
yFit = normMTSIG(xFit, *pOpt)
plt.plot(xFit, yFit, color='green', label=f'sigma: {pOpt[1]:.2f}, L: {pOpt[0]:.2f}')
plt.scatter(contrasts, prefCenterMean, color='green')

pOpt, pCov = curve_fit(normMTSIG, contrasts, prefPeriMean,
                       bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
xFit = np.logspace(-1, 2, 100)
yFit = normMTSIG(xFit, *pOpt)
plt.plot(xFit, yFit, color='green', alpha=0.5, label=f'sigma: {pOpt[1]:.2f}, L: {pOpt[0]:.2f}')
plt.scatter(contrasts, prefPeriMean, color='green', alpha=0.5)


pOpt, pCov = curve_fit(normMTSIG, contrasts, nonprefCenterMean,
                       bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
xFit = np.logspace(-1, 2, 100)
yFit = normMTSIG(xFit, *pOpt)
plt.plot(xFit, yFit, color='red', label=f'sigma: {pOpt[1]:.2f}, L: {pOpt[0]:.2f}')
plt.scatter(contrasts, nonprefCenterMean, color='red')

pOpt, pCov = curve_fit(normMTSIG, contrasts, nonprefPeriMean,
                       bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
xFit = np.logspace(-1, 2, 100)
yFit = normMTSIG(xFit, *pOpt)
plt.plot(xFit, yFit, color='red', alpha=0.5, label=f'sigma: {pOpt[1]:.2f}, L: {pOpt[0]:.2f}')
plt.scatter(contrasts, nonprefPeriMean, color='red', alpha=0.5)

plt.xscale('symlog', linthresh=0.1)
plt.legend()

plt.show()


#### WORK STATION BELOW

# Define the conditional model
def weightedNormMTSIG(C_terms, Lp, Lnp, w, sigma, b, condition, conditionPorNP):
    result = np.zeros_like(C_terms)  # Initialize result array of same shape as C_terms
    for i, C in enumerate(C_terms):
        if conditionPorNP[i]:
            if condition[i]:
                result[i] = (C * Lp * w**1) / ((C * w) + sigma) + b
            else:
                result[i] = (C * Lp) / (C + sigma) + b
        else:
            if condition[i]:
                result[i] = (C * Lnp * w**1) / ((C * w) + sigma) + b
            else:
                result[i] = (C * Lnp) / (C + sigma) + b
    return result


def fit_data(C_terms, y_data, condition, conditionPorNP):
    # Define the objective function for curve fitting
    def objective(C_terms, Lp, Lnp, w, sigma, b):
        return weightedNormMTSIG(C_terms, Lp, Lnp, w, sigma, b, condition, conditionPorNP)

    # Initial guess for L, w, and sigma
    initial_guess = [1, .5, 0.5, 0.05, resp[0]]

    # Perform the curve fitting
    popt, pcov = curve_fit(objective, C_terms, y_data, p0=initial_guess)

    return popt, pcov


resp = np.concatenate((prefCenterMean, prefPeriMean, nonprefCenterMean, nonprefPeriMean), axis=0)
cont = np.concatenate((contrasts, contrasts, contrasts, contrasts), axis=0)
conditionArray = np.concatenate((np.full(7, False), np.full(7, True),
                                 np.full(7, False), np.full(7, True)),
                                axis=0)
condition2Array = np.concatenate((np.full(14, True), np.full(14, False)), axis=0)

# plot population RF weighted normalization (original)
for filler in range(1):
    popt, pcov = fit_data(cont, resp, conditionArray, condition2Array)
    yPred = weightedNormMTSIG(cont, *popt, conditionArray, condition2Array)
    print(r2_score(resp, yPred))

    # plotting
    xFit = np.logspace(-1, 2, 100)

    # pref center
    yFit = weightedNormMTSIG(xFit, *popt, np.full(100, False), np.full(100, True))
    plt.plot(xFit, yFit, color='green', linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[:7], *popt, conditionArray[:7])
    # plt.plot(cont[:7], yFit, color='green', linestyle='dashed')
    plt.scatter(contrasts, prefCenterMean, color='green')

    y_min, y_max = resp[0], np.max(yFit)
    y_halfway = (y_min + y_max) / 2
    # closest_index = np.argmin(np.abs(yFit - y_halfway))
    closest_index = np.argmin(np.abs(yFit - 0.5))
    x_halfway = xFit[closest_index]
    # plt.axhline(y=y_halfway, color='g')
    # plt.axvline(x=x_halfway, color='g', label=f'c50 = {x_halfway:.2f}')
    # plt.plot([xFit[0], x_halfway], [y_halfway, y_halfway], color='g')
    # plt.plot([x_halfway, x_halfway], [0, y_halfway], color='g', label=f'c50 = {x_halfway:.2f}')
    plt.plot([xFit[0], x_halfway], [0.5, 0.5], color='g')
    plt.plot([x_halfway, x_halfway], [0, 0.5], color='g', label=f'c50 = {x_halfway:.2f}')

    # pref peri
    yFit = weightedNormMTSIG(xFit, *popt, np.full(100, True), np.full(100, True))
    plt.plot(xFit, yFit, color='green', alpha=0.5, linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[7:], *popt, conditionArray[7:])
    # plt.plot(cont[7:], yFit, color='green', alpha=0.5, linestyle='dashed')
    plt.scatter(contrasts, prefPeriMean, color='green', alpha=0.5)

    y_min, y_max = resp[0], np.max(yFit)
    y_halfway = (y_min + y_max) / 2
    # closest_index = np.argmin(np.abs(yFit - y_halfway))
    closest_index = np.argmin(np.abs(yFit - 0.5))
    x_halfway = xFit[closest_index]
    # plt.axhline(y=y_halfway, color='g', alpha=0.5, linestyle='-')
    # plt.axvline(x=x_halfway, color='g', alpha=0.5, linestyle='-', label=f'c50 = {x_halfway:.2f}')
    # plt.plot([xFit[0], x_halfway], [y_halfway, y_halfway], color='g', alpha=0.5)
    # plt.plot([x_halfway, x_halfway], [0, y_halfway], color='g', alpha=0.5, label=f'c50 = {x_halfway:.2f}')
    plt.plot([xFit[0], x_halfway], [0.5, 0.5], color='g', alpha=0.5)
    plt.plot([x_halfway, x_halfway], [0, 0.5], color='g', alpha=0.5, label=f'c50 = {x_halfway:.2f}')

    # nonpref center
    yFit = weightedNormMTSIG(xFit, *popt, np.full(100, False), np.full(100, False))
    plt.plot(xFit, yFit, color='red', linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[:7], *popt, conditionArray[:7])
    # plt.plot(cont[:7], yFit, color='green', linestyle='dashed')
    plt.scatter(contrasts, nonprefCenterMean, color='red')

    y_min, y_max = resp[0], np.max(yFit)
    y_halfway = (y_min + y_max) / 2
    # closest_index = np.argmin(np.abs(yFit - y_halfway))
    closest_index = np.argmin(np.abs(yFit - 0.5))
    x_halfway = xFit[closest_index]
    plt.plot([xFit[0], x_halfway], [0.5, 0.5], color='r')
    plt.plot([x_halfway, x_halfway], [0, 0.5], color='r', label=f'c50 = {x_halfway:.2f}')

    # nonpref peri
    yFit = weightedNormMTSIG(xFit, *popt, np.full(100, True), np.full(100, False))
    plt.plot(xFit, yFit, color='red', alpha=0.5, linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[7:], *popt, conditionArray[7:])
    # plt.plot(cont[7:], yFit, color='green', alpha=0.5, linestyle='dashed')
    plt.scatter(contrasts, nonprefPeriMean, color='red', alpha=0.5)

    y_min, y_max = resp[0], np.max(yFit)
    y_halfway = (y_min + y_max) / 2
    # closest_index = np.argmin(np.abs(yFit - y_halfway))
    closest_index = np.argmin(np.abs(yFit - 0.5))
    x_halfway = xFit[closest_index]
    # plt.axhline(y=y_halfway, color='g', alpha=0.5, linestyle='-')
    # plt.axvline(x=x_halfway, color='g', alpha=0.5, linestyle='-', label=f'c50 = {x_halfway:.2f}')
    # plt.plot([xFit[0], x_halfway], [y_halfway, y_halfway], color='g', alpha=0.5)
    # plt.plot([x_halfway, x_halfway], [0, y_halfway], color='g', alpha=0.5, label=f'c50 = {x_halfway:.2f}')
    plt.plot([xFit[0], x_halfway], [0.5, 0.5], color='r', alpha=0.5)
    plt.plot([x_halfway, x_halfway], [0, 0.5], color='r', alpha=0.5, label=f'c50 = {x_halfway:.2f}')

    plt.xscale('symlog', linthresh=0.1)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(loc='lower right')
    plt.show()


# plot population RF weighted normalization (center)
for filler in range(1):
    popt, pcov = fit_data(cont, resp, conditionArray, condition2Array)
    yPred = weightedNormMTSIG(cont, *popt, conditionArray, condition2Array)
    print(r2_score(resp, yPred))

    # plotting
    xFit = np.logspace(-1, 2, 100)

    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')

    # pref center
    yFit_pref = weightedNormMTSIG(xFit, *popt, np.full(100, False), np.full(100, True))
    ax1.plot(xFit, yFit_pref, color='black', linestyle='dashed')
    ax1.scatter(contrasts, prefCenterMean, color='black')

    # nonpref center
    yFit_nonpref = weightedNormMTSIG(xFit, *popt, np.full(100, False), np.full(100, False))
    ax1.plot(xFit, yFit_nonpref, color='gray', linestyle='dashed')
    ax1.scatter(contrasts, nonprefCenterMean, color='gray')

    # Determine halfway point (50% of max response) for pref and non-pref curves
    half_y_pref = max(yFit_pref) / 2
    half_y_nonpref = max(yFit_nonpref) / 2

    # Find x-values closest to the halfway point for each curve
    x_half_pref = xFit[np.abs(yFit_pref - half_y_pref).argmin()]
    x_half_nonpref = xFit[np.abs(yFit_nonpref - half_y_nonpref).argmin()]

    # Plot truncated vertical lines
    ax1.plot([x_half_pref, x_half_pref], [0, half_y_pref], color='black',
             linestyle=':', label=f'{x_half_pref:.2f}')
    ax1.plot([x_half_nonpref, x_half_nonpref], [0, half_y_nonpref], color='gray',
             linestyle=':', label=f'{x_half_nonpref:.2f}')

    # Plot truncated horizontal lines that stop at the intersection
    ax1.plot([0, x_half_pref], [half_y_pref, half_y_pref], color='black', linestyle=':')
    ax1.plot([0, x_half_nonpref], [half_y_nonpref, half_y_nonpref], color='gray', linestyle=':')

    ax1.set_xscale('symlog', linthresh=0.1)
    ax1.set_xlim(left=0)
    ax1.set_ylim([0, 1])
    ax1.set_xlabel('Contrast (%)', **hfont, fontsize=25)
    ax1.set_ylabel('Normalized Response', **hfont, fontsize=25)

    plt.legend()

    plt.show()


# plot population RF weighted normalization (peri)
for filler in range(1):
    popt, pcov = fit_data(cont, resp, conditionArray, condition2Array)
    yPred = weightedNormMTSIG(cont, *popt, conditionArray, condition2Array)
    print(r2_score(resp, yPred))

    # plotting
    xFit = np.logspace(-1, 2, 100)

    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    # pref peri
    yFit_pref = weightedNormMTSIG(xFit, *popt, np.full(100, True), np.full(100, True))
    ax1.plot(xFit, yFit_pref, color='black', alpha=0.5, linestyle='dashed')
    ax1.scatter(contrasts, prefPeriMean, color='black', alpha=0.5)

    # nonpref peri
    yFit_nonpref = weightedNormMTSIG(xFit, *popt, np.full(100, True), np.full(100, False))
    ax1.plot(xFit, yFit_nonpref, color='gray', alpha=0.5, linestyle='dashed')
    ax1.scatter(contrasts, nonprefPeriMean, color='gray', alpha=0.5)

    # Determine halfway point (50% of max response) for pref and non-pref curves
    half_y_pref = max(yFit_pref) / 2
    half_y_nonpref = max(yFit_nonpref) / 2

    # Find x-values closest to the halfway point for each curve
    x_half_pref = xFit[np.abs(yFit_pref - half_y_pref).argmin()]
    x_half_nonpref = xFit[np.abs(yFit_nonpref - half_y_nonpref).argmin()]

    # Plot truncated vertical lines
    ax1.plot([x_half_pref, x_half_pref], [0, half_y_pref], color='black',
             linestyle=':', label=f'{x_half_pref:.2f}')
    ax1.plot([x_half_nonpref, x_half_nonpref], [0, half_y_nonpref], color='gray',
             linestyle=':', label=f'{x_half_nonpref:.2f}')

    # Plot truncated horizontal lines that stop at the intersection
    ax1.plot([0, x_half_pref], [half_y_pref, half_y_pref], color='black',
             alpha=0.5, linestyle=':')
    ax1.plot([0, x_half_nonpref], [half_y_nonpref, half_y_nonpref], color='gray',
             alpha=0.5, linestyle=':')

    ax1.set_xscale('symlog', linthresh=0.1)
    ax1.set_xlim(left=0)
    ax1.set_ylim([0, 1])
    ax1.set_xlabel('Contrast (%)', **hfont, fontsize=25)
    ax1.set_ylabel('Normalized Response', **hfont, fontsize=25)
    plt.legend()

    plt.show()


# individual neuron fits
r2Scores = []
for unit in range(len(unitList)):
    resp = np.concatenate((prefCenterPopResp[unit], prefPeriPopResp[unit],
                           nonprefCenterPopResp[unit], nonprefPeriPopResp[unit]),
                          axis=0)
    popt, pcov = fit_data(cont, resp, conditionArray, condition2Array)
    yPred = weightedNormMTSIG(cont, *popt, conditionArray, condition2Array)
    r2Scores.append((r2_score(resp, yPred)))


wInNum = np.copy(r2Scores)


plt.scatter(r2Scores, wInNum)
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.show()




###### OLD WORKING VERSION
#
# # Define the conditional model
# def weightedNormMTSIG(C_terms, L, w, sigma, b, condition):
#     result = np.zeros_like(C_terms)  # Initialize result array of same shape as C_terms
#     for i, C in enumerate(C_terms):
#         if condition[i]:  # Check condition for each C term
#             result[i] = (C * L * w**1) / ((C * w) + sigma) + b
#         else:
#             result[i] = (C * L) / (C + sigma) + b
#     return result
#
#
# # Example of the fitting procedure using curve_fit
# def fit_data(C_terms, y_data, condition):
#     # Define the objective function for curve fitting
#     def objective(C_terms, L, w, sigma, b):
#         return weightedNormMTSIG(C_terms, L, w, sigma, b, condition)
#
#     # Initial guess for L, w, and sigma
#     initial_guess = [.5, 0.5, 0.05, resp[0]]
#
#     # Perform the curve fitting
#     popt, pcov = curve_fit(objective, C_terms, y_data, p0=initial_guess)
#
#     return popt, pcov
#
#
# resp = np.concatenate((nonprefCenterMean, nonprefPeriMean), axis=0)
# cont = np.concatenate((contrasts, contrasts), axis=0)
# conditionArray = np.concatenate((np.full(7, False), np.full(7, True)), axis=0)
#
# popt, pcov = fit_data(cont, resp, conditionArray)
#
#
# # plotting
# xFit = np.logspace(-1, 2, 100)
#
# yFit = weightedNormMTSIG(xFit, *popt, np.full(100, False))
# plt.plot(xFit, yFit, color='green', linestyle='dashed')
# # yFit = weightedNormMTSIG(cont[:7], *popt, conditionArray[:7])
# # plt.plot(cont[:7], yFit, color='green', linestyle='dashed')
# plt.plot(contrasts, prefCenterMean, color='green')
#
# y_min, y_max = resp[0], np.max(yFit)
# y_halfway = (y_min + y_max) / 2
# # closest_index = np.argmin(np.abs(yFit - y_halfway))
# closest_index = np.argmin(np.abs(yFit - 0.5))
# x_halfway = xFit[closest_index]
# # plt.axhline(y=y_halfway, color='g')
# # plt.axvline(x=x_halfway, color='g', label=f'c50 = {x_halfway:.2f}')
# # plt.plot([xFit[0], x_halfway], [y_halfway, y_halfway], color='g')
# # plt.plot([x_halfway, x_halfway], [0, y_halfway], color='g', label=f'c50 = {x_halfway:.2f}')
# plt.plot([xFit[0], x_halfway], [0.5, 0.5], color='g')
# plt.plot([x_halfway, x_halfway], [0, 0.5], color='g', label=f'c50 = {x_halfway:.2f}')
#
#
# yFit = weightedNormMTSIG(xFit, *popt, np.full(100, True))
# plt.plot(xFit, yFit, color='green', alpha=0.5, linestyle='dashed')
# # yFit = weightedNormMTSIG(cont[7:], *popt, conditionArray[7:])
# # plt.plot(cont[7:], yFit, color='green', alpha=0.5, linestyle='dashed')
# plt.plot(contrasts, prefPeriMean, color='green', alpha=0.5)
#
# y_min, y_max = resp[0], np.max(yFit)
# y_halfway = (y_min + y_max) / 2
# # closest_index = np.argmin(np.abs(yFit - y_halfway))
# closest_index = np.argmin(np.abs(yFit - 0.5))
# x_halfway = xFit[closest_index]
# # plt.axhline(y=y_halfway, color='g', alpha=0.5, linestyle='-')
# # plt.axvline(x=x_halfway, color='g', alpha=0.5, linestyle='-', label=f'c50 = {x_halfway:.2f}')
# # plt.plot([xFit[0], x_halfway], [y_halfway, y_halfway], color='g', alpha=0.5)
# # plt.plot([x_halfway, x_halfway], [0, y_halfway], color='g', alpha=0.5, label=f'c50 = {x_halfway:.2f}')
# plt.plot([xFit[0], x_halfway], [0.5, 0.5], color='g', alpha=0.5)
# plt.plot([x_halfway, x_halfway], [0, 0.5], color='g', alpha=0.5, label=f'c50 = {x_halfway:.2f}')
#
# # plt.xscale('symlog', linthresh=0.1)
# plt.xlim(left=0)
# plt.ylim(bottom=0)
# plt.legend(loc='lower right')
# plt.show()





