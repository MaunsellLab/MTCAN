"""
MTSigma main analysis script. This script will extract spike counts during
valid trials from units active in a session and put them in the appropriate
matrix. From there, this script plots the contrast response function for pref/non-pref
at two points in the RF (center vs peri).

Chery - May 2024
"""

# Imports
from usefulFns import *
from MTSIGFunctions import *
import numpy as np
import psignifit as ps
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
import warnings

# # # # # 240903 has two files (MTSIG and MTSIG_unit2) with two units 3 and 8
# fileList = ['Akshan_240603', 'Akshan_240606', 'Akshan_240610', 'Akshan_240701',
#             'Akshan_240903', 'Akshan_240904', 'Akshan_240905', 'Akshan_240906',
#             'Akshan_240909', 'Akshan_240911', 'Akshan_240913', 'Akshan_240916',
#             'Akshan_240917', 'Akshan_240922', 'Akshan_240924', 'Akshan_241029',
#             'Akshan_241104']

fileList = ['Meetz_241115']

plotSingleDay = 0
if len(fileList) == 1:
    plotSingleDay = 1

potentialGoodUnits = []
potentialGoodUnitsCI = []

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

    # create folder and change directory to save PDFs and np.array
    if not os.path.exists('CRFgoodUnits1'):
        os.makedirs('CRFgoodUnits1')
    os.chdir('CRFgoodUnits1/')

    # list of indices of correctTrials (non-instruct, valid trialCertify)
    corrTrials = correctTrialsMTX(allTrials)

    # generate list of unique active units, and their channel
    units = activeUnits('spikeData', allTrials)
    unitCluster = allTrials[corrTrials[0]]['spikeTempInfo']['cgs']
    unitsChannel = unitsInfo(units, corrTrials, allTrials)

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

    # plot CRF
    unit = 42

    # plot single day all units
    if plotSingleDay == 1:
        for unit in units:
            t1 = time.time()
            unitID = np.where(units == unit)[0][0]
            baselineResp = np.mean(sponRate[unitID][:numContrasts * 4 * blocksDone]) * 1000 / trueStimDurMS

            colorID = 0
            plotColors = ['green', 'red', 'green', 'red']
            plotAlphas = [1, 1, 0.5, 0.5]
            plt.figure(figsize=(12, 8))
            warnings.simplefilter('ignore', OptimizeWarning)
            for i in range(2):
                for j in range(2):
                    response = meanSpikeReshaped[unitID][i][j].copy()
                    response = np.insert(response, 0, baselineResp)
                    try:
                        plt.scatter(contrasts, response, color=plotColors[colorID],
                                    alpha=plotAlphas[colorID])
                        initialGuess = [baselineResp, max(response), np.median(contrasts), 2.0]
                        pOpt, pCov = curve_fit(contrastFn, contrasts, response,
                                               bounds=([baselineResp, 0, 0, 0],
                                                       [np.inf, np.inf, np.inf, np.inf]))
                        xFit = np.logspace(-1, 2, 100)
                        yFit = contrastFn(xFit, *pOpt)
                        lower, upper = confidenceIntervalCRF(pOpt, pCov, xFit)
                        plt.plot(xFit, yFit, color=plotColors[colorID],
                                 alpha=plotAlphas[colorID], label=f'{pOpt[2]:.2f}')
                        if j == 2:
                            plt.fill_between(xFit, lower, upper, color=plotColors[colorID], alpha=0.2)
                    except (RuntimeError, ValueError) as e:
                        plt.scatter(contrasts, response, color=plotColors[colorID], alpha=plotAlphas[colorID])
                    colorID += 1

            plt.xscale('symlog', linthresh=0.1)
            plt.xlabel('Contrast (%)')
            plt.ylabel('Spikes/s')
            plt.title(f'Contrast Response Function, unit {unit}')
            plt.legend()
            sns.despine(offset=5)

            # save figs of good units
            plt.savefig(f'{unit}.pdf')
            plt.close('all')
    else:
        for unit in units:
            t1 = time.time()
            unitID = np.where(units == unit)[0][0]
            baselineResp = np.mean(sponRate[unitID][:numContrasts*4*blocksDone]) * 1000/trueStimDurMS

            # # check to see if unit meets inclusion criterion, if so add to good units
            # condition = True
            #
            # bonfCorrectedPVal = 0.05/6
            #
            # # check if 100% center pref resp is significantly greater than baseline
            # prefCenterSpikeCounts = spikeCountMat[unitID][:blocksDone, 5]
            # sponSpikeCounts = sponRate[unitID][:numContrasts*4*blocksDone]
            # # Perform the one-sided Mann-Whitney U test (pref center > baseline)
            # stat, p_value = mannwhitneyu(prefCenterSpikeCounts,
            #                              sponSpikeCounts, alternative='greater')
            # if p_value > bonfCorrectedPVal:
            #     condition = False
            #
            # # check if 100% peri pref resp is significantly greater than baseline
            # prefPeriSpikeCounts = spikeCountMat[unitID][:blocksDone, 17]
            # sponSpikeCounts = sponRate[unitID][:numContrasts*4*blocksDone]
            # # Perform the one-sided Mann-Whitney U test (pref center > baseline)
            # stat, p_value = mannwhitneyu(prefPeriSpikeCounts,
            #                              sponSpikeCounts, alternative='greater')
            # if p_value > bonfCorrectedPVal:
            #     condition = False
            #
            # # check if 100% center nonpref resp is significantly greater than baseline
            # nonprefCenterSpikeCounts = spikeCountMat[unitID][:blocksDone, 11]
            # stat, p_value = mannwhitneyu(nonprefCenterSpikeCounts,
            #                              sponSpikeCounts, alternative='greater')
            # if p_value > bonfCorrectedPVal:
            #     condition = False
            #
            # # check if 100% peri nonpref resp is significantly greater than baseline
            # nonprefPeriSpikeCounts = spikeCountMat[unitID][:blocksDone, 23]
            # stat, p_value = mannwhitneyu(nonprefPeriSpikeCounts,
            #                              sponSpikeCounts, alternative='greater')
            # if p_value > bonfCorrectedPVal:
            #     condition = False
            #
            # # check if 100% center pref significantly greater than 100% peri pref
            # stat, p_value = mannwhitneyu(prefCenterSpikeCounts,
            #                              prefPeriSpikeCounts, alternative='greater')
            # if p_value > bonfCorrectedPVal:
            #     condition = False
            #
            # # check if 100% center nonpref significantly greater than 100% peri nonpref
            # stat, p_value = mannwhitneyu(nonprefCenterSpikeCounts,
            #                              nonprefPeriSpikeCounts, alternative='greater')
            # if p_value > bonfCorrectedPVal:
            #     condition = False
            #
            # if condition:
            #     print(f'{seshDate}, {unit}')

            condition = True

            bonfCorrectedPVal = 0.05/6

            # check if at least 2 center pref resp is significantly greater than baseline
            prefCenterSpikeCounts = spikeCountMat[unitID][:blocksDone, :6]
            sponSpikeCounts = sponRate[unitID][:numContrasts*4*blocksDone]
            # Perform the one-sided Mann-Whitney U test (pref center > baseline)
            stat, p_value = mannwhitneyu(prefCenterSpikeCounts,
                                         sponSpikeCounts, alternative='greater')

            if np.sum([p_value < bonfCorrectedPVal]) < 2:
                condition = False

            # check if 100% peri pref resp is significantly greater than baseline
            prefPeriSpikeCounts = spikeCountMat[unitID][:blocksDone, 12:18]
            # Perform the one-sided Mann-Whitney U test (pref center > baseline)
            stat, p_value = mannwhitneyu(prefPeriSpikeCounts,
                                         sponSpikeCounts, alternative='greater')
            if np.sum([p_value < bonfCorrectedPVal]) < 2:
                condition = False

            # check if 100% center nonpref resp is significantly greater than baseline
            nonprefCenterSpikeCounts = spikeCountMat[unitID][:blocksDone, 6:12]
            stat, p_value = mannwhitneyu(nonprefCenterSpikeCounts,
                                         sponSpikeCounts, alternative='greater')
            if np.sum([p_value < bonfCorrectedPVal]) < 2:
                condition = False

            # check if 100% peri nonpref resp is significantly greater than baseline
            nonprefPeriSpikeCounts = spikeCountMat[unitID][:blocksDone, 18:]
            stat, p_value = mannwhitneyu(nonprefPeriSpikeCounts,
                                         sponSpikeCounts, alternative='greater')
            if np.sum([p_value < bonfCorrectedPVal]) < 2:
                condition = False

            print(unit, condition)

            numFailures = 0
            colorID = 0
            plotColors = ['green', 'red', 'green', 'red']
            plotAlphas = [1, 1, 0.5, 0.5]
            plt.figure(figsize=(12, 8))
            warnings.simplefilter('ignore', OptimizeWarning)
            for i in range(2):
                for j in range(2):
                    response = meanSpikeReshaped[unitID][i][j].copy()
                    response = np.insert(response, 0, baselineResp)
                    try:
                        plt.scatter(contrasts, response, color=plotColors[colorID],
                                    alpha=plotAlphas[colorID])
                        initialGuess = [baselineResp, max(response), np.median(contrasts), 2.0]
                        pOpt, pCov = curve_fit(contrastFn, contrasts, response,
                                               bounds=([baselineResp, 0, 0, 0],
                                                       [np.inf, np.inf, np.inf, np.inf]))
                        xFit = np.logspace(-1, 2, 100)
                        yFit = contrastFn(xFit, *pOpt)
                        lower, upper = confidenceIntervalCRF(pOpt, pCov, xFit)
                        plt.plot(xFit, yFit, color=plotColors[colorID],
                                 alpha=plotAlphas[colorID], label=f'{pOpt[2]:.2f}')
                        if j == 2:
                            plt.fill_between(xFit, lower, upper, color=plotColors[colorID], alpha=0.2)
                    except (RuntimeError, ValueError) as e:
                        plt.scatter(contrasts, response, color=plotColors[colorID], alpha=plotAlphas[colorID])
                        numFailures += 1
                    colorID += 1

            plt.xscale('symlog', linthresh=0.1)
            plt.xlabel('Contrast (%)')
            plt.ylabel('Spikes/s')
            plt.title(f'Contrast Response Function, unit {unit}')
            plt.legend()
            sns.despine(offset=5)

            if numFailures == 0:
                t1 = time.time()
                potentialGoodUnits.append(f'{seshDate}_{unit}')
                # bootstrap to estimate CI on c50
                prefCentRawSpikeCounts = spikeCountMat[unitID][:blocksDone, :6] * 1000 / trueStimDurMS
                prefCentRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
                                                            prefCentRawSpikeCounts))
                prefCentResp = prefCentRawSpikeCountsWithBase - baselineResp

                nonprefCentRawSpikeCounts = spikeCountMat[unitID][:blocksDone, 6:12] * 1000 / trueStimDurMS
                nonprefCentRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
                                                               nonprefCentRawSpikeCounts))
                nonprefCentResp = nonprefCentRawSpikeCountsWithBase - baselineResp

                prefPeriRawSpikeCounts = spikeCountMat[unitID][:blocksDone, 12:18] * 1000 / trueStimDurMS
                prefPeriRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
                                                            prefPeriRawSpikeCounts))
                prefPeriResp = prefPeriRawSpikeCountsWithBase - baselineResp

                nonprefPeriRawSpikeCounts = spikeCountMat[unitID][:blocksDone, 18:] * 1000 / trueStimDurMS
                nonprefPeriRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
                                                               nonprefPeriRawSpikeCounts))
                nonprefPeriResp = nonprefPeriRawSpikeCountsWithBase - baselineResp

                n_bootstraps = 500

                # Compute bootstrap c50s for each condition using vector approach
                c50bootstrapPrefCent = bootstrap_c50(prefCentResp, contrasts, n_bootstraps)
                c50bootstrapNonprefCent = bootstrap_c50(nonprefCentResp, contrasts, n_bootstraps)
                c50bootstrapPrefPeri = bootstrap_c50(prefPeriResp, contrasts, n_bootstraps)
                c50bootstrapNonprefPeri = bootstrap_c50(nonprefPeriResp, contrasts, n_bootstraps)

                print(f'done bootstrap {unit}')
                print(time.time() - t1)

                # Calculate confidence interval (e.g., 95%)
                # confidence_interval = np.percentile(c50bootstrapPrefCent, [2.5, 97.5])

                prefCentCI = np.percentile(c50bootstrapPrefCent, [2.5, 97.5])
                nonprefCentCI = np.percentile(c50bootstrapNonprefCent, [2.5, 97.5])
                prefPeriCI = np.percentile(c50bootstrapPrefPeri, [2.5, 97.5])
                nonprefPeriCI = np.percentile(c50bootstrapNonprefPeri, [2.5, 97.5])

                print(prefCentCI)
                print(nonprefCentCI)
                print(prefPeriCI)
                print(nonprefPeriCI)



                potentialGoodUnitsCI.append([prefCentCI, nonprefCentCI, prefPeriCI, nonprefPeriCI])

                randVar = 0
                threshValue = 30
                if prefCentCI[1] / prefCentCI[0] > threshValue:
                    randVar = 1
                if nonprefCentCI[1] / nonprefCentCI[0] > threshValue:
                    randVar = 1
                if prefPeriCI[1] / prefPeriCI[0] > threshValue:
                    randVar = 1
                if nonprefPeriCI[1] / nonprefPeriCI[0] > threshValue:
                    randVar = 1

                if randVar == 0:
                    print(f'unit meets threshold cutoff')

                # save figs of good units
                plt.savefig(f'{unit}.pdf')
                filePath = f'/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/Akshan/MTSIGIndividualNeuronPlots/{seshDate}_{unit}.pdf'
                plt.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
                plt.close('all')

    # to open another file in the loop
    os.chdir('../../../Python Code')


#################################################### EXTRA CODE ########################################################



import numpy as np
import matplotlib.pyplot as plt

# Define constants
L = 2  # Fixed value for L
sigma_values = [0.5]  # Different values of sigma to overlay
colorCode = ['blue', 'green', 'red']
weights = [1, 0.5, 0.25]

# Define a range of C values (contrast) on a logarithmic scale
C = np.logspace(-2, 2, 100)  # From 0.01 to 100 with 100 points

# Plotting
plt.figure(figsize=(8, 6))

# Loop over each sigma and calculate y values, then plot
for i, w in enumerate(weights):
    y = (C * w * L) / ((C * w) + sigma_values)
    plt.plot(C, y, label=w, color=colorCode[i])

# Log scale on x-axis
plt.xscale('log')

# Add labels and title
plt.xlabel("Contrast (C)")
plt.ylabel("Response (y)")
plt.title("Plot of (C * L) / (C + sigma) with Logarithmic X-axis for Different Sigma Values")
plt.legend(title="Sigma Values")
plt.grid(True, which="both", ls="--", lw=0.5)

plt.show()


##############

potentialGoodUnits = ['240603_56',
                         '240603_127',
                         '240603_128',
                         '240606_84',
                         '240606_85',
                         '240606_95',
                         '240606_104',
                         '240610_2',
                         '240610_11',
                         '240610_20',
                         '240610_22',
                         '240610_31',
                         '240610_33',
                         '240610_37',
                         '240610_50',
                         '240610_53',
                         '240610_58',
                         '240610_63',
                         '240610_68',
                         '240610_76',
                         '240610_90',
                         '240610_92',
                         '240610_97',
                         '240610_98',
                         '240610_103',
                         '240701_5',
                         '240903_3',
                         '240903_10',
                         '240903_11',
                         '240903_12',
                         '240903_13',
                         '240903_15',
                         '240904_2',
                         '240904_3',
                         '240904_8',
                         '240904_11',
                         '240904_15',
                         '240905_4',
                         '240906_7',
                         '240906_25',
                         '240906_26',
                         '240909_1',
                         '240909_4',
                         '240909_5',
                         '240909_7',
                         '240909_10',
                         '240909_12',
                         '240909_20',
                         '240909_25',
                         '240911_3',
                         '240911_16',
                         '240911_17',
                         '240913_1',
                         '240913_24',
                         '240913_27',
                         '240913_36',
                         '240913_43',
                         '240913_47',
                         '240913_48',
                         '240913_49',
                         '240913_50',
                         '240916_1',
                         '240916_7',
                         '240916_9',
                         '240916_10',
                         '240916_13',
                         '240917_13',
                         '240917_25',
                         '240917_44',
                         '240917_48',
                         '240917_49',
                         '240922_0',
                         '240922_4',
                         '240922_6',
                         '240922_9',
                         '240924_0',
                         '240924_6',
                         '240924_8',
                         '240924_19',
                         '240924_21',
                         '240924_22',
                         '240924_30',
                         '240924_35',
                         '241029_4',
                         '241029_22',
                         '241029_25',
                         '241029_29',
                         '241029_33',
                         '241104_35',
                         '241104_40',
                         '241104_41',
                         '241104_43']

potentialGoodUnitsCI = [[np.array([5.91795938e-01, 1.55342980e+12]),
                         np.array([4.09921058e-03, 4.52325208e+00]),
                         np.array([3.35897415e-03, 1.60030635e+12]),
                         np.array([1.87920026e-01, 7.19390378e+11])],
                        [np.array([2.78745115, 670.20364292]),
                          np.array([1.87201043, 3.62325889]),
                          np.array([3.04311842, 27.93791012]),
                          np.array([2.33837532e+00, 9.68002107e+11])],
                         [np.array([12.6302946, 39.21439745]),
                          np.array([1.26687548e+01, 1.00910985e+08]),
                          np.array([2.62178319e+01, 1.63618663e+08]),
                          np.array([2.45169611e+01, 6.20303540e+11])],
                         [np.array([4.788714, 18.81480995]),
                          np.array([6.24776448, 557.58970907]),
                          np.array([2.0797236e+01, 1.0000000e+06]),
                          np.array([1.23774212e+01, 3.28849989e+11])],
                         [np.array([1.68482014e+00, 8.07558632e+11]),
                          np.array([0.23639236, 2.98050607]),
                          np.array([2.65272676e+00, 1.75573847e+12]),
                          np.array([3.79382023e-02, 1.28893188e+12])],
                         [np.array([4.74478379e+00, 7.88150179e+05]),
                          np.array([4.71506486e+00, 1.19456934e+12]),
                          np.array([8.15175320e+00, 9.69474948e+11]),
                          np.array([1.49730112e+00, 2.29284541e+12])],
                         [np.array([7.93340496, 2944.30207565]),
                          np.array([10.20138045, 58.64161311]),
                          np.array([5.05694408e+00, 1.29782212e+09]),
                          np.array([6.91166302, 70.523994])],
                         [np.array([8.32155940e-02, 1.85357038e+12]),
                          np.array([0.02378062, 0.35782338]),
                          np.array([7.52912479e-02, 4.98167463e+11]),
                          np.array([2.04162874e-03, 9.50917974e+11])],
                         [np.array([2.08534158e-02, 2.33680374e+12]),
                          np.array([8.95943669e-02, 2.43727127e+12]),
                          np.array([0.00145322, 0.24052079]),
                          np.array([0.01011567, 3.24521127])],
                         [np.array([0.00279884, 0.22794315]),
                          np.array([0.00481064, 3.36527434]),
                          np.array([1.68915958e-02, 3.21525222e+12]),
                          np.array([1.73973918e-02, 2.85306942e+12])],
                         [np.array([3.13443595e+00, 1.12628302e+12]),
                          np.array([4.34350436, 6.29955188]),
                          np.array([2.27398651e-03, 5.40260129e+12]),
                          np.array([4.45125064e-03, 5.81276420e+12])],
                         [np.array([1.69042686e+00, 1.67790836e+12]),
                          np.array([1.33545207e-03, 2.97195048e+12]),
                          np.array([2.60900667e+00, 3.61285581e+12]),
                          np.array([5.30139560e-03, 4.26833513e+11])],
                         [np.array([0.0290039, 3.09745061]),
                          np.array([0.05522573, 5.32340575]),
                          np.array([0.02168106, 2.62443653]),
                          np.array([7.6210144e-02, 2.9199053e+11])],
                         [np.array([5.99208263, 10.37713362]),
                          np.array([5.97628002, 12.30153137]),
                          np.array([6.88745037e+00, 4.13474599e+11]),
                          np.array([6.23958835e-05, 2.50890708e+12])],
                         [np.array([0.02343651, 3.26744534]),
                          np.array([0.01539418, 3.16841556]),
                          np.array([1.17798472e-02, 2.57419912e+12]),
                          np.array([8.71656450e-03, 1.95976581e+12])],
                         [np.array([4.22371964e-07, 4.71144388e+12]),
                          np.array([2.09617774e+00, 7.09759337e+11]),
                          np.array([0.00409775, 3.25126429]),
                          np.array([0.00798083, 3.27267792])],
                         [np.array([2.09029375, 6.43343787]),
                          np.array([3.28881698, 6.39021692]),
                          np.array([1.79830426e-03, 2.11275447e+12]),
                          np.array([3.40682257e-03, 6.91781003e+00])],
                         [np.array([1.67495346e-04, 1.26093428e+01]),
                          np.array([1.73991731e-03, 3.02280612e+00]),
                          np.array([2.98217435e-04, 2.44487364e+12]),
                          np.array([3.74285352e-05, 3.12848310e+00])],
                         [np.array([0.00140212, 0.06842182]),
                          np.array([7.07799337e-06, 1.35867226e+06]),
                          np.array([5.85512492e-07, 2.35524164e+12]),
                          np.array([3.59786082e-02, 1.65514153e+12])],
                         [np.array([2.39872659e-03, 3.00728461e+00]),
                          np.array([0.00488512, 3.30397497]),
                          np.array([0.00388046, 3.31183235]),
                          np.array([0.00270616, 0.65395989])],
                         [np.array([3.23848889e-02, 2.07737049e+12]),
                          np.array([3.56804241e-02, 2.80739312e+12]),
                          np.array([3.87463702e-02, 1.89878093e+12]),
                          np.array([2.21495928, 5.09658935])],
                         [np.array([6.43046624, 8.79756347]),
                          np.array([5.98866486, 12.09232662]),
                          np.array([7.05270387, 22.72543583]),
                          np.array([5.2737344, 7.25586952])],
                         [np.array([6.79504115e-03, 1.78195590e+12]),
                          np.array([0.00445342, 3.06390002]),
                          np.array([1.07399769e-02, 7.75363908e+11]),
                          np.array([2.41076386e+00, 2.72793659e+12])],
                         [np.array([9.23083881, 16.44096856]),
                          np.array([8.37996063, 16.37604913]),
                          np.array([1.01359913e+01, 5.12106543e+11]),
                          np.array([3.67418975e+00, 1.15339435e+12])],
                         [np.array([4.0206162e-03, 3.1904042e+12]),
                          np.array([1.73101612e-03, 1.87829022e+11]),
                          np.array([0.0037516, 3.27648854]),
                          np.array([6.74458128e-03, 1.10609827e+12])],
                         [np.array([6.26627434, 8.24558974]),
                          np.array([7.8293776, 9.11201609]),
                          np.array([16.70782458, 69.71868602]),
                          np.array([7.97278381, 74.49656454])],
                         [np.array([3.20464465, 3.8825118]),
                          np.array([4.14502073, 4.75665163]),
                          np.array([4.24120996, 6.54245311]),
                          np.array([3.69049407, 4.56149131])],
                         [np.array([2.78830604e-06, 2.03549618e+12]),
                          np.array([1.02204855e-10, 2.77512661e+12]),
                          np.array([1.31328491e-11, 1.75482137e+01]),
                          np.array([1.46604129e-06, 4.09424221e+12])],
                         [np.array([4.27574841e-03, 5.00456821e+12]),
                          np.array([7.43540394e-03, 1.23469707e+12]),
                          np.array([6.12664880e-03, 3.01706451e+12]),
                          np.array([8.14181335e-04, 2.89017446e+12])],
                         [np.array([1.19149578e-03, 2.17616265e+12]),
                          np.array([5.07151451e-06, 2.14906484e+12]),
                          np.array([8.91431470e-09, 9.55244088e+01]),
                          np.array([6.71585014e-15, 4.18652658e+12])],
                         [np.array([1.03453272e-03, 2.81217802e+00]),
                          np.array([0.06352958, 0.18608927]),
                          np.array([2.02159080e-02, 8.71230539e+11]),
                          np.array([0.05614381, 3.15500334])],
                         [np.array([1.11673308e-02, 3.02752227e+11]),
                          np.array([1.48571231e-07, 2.03156452e+12]),
                          np.array([2.66485574e-04, 2.96061667e+00]),
                          np.array([1.68767251e-02, 1.43853640e+12])],
                         [np.array([4.01444065e-03, 9.54261948e+11]),
                          np.array([5.88627068e-03, 3.94566062e+12]),
                          np.array([8.71076915e-04, 3.56964549e+12]),
                          np.array([3.36581050e-03, 1.35449228e+12])],
                         [np.array([7.16580784e-03, 3.26181557e+12]),
                          np.array([0.03580352, 3.10375443]),
                          np.array([0.02438818, 1.81455495]),
                          np.array([8.63530224e-03, 2.35433928e+12])],
                         [np.array([5.6456856, 7.67568749]),
                          np.array([5.41656029, 7.10325149]),
                          np.array([5.16133264, 7.11529888]),
                          np.array([7.13098112, 9.98507375])],
                         [np.array([0.0202552, 0.4174526]),
                          np.array([2.19593035e-02, 1.78502342e+12]),
                          np.array([0.03503363, 3.02951741]),
                          np.array([1.88679805e-02, 1.07614191e+12])],
                         [np.array([1.26579519e-14, 3.24773702e+12]),
                          np.array([3.59922711e-04, 4.24435553e+00]),
                          np.array([1.08369464e-07, 9.49159986e+00]),
                          np.array([1.74633947e-07, 2.11561114e+02])],
                         [np.array([3.80742058, 5.463251]),
                          np.array([7.03782658, 16.63562293]),
                          np.array([6.54249461, 10.8866857]),
                          np.array([1.20157934e+01, 2.78907481e+08])],
                         [np.array([0.06989467, 3.00999127]),
                          np.array([2.3095101, 3.50501535]),
                          np.array([7.15400306e-03, 1.01469652e+12]),
                          np.array([4.83459153e-02, 2.08219640e+12])],
                         [np.array([2.62081461, 4.24510016]),
                          np.array([3.63126444, 4.77316951]),
                          np.array([4.44870825, 35.55160829]),
                          np.array([3.91036154, 48.85116979])],
                         [np.array([0.02486275, 3.03337505]),
                          np.array([0.07151637, 2.93684833]),
                          np.array([0.3250701, 3.03756484]),
                          np.array([0.01153511, 2.88472509])],
                         [np.array([0.02650174, 3.13403384]),
                          np.array([0.0215162, 3.0865668]),
                          np.array([0.01286833, 2.64257412]),
                          np.array([4.40959847e-02, 2.82561710e+12])],
                         [np.array([0.01087976, 0.22366487]),
                          np.array([1.35187151e-02, 7.69942649e+11]),
                          np.array([4.41465344e-02, 2.54450513e+12]),
                          np.array([1.18748005e-02, 1.78383389e+12])],
                         [np.array([2.061746, 5.98215278]),
                          np.array([2.00632073e-03, 3.37591712e+00]),
                          np.array([1.85217585e+00, 1.68988441e+12]),
                          np.array([1.53876470e-02, 2.31982896e+12])],
                         [np.array([4.6606173, 13.269681]),
                          np.array([1.18916649e+01, 7.49399614e+11]),
                          np.array([3.67912313e-03, 4.41170816e+00]),
                          np.array([5.96601200e-02, 1.10621551e+12])],
                         [np.array([4.37720966, 5.47344098]),
                          np.array([5.89297968, 73.85459234]),
                          np.array([5.10084485, 5.79850179]),
                          np.array([1.32904961e+01, 9.53206196e+09])],
                         [np.array([4.70824116e-01, 3.05109833e+12]),
                          np.array([6.28737654e-03, 4.97371225e+01]),
                          np.array([5.34078912e-03, 3.24880275e+12]),
                          np.array([9.08341458e-04, 2.01150877e+12])],
                         [np.array([0.03018848, 3.1115513]),
                          np.array([1.43514746e+00, 1.67352590e+12]),
                          np.array([5.15028133e-02, 6.00415924e+11]),
                          np.array([0.04177843, 8.959634])],
                         [np.array([0.01820455, 3.01789195]),
                          np.array([0.00655359, 0.25020116]),
                          np.array([1.32684192e+00, 1.02875983e+12]),
                          np.array([1.73460836e+00, 5.03122158e+11])],
                         [np.array([0.05284456, 3.15443766]),
                          np.array([1.87055393e-01, 2.20024773e+12]),
                          np.array([1.37506560e-02, 2.13425881e+12]),
                          np.array([0.01658045, 3.00280825])],
                         [np.array([1.65511347e+00, 3.26329708e+03]),
                          np.array([2.96991620e-02, 1.24642628e+12]),
                          np.array([4.64265344e-02, 7.81585094e+11]),
                          np.array([0.04674361, 3.17127755])],
                         [np.array([4.85323289, 5.5512667]),
                          np.array([6.84606725, 8.09161383]),
                          np.array([5.93558101, 7.50041515]),
                          np.array([5.62720582, 7.34175476])],
                         [np.array([4.06748305, 6.61002687]),
                          np.array([3.93662846e+00, 1.25620130e+12]),
                          np.array([7.19608673, 65.01038031]),
                          np.array([2.08477623, 11.31697458])],
                         [np.array([4.04186723, 13.07970239]),
                          np.array([2.59272955, 32.03193838]),
                          np.array([6.12259850e+00, 7.31893435e+07]),
                          np.array([5.33194388e+00, 1.00000000e+06])],
                         [np.array([10.28439443, 16.07475503]),
                          np.array([9.26701179, 35.94803934]),
                          np.array([12.82183585, 24.82876342]),
                          np.array([1.47535364e+01, 8.01151424e+07])],
                         [np.array([5.71029629, 13.54966215]),
                          np.array([3.20104281e+00, 1.28636457e+12]),
                          np.array([7.90779360e+00, 6.32138122e+11]),
                          np.array([4.41075632e+00, 8.56692274e+11])],
                         [np.array([8.32527659, 155.09028454]),
                          np.array([3.73480702e+00, 1.04467708e+04]),
                          np.array([6.20987728, 19.99710408]),
                          np.array([4.37721422e+00, 1.01348590e+10])],
                         [np.array([7.1210537, 11.11066684]),
                          np.array([6.07662699, 51.12200114]),
                          np.array([13.28821151, 24.73820064]),
                          np.array([9.29498734, 38.87128898])],
                         [np.array([4.77941611, 13.44166733]),
                          np.array([5.47737677, 11.96866407]),
                          np.array([1.29574290e+01, 5.34986851e+05]),
                          np.array([8.88321323, 3858.66484113])],
                         [np.array([1.54054029, 3.05985391]),
                          np.array([0.03175636, 0.28081744]),
                          np.array([6.81118216e-02, 1.81293770e+12]),
                          np.array([0.08983594, 0.28157716])],
                         [np.array([5.99825409e-03, 6.79043044e+00]),
                          np.array([2.89254468, 4.82569973]),
                          np.array([6.56905826, 20.63131318]),
                          np.array([3.24339616e+00, 6.34568276e+11])],
                         [np.array([1.2591024e+01, 7.4554260e+10]),
                          np.array([8.20542711, 16.74970851]),
                          np.array([1.1902828e+01, 8.2831720e+08]),
                          np.array([1.17080546e+01, 5.35198964e+09])],
                         [np.array([6.7817362, 8.87746687]),
                          np.array([7.83155625, 19.39980439]),
                          np.array([15.57137548, 23.81934829]),
                          np.array([14.77879726, 69.1990983])],
                         [np.array([13.08530377, 79.97717657]),
                          np.array([20.28179312, 88.69006887]),
                          np.array([2.59180456e+01, 5.44426888e+09]),
                          np.array([2.24706110e+01, 1.83904941e+09])],
                         [np.array([4.45685875, 25.39082914]),
                          np.array([3.29754692, 5.11031137]),
                          np.array([8.67837668e+00, 2.17952741e+12]),
                          np.array([6.77619347e+00, 7.53545726e+07])],
                         [np.array([3.52607178, 5.38935422]),
                          np.array([5.95702303, 101.80700512]),
                          np.array([4.02857728, 6.01561493]),
                          np.array([6.37982182, 961.49499859])],
                         [np.array([8.53994154e-03, 6.35563746e+11]),
                          np.array([1.50038638e-02, 1.74441998e+12]),
                          np.array([8.10015118e-03, 3.23633008e+12]),
                          np.array([2.78155063e-02, 5.55956007e+11])],
                         [np.array([0.01552931, 2.8445541]),
                          np.array([0.02823941, 3.07820388]),
                          np.array([2.04185155e-02, 6.05313799e+11]),
                          np.array([2.07187634e-02, 1.08494219e+12])],
                         [np.array([3.8427751, 4.34096662]),
                          np.array([4.76179933, 5.78289559]),
                          np.array([7.19170145, 12.52839962]),
                          np.array([10.27655839, 30.62224458])],
                         [np.array([1.03122369e-02, 3.02982533e+12]),
                          np.array([2.87078873e+00, 2.44925985e+12]),
                          np.array([7.42309392e-02, 3.97478960e+12]),
                          np.array([0.01739307, 3.07672])],
                         [np.array([2.34508616e-03, 6.43300132e+12]),
                          np.array([2.58615213e-03, 4.62220542e+00]),
                          np.array([4.83213925e-03, 4.39829607e+12]),
                          np.array([3.72800234e-03, 4.88448124e+00])],
                         [np.array([0.00223493, 0.19975973]),
                          np.array([4.93844266e-01, 2.14485390e+11]),
                          np.array([3.62050461e-03, 3.16829677e+12]),
                          np.array([6.04165078e-02, 4.51357188e+12])],
                         [np.array([0.00442602, 0.22763587]),
                          np.array([0.02695075, 0.29076182]),
                          np.array([3.16028193e-02, 6.20141572e+02]),
                          np.array([6.37265034e-01, 1.86017823e+12])],
                         [np.array([3.69580262, 9.99542773]),
                          np.array([3.75032234, 10.15160575]),
                          np.array([6.28439126, 22.48722499]),
                          np.array([8.329369, 21.56686908])],
                         [np.array([0.00038162, 0.04881938]),
                          np.array([1.13218259e-05, 4.62010546e+00]),
                          np.array([4.60679100e-07, 4.51853102e+12]),
                          np.array([4.72094178e-04, 1.00000000e+06])],
                         [np.array([0.07094972, 2.97606792]),
                          np.array([0.014052, 2.5995031]),
                          np.array([1.28782677e-02, 3.18536334e+11]),
                          np.array([1.48434533e-02, 2.50215179e+11])],
                         [np.array([2.27561162e-03, 3.10225124e+00]),
                          np.array([1.31699340e-03, 6.34780776e+12]),
                          np.array([1.05432326e-03, 2.97088229e+00]),
                          np.array([0.00644191, 3.52707766])],
                         [np.array([0.01565708, 3.97994984]),
                          np.array([1.47250040e-02, 1.95280462e+12]),
                          np.array([4.63867529e-03, 2.96723862e+12]),
                          np.array([3.17635583e+00, 2.14561440e+12])],
                         [np.array([3.58677549e-03, 1.81148948e+12]),
                          np.array([2.82945265e-02, 3.59457645e+12]),
                          np.array([2.10941873e-02, 5.38403754e+12]),
                          np.array([3.91260249e-03, 4.22660281e+12])],
                         [np.array([3.29849552, 4.55612954]),
                          np.array([4.39164114, 5.34840981]),
                          np.array([5.25317058, 8.40694956]),
                          np.array([4.98133503, 5.77675179])],
                         [np.array([6.68930434, 8.83867994]),
                          np.array([12.66830781, 31.48472501]),
                          np.array([14.08143924, 28.59858841]),
                          np.array([2.64368018e+01, 8.63666686e+07])],
                         [np.array([6.14775256e-03, 6.96044180e+11]),
                          np.array([1.03125856e-02, 2.14997859e+12]),
                          np.array([4.47457503e-01, 1.89014638e+12]),
                          np.array([6.34438913e-03, 1.18826050e+12])],
                         [np.array([5.94035604, 8.11616528]),
                          np.array([5.67625661, 8.61942032]),
                          np.array([16.12110021, 181.59529554]),
                          np.array([8.04909333, 61.25041689])],
                         [np.array([2.43958180e+01, 1.54312877e+12]),
                          np.array([2.79648602e-03, 3.72092690e+12]),
                          np.array([0.00464098, 3.4116716]),
                          np.array([5.04353896e-03, 1.33142060e+12])],
                         [np.array([4.32690857e+00, 7.12767815e+11]),
                          np.array([4.89292975, 30.96000801]),
                          np.array([3.48284018e+00, 5.73844587e+10]),
                          np.array([0.00512019, 4.03060692])],
                         [np.array([7.03659768, 22.83684374]),
                          np.array([9.15997175e+00, 1.24968534e+10]),
                          np.array([5.39374591e+00, 3.61108598e+11]),
                          np.array([5.04805910e+00, 1.89919632e+12])],
                         [np.array([1.30867253e+01, 1.38216997e+05]),
                          np.array([1.57130158e+01, 3.10200253e+10]),
                          np.array([4.68906665e+00, 6.94503559e+11]),
                          np.array([4.89217530e+00, 5.31395984e+10])],
                         [np.array([9.69017717, 26.44982746]),
                          np.array([1.59932115e+01, 2.29728877e+11]),
                          np.array([3.25829152, 11.48692901]),
                          np.array([3.10016499e+00, 3.37985481e+11])],
                         [np.array([12.4115855, 27.33485503]),
                          np.array([8.53069983, 16.07253588]),
                          np.array([7.92285466, 145.46867793]),
                          np.array([2.31747505e+01, 8.48260937e+08])],
                         [np.array([5.25900571, 7.23170626]),
                          np.array([6.26891976, 13.59683722]),
                          np.array([4.31840052, 5.97800111]),
                          np.array([1.18328494e+01, 8.11253963e+11])],
                         [np.array([5.19294268, 7.06314957]),
                          np.array([10.80611961, 91.0104407]),
                          np.array([4.18695531, 5.99901771]),
                          np.array([3.2766317, 8.44872007])],
                         [np.array([5.07969014, 6.21345736]),
                          np.array([6.76787178, 11.49602165]),
                          np.array([4.42463312, 5.11127185]),
                          np.array([5.32242383, 12.93046615])]]

# Find max ratio for each sublist in a
max_ratios = []
for sublist in potentialGoodUnitsCI:
    ratios = [arr[1] / arr[0] for arr in sublist]  # Calculate ratios
    max_ratios.append(max(ratios))  # Get max ratio for each sublist

# Convert max_ratios to a numpy array if needed
max_ratios = np.array(max_ratios)
potentialGoodUnits = np.array(potentialGoodUnits)


potentialGoodUnits[np.where(max_ratios < 1000)[0]]
potentialGoodUnits[np.where(max_ratios < 50)[0]]


### Simulation of c50 estimation as a function of Rmax (different neuron rates)
# Naka-Rushton function without baseline
def naka_rushton(C, Rmax, c50, gamma=2):
    return (Rmax * C**gamma) / (C**gamma + c50**gamma)

# Simulate spike counts using Poisson process
def simulate_spikes(contrasts, Rmax, c50, gamma, num_trials=100):
    mean_rates = naka_rushton(contrasts, Rmax, c50, gamma)
    spike_counts = np.random.poisson(mean_rates[:, None], (len(contrasts), num_trials))
    return spike_counts

# Fit the Naka-Rushton function to mean spike rates
def fit_naka_rushton(contrasts, mean_rates):
    initial_guess = [np.max(mean_rates), 0.5]
    try:
        popt, _ = curve_fit(
            lambda C, Rmax, c50: naka_rushton(C, Rmax, c50),
            contrasts,
            mean_rates,
            p0=initial_guess,
            maxfev=10000  # Increase max function evaluations
        )
        return popt[1]  # Return only the estimated c50
    except RuntimeError:
        return np.nan  # Return NaN if the fitting fails

# Parameters
contrasts = np.logspace(np.log10(0.03), np.log10(1.0), 6)  # 6 contrasts from ~3% to 100%
true_c50 = 0.3  # True c50
gamma = 2       # Slope parameter
num_trials = 100  # Number of trials per contrast
Rmax_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100]  # Specific Rmax values to simulate
Rmax_values = np.logspace(0, 1.5, 20)
num_simulations = 5000  # Number of repetitions for each Rmax

# Store results
results = {}

for Rmax in Rmax_values:
    c50_estimates = []
    for _ in range(num_simulations):
        # Simulate spike counts
        spike_counts = simulate_spikes(contrasts, Rmax, true_c50, gamma, num_trials)
        # Calculate mean rates across trials
        mean_rates = np.mean(spike_counts, axis=1)
        # Fit the Naka-Rushton function and get estimated c50
        estimated_c50 = fit_naka_rushton(contrasts, mean_rates)
        c50_estimates.append(estimated_c50)
    # Filter out NaN values
    valid_c50_estimates = [c for c in c50_estimates if not np.isnan(c)]
    # Compute mean and SEM of estimated c50
    c50_mean = np.mean(valid_c50_estimates)
    c50_sem = np.std(valid_c50_estimates) / np.sqrt(len(valid_c50_estimates))
    results[Rmax] = (c50_mean, c50_sem)

# Plot results
fig, ax = plt.subplots(figsize=(8, 6))
Rmax_log = np.log10(Rmax_values)
means = [results[Rmax][0] for Rmax in Rmax_values]
sems = [results[Rmax][1] for Rmax in Rmax_values]

ax.errorbar(Rmax_log, means, yerr=sems, fmt='o', label='Estimated c50 ± SEM')
ax.axhline(y=true_c50, color='r', linestyle='--', label='True c50')
ax.set_xticks(Rmax_log)
ax.set_xticklabels([f"{Rmax:.1f}" for Rmax in Rmax_values])
ax.set_xlabel('Rmax (Hz, log scale)')
ax.set_ylabel('Estimated c50')
ax.set_title('Estimated c50 for Different Rmax Values')
ax.legend()
ax.grid(True)
plt.show()

