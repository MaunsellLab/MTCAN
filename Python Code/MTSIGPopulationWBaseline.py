"""
MTSigma population main analysis script. This script will extract spike counts during
valid trials from units active in a session and put them in the appropriate
matrix. Then we will add this to a population matrix to measure how c50 changes
with rMax at a population level.

This script differs rom MTSIGPopulation in that it doesnt subtract the baseline out. Helps with testing
normalization equivalence of Lnaught and sigma

Chery - November 2024

Notes: 240910 and 240911 are the same unit (just labelled different files to help with analysis)
       240912 and 240913 are the same units (240913 has better stimulus optimization to drive more units in an optimized way)
"""

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
            'Akshan_241104', 'Akshan_241112', 'Meetz_241115']

# # unit list with 50x difference b/w CI (for inclusion criteria) 241104, 240913, 240924 resorted
unitList = ['240610_92', '240701_5', '240903_3', '240904_8', '240906_25',
            '240911_17', '240913_0', '240913_4', '240913_13', '240913_44',
            '240916_7', '240917_44', '240922_6', '240924_25', '240924_37',
            '241104_3', '241104_4', '241104_23', '241112_18', '241112_19',
            '241115_2', '241115_3', '241115_15', '241115_18', '241115_32',
            '241115_34']

prefCenterPopResp = []
prefPeriPopResp = []
nonprefCenterPopResp = []
nonprefPeriPopResp = []
popBlocksDone = []
popUnitCluster = []
popGaborSD = []

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

    # RF Gabor SD
    popGaborSD.append(allTrials[0]['rfGabor']['data']['sigmaDeg'])

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
    if seshDate == '241115':
        indx = np.where(units == 18)[0][0]
        meanSpikeReshaped[indx] = meanSpikeReshaped[indx][:, ::-1, :]
        indx = np.where(units == 34)[0][0]
        meanSpikeReshaped[indx] = meanSpikeReshaped[indx][:, ::-1, :]

    # fit CRF
    for unit in sessionGoodUnits:
        count = np.where(units == unit)[0][0]
        popUnitCluster.append(unitCluster[count])
        baselineResp = np.mean(sponRate[count][:numContrasts * 4 * blocksDone]) * 1000 / trueStimDurMS
        warnings.simplefilter('ignore', OptimizeWarning)

        # create population avg CRF curve for pref resp at center and peri
        # response matrix is baseline subtracted and normalized to center resp
        # at 100% contrast
        prefCenterResp = meanSpikeReshaped[count][0][0]
        prefCenterResp = np.insert(prefCenterResp, 0, baselineResp)
        # prefCenterResp = prefCenterResp - baselineResp
        prefPeriResp = meanSpikeReshaped[count][1][0]
        prefPeriResp = np.insert(prefPeriResp, 0, baselineResp)
        # prefPeriResp = prefPeriResp - baselineResp
        prefCenterPopResp.append(prefCenterResp/prefCenterResp[-1])
        prefPeriPopResp.append(prefPeriResp/prefCenterResp[-1])

        # create population avg CRF curve for nonpref resp at center and peri
        # response matrix is baseline subtracted and normalized to center resp
        # at 100% contrast
        nonprefCenterResp = meanSpikeReshaped[count][0][1]
        nonprefCenterResp = np.insert(nonprefCenterResp, 0, baselineResp)
        # nonprefCenterResp = nonprefCenterResp - baselineResp
        nonprefPeriResp = meanSpikeReshaped[count][1][1]
        nonprefPeriResp = np.insert(nonprefPeriResp, 0, baselineResp)
        # nonprefPeriResp = nonprefPeriResp - baselineResp
        nonprefCenterPopResp.append(nonprefCenterResp/prefCenterResp[-1])
        nonprefPeriPopResp.append(nonprefPeriResp/prefCenterResp[-1])
        popBlocksDone.append(blocksDone)

    # to open another file in the loop
    os.chdir('../../Python Code')

print(time.time()-t0)

# population plots
prefCenterPopResp = np.array(prefCenterPopResp)
prefPeriPopResp = np.array(prefPeriPopResp)
nonprefCenterPopResp = np.array(nonprefCenterPopResp)
nonprefPeriPopResp = np.array(nonprefPeriPopResp)
popBlocksDone = np.array(popBlocksDone)
popUnitCluster = np.array(popUnitCluster)

prefCenterMean = np.mean(prefCenterPopResp, axis=0)
prefPeriMean = np.mean(prefPeriPopResp, axis=0)
prefCenterSEM = stats.sem(prefCenterPopResp, axis=0)
prefPeriSEM = stats.sem(prefPeriPopResp, axis=0)

nonprefCenterMean = np.mean(nonprefCenterPopResp, axis=0)
nonprefPeriMean = np.mean(nonprefPeriPopResp, axis=0)
nonprefCenterSEM = stats.sem(nonprefCenterPopResp, axis=0)
nonprefPeriSEM = stats.sem(nonprefPeriPopResp, axis=0)

hfont = {'fontname': 'Arial'}
lineWidth = 1
fontSize = 8
mSize = 4

# contrastFN params = r0, rMax, c50, n)

# plot pref and non pref center CRF using Naka-Rushton
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    x_min = ax1.get_xlim()[0]
    pOpt, pCov = curve_fit(contrastFn, contrasts, prefCenterMean,
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='black', label=f'c50={pOpt[2]:.2f}, rMax={pOpt[1]:.2f}')

    # Add truncated lines for pref curve
    c50_pref = pOpt[2]
    rMax_pref = pOpt[1]
    half_response_pref = 0.5 * rMax_pref
    ax1.plot([c50_pref, c50_pref], [0, half_response_pref],
             color='black', linestyle=':', label='Pref c50')
    ax1.plot([x_min, c50_pref], [half_response_pref, half_response_pref],
             color='black', linestyle=':')

    # plot nonpref Center
    pOpt, pCov = curve_fit(contrastFn, contrasts, nonprefCenterMean,
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='gray', label=f'c50={pOpt[2]:.2f}, rMax={pOpt[1]:.2f}')

    # Add truncated lines for non-pref curve
    c50_nonpref = pOpt[2]
    rMax_nonpref = pOpt[1]
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

    pOpt, pCov = curve_fit(contrastFn, contrasts, prefPeriMean,
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='black',
             label=f'c50={pOpt[2]:.2f}, rMax={pOpt[1]:.2f}')

    # Add truncated lines for pref curve
    c50_pref = pOpt[2]
    rMax_pref = pOpt[1]
    half_response_pref = 0.5 * rMax_pref
    ax1.plot([c50_pref, c50_pref], [0, half_response_pref],
             color='black', linestyle=':', label='Pref c50')
    ax1.plot([0, c50_pref], [half_response_pref, half_response_pref],
             color='black', linestyle=':')

    # plot nonpref peri
    pOpt, pCov = curve_fit(contrastFn, contrasts, nonprefPeriMean,
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='gray',
             label=f'c50={pOpt[2]:.2f}, rMax={pOpt[1]:.2f}')

    # Add truncated lines for non-pref curve
    c50_nonpref = pOpt[2]
    rMax_nonpref = pOpt[1]
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

# combined plot of pref/non pref at center and peri
for filler in range(1):
    fig, ax1 = plt.subplots(figsize=(3, 3))

    # plot pref center
    x_min = ax1.get_xlim()[0]
    pOpt, pCov = curve_fit(contrastFn, contrasts, prefCenterMean,
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='black', label=f'c50={pOpt[2]:.2f}', linewidth=lineWidth)

    # Add truncated lines for pref curve
    c50_pref = pOpt[2]
    rMax_pref = pOpt[1]
    half_response_pref = (0.5 * (rMax_pref - pOpt[0])) + (2*pOpt[0])
    # ax1.plot([c50_pref, c50_pref], [0, half_response_pref],
    #          color='black', linestyle=':')
    ax1.scatter(c50_pref, half_response_pref, marker="x", color='black', s=mSize**2)

    # plot nonpref Center
    pOpt, pCov = curve_fit(contrastFn, contrasts, nonprefCenterMean,
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='black', label=f'c50={pOpt[2]:.2f}', linewidth=lineWidth)

    # Add truncated lines for non-pref curve
    c50_nonpref = pOpt[2]
    rMax_nonpref = pOpt[1]
    half_response_nonpref = (0.5 * (rMax_nonpref - pOpt[0])) + (2*pOpt[0])
    # ax1.plot([c50_nonpref, c50_nonpref], [0, half_response_nonpref],
    #          color='black', linestyle=':')
    ax1.scatter(c50_nonpref, half_response_nonpref, marker="x", color='black', s=mSize**2)

    # plot pref peri
    pOpt, pCov = curve_fit(contrastFn, contrasts, prefPeriMean,
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='gray', linestyle='--',
             label=f'c50={pOpt[2]:.2f}', linewidth=lineWidth)

    # Add truncated lines for pref curve
    c50_pref = pOpt[2]
    rMax_pref = pOpt[1]
    half_response_pref = (0.5 * (rMax_pref - pOpt[0])) + (2*pOpt[0])
    # ax1.plot([c50_pref, c50_pref], [0, half_response_pref],
    #          color='gray', linestyle=':')
    ax1.scatter(c50_pref, half_response_pref, marker="x", color='gray',
                s=mSize**2)

    # plot nonpref peri
    initialGuess = [max(nonprefPeriMean), np.median(contrasts), 2.0]
    pOpt, pCov = curve_fit(contrastFn, contrasts, nonprefPeriMean,
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]))
    xFit = np.logspace(-1, 2, 100)
    yFit = contrastFn(xFit, *pOpt)
    ax1.plot(xFit, yFit, color='gray', linestyle='--',
             label=f'c50={pOpt[2]:.2f}', linewidth=lineWidth)

    # Add truncated lines for non-pref curve
    c50_nonpref = pOpt[2]
    rMax_nonpref = pOpt[1]
    half_response_nonpref = (0.5 * (rMax_nonpref - pOpt[0])) + (2*pOpt[0])
    # ax1.plot([c50_nonpref, c50_nonpref], [0, half_response_nonpref],
    #          color='gray', linestyle=':')
    ax1.scatter(c50_nonpref, half_response_nonpref, marker="x", color='gray',
                s=mSize**2)

    # error bar scatters
    ax1.errorbar(contrasts[1:], prefPeriMean[1:], yerr=prefPeriSEM[1:],
                 fmt='o', color='gray', markersize=mSize)
    ax1.errorbar(contrasts[1:], nonprefPeriMean[1:], yerr=nonprefPeriSEM[1:],
                 fmt='o', mfc='w', color='gray', markersize=mSize)
    ax1.errorbar(contrasts[1:], prefCenterMean[1:], yerr=prefCenterSEM[1:], fmt='o',
                 color='black', markersize=mSize)
    ax1.errorbar(contrasts[1:], nonprefCenterMean[1:], yerr=nonprefCenterSEM[1:], fmt='o',
                 mfc='w', color='black', markersize=mSize)

    ax1.set_xscale('symlog', linthresh=0.1)
    ax1.set_xlabel('Contrast (%)', **hfont, fontsize=fontSize)
    ax1.set_ylabel('Normalized\nRate', **hfont, fontsize=fontSize)
    ax1.yaxis.label.set_rotation(0)
    # ax1.set_title(f'Population Contrast Response Function n={len(prefCenterPopResp)}')

    ax1.set_ylim([0, 1.1])
    # Explicitly set y-ticks at 0, 1, and the others from get_yticks()
    y_ticks = ax1.get_yticks()  # Get all existing y-ticks
    if 1 not in y_ticks:
        y_ticks = list(y_ticks) + [1]  # Add 1 explicitly if missing
    # create labels: Only 0 and 1 get labels, others are empty
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1 else
        '' for tick in y_ticks
    ]
    # Apply the y-ticks and labels
    ax1.set_yticks(y_ticks)  # Set ticks (including 0 and 1)
    ax1.set_yticklabels(y_tick_labels, **hfont, fontsize=fontSize)  # Label only 0 and 1

    # Explicitly set x-ticks and labels
    x_ticks = [0.1, 1, 10, 100]  # Define tick positions
    x_tick_labels = [f"{tick:.1f}" if tick < 1 else f"{int(tick)}" for tick in x_ticks]  # Format labels
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_tick_labels, fontsize=fontSize, **hfont)
    ax1.tick_params(axis='both', width=1, length=4, direction='out', pad=3, labelsize=fontSize)

    ax1.legend()
    ax1.set_xlim(left=0.3)
    ax1.xaxis.set_minor_locator(mticker.LogLocator(numticks=99, subs="auto"))
    # sns.despine(offset=5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(1)
    ax1.spines['left'].set_linewidth(1)

# save fig
plt.savefig('/Users/chery/Documents/Grad School/Maunsell Lab/Manuscripts/PythonMainFigures/CRF_plot.pdf',
            bbox_inches='tight')
plt.close('all')


# plot box plots of center/peri pref/non-pref fit c50 naka-rushton
centerPrefPopC50 = []
centerPrefPopRMax = []
centerPrefPopN = []
centerNonprefPopC50 = []
centerNonprefPopRMax = []
centerNonprefPopN = []
periPrefPopC50 = []
periPrefPopRMax = []
periPrefPopN = []
periNonprefPopC50 = []
periNonprefPopRMax = []
periNonprefPopN = []
for i in range(len(prefCenterPopResp)):
    # pref center
    initialGuess = [max(prefCenterPopResp[i]), np.median(contrasts),
                    2.0, prefCenterPopResp[i][0]]
    pOpt, pCov = curve_fit(contrastFn, contrasts, prefCenterPopResp[i],
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]),
                           p0=initialGuess)
    centerPrefPopC50.append(pOpt[2])
    centerPrefPopRMax.append(pOpt[1])
    centerPrefPopN.append(pOpt[3])

    # nonpref center
    initialGuess = [max(nonprefCenterPopResp[i]), np.median(contrasts),
                    2.0, nonprefCenterPopResp[i][0]]
    pOpt, pCov = curve_fit(contrastFn, contrasts, nonprefCenterPopResp[i],
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]),
                           p0=initialGuess)
    centerNonprefPopC50.append(pOpt[2])
    centerNonprefPopRMax.append(pOpt[1])
    centerNonprefPopN.append(pOpt[3])

    # pref peri
    initialGuess = [max(prefPeriPopResp[i]), np.median(contrasts),
                    2.0, prefPeriPopResp[i][0]]
    pOpt, pCov = curve_fit(contrastFn, contrasts, prefPeriPopResp[i],
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]),
                           p0=initialGuess)
    periPrefPopC50.append(pOpt[2])
    periPrefPopRMax.append(pOpt[1])
    periPrefPopN.append(pOpt[3])

    # nonpref peri
    initialGuess = [max(nonprefPeriPopResp[i]), np.median(contrasts),
                    2.0, nonprefPeriPopResp[i][0]]
    pOpt, pCov = curve_fit(contrastFn, contrasts, nonprefPeriPopResp[i],
                           bounds=([0, 0, 0, 0],
                                   [np.inf, np.inf, np.inf, np.inf]),
                           p0=initialGuess)
    periNonprefPopC50.append(pOpt[2])
    periNonprefPopRMax.append(pOpt[1])
    periNonprefPopN.append(pOpt[3])

centerPrefPopC50 = np.array(centerPrefPopC50)
centerPrefPopRMax = np.array(centerPrefPopRMax)
centerPrefPopN = np.array(centerPrefPopN)
centerNonprefPopRMax = np.array(centerNonprefPopRMax)
centerNonprefPopC50 = np.array(centerNonprefPopC50)
centerNonprefPopN = np.array(centerNonprefPopN)
periPrefPopRMax = np.array(periPrefPopRMax)
periPrefPopC50 = np.array(periPrefPopC50)
periPrefPopN = np.array(periPrefPopN)
periNonprefPopC50 = np.array(periNonprefPopC50)
periNonprefPopRMax = np.array(periNonprefPopRMax)
periNonprefPopN = np.array(periNonprefPopN)

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
# alpha_adjusted = alpha / len(data_pairs)
alpha_adjusted = alpha / 2

# Conduct Mann-Whitney U tests for each pair
for data1, data2, comparison_name in data_pairs:
    # stat, p = mannwhitneyu(data1, data2, alternative='two-sided')
    stat, p = mannwhitneyu(data2, data1, alternative='greater')
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
        pOpt, _ = curve_fit(contrastFn, contrasts, meanCount,
                            bounds=([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]),
                            maxfev=10000000)
        c50bootstrapPrefCent.append(pOpt[2])

        # nonpref cent
        bootstrap_sample = nonprefCenterPopResp[
            np.random.choice(nonprefCenterPopResp.shape[0],
                             size=nonprefCenterPopResp.shape[0],
                             replace=True)]
        meanCount = np.mean(bootstrap_sample, axis=0)
        pOpt, _ = curve_fit(contrastFn, contrasts, meanCount,
                            bounds=([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]),
                            maxfev=10000000)
        c50bootstrapNonprefCent.append(pOpt[2])

        # pref peri
        bootstrap_sample = prefPeriPopResp[
            np.random.choice(prefPeriPopResp.shape[0],
                             size=prefPeriPopResp.shape[0],
                             replace=True)]
        meanCount = np.mean(bootstrap_sample, axis=0)
        pOpt, _ = curve_fit(contrastFn, contrasts, meanCount,
                            bounds=([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]),
                            maxfev=10000000)
        c50bootstrapPrefPeri.append(pOpt[2])

        # nonpref peri
        bootstrap_sample = nonprefPeriPopResp[
            np.random.choice(nonprefPeriPopResp.shape[0],
                             size=nonprefPeriPopResp.shape[0],
                             replace=True)]
        meanCount = np.mean(bootstrap_sample, axis=0)
        pOpt, _ = curve_fit(contrastFn, contrasts, meanCount,
                            bounds=([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]),
                            maxfev=10000000)
        c50bootstrapNonprefPeri.append(pOpt[2])

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

########################################################################################################################
########################################################################################################################


# normalization fit.
def weightedNormMTSIG(C_terms, Lp, Lnp, W0, W1, sigma, b, condition, conditionPorNP):
    result = np.zeros_like(C_terms)  # Initialize result array of same shape as C_terms
    for i, C in enumerate(C_terms):
        if conditionPorNP[i]:
            if condition[i]:
                result[i] = ((C * Lp * W1) / ((C * W1) + (1 * sigma))) + b
            else:
                result[i] = ((C * Lp * W0) / ((C * W0) + sigma)) + b
        else:
            if condition[i]:
                result[i] = ((C * Lnp * W1) / ((C * W1) + (1 * sigma))) + b
            else:
                result[i] = ((C * Lnp * W0) / ((C * W0) + sigma)) + b
    return result


def fit_data(C_terms, y_data, condition, conditionPorNP):
    # Define the objective function for curve fitting
    def objective(C_terms, Lp, Lnp, W0, W1, sigma, b):
        return weightedNormMTSIG(C_terms, Lp, Lnp, W0, W1, sigma, b, condition, conditionPorNP)

    # Initial guess for L, w, and sigma
    initial_guess = [100, 50, 1, 0.5, 0.07, y_data[0]]
    bou = [[0, 0, 0, 0, 0, 0],
           [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]]

    # Perform the curve fitting
    popt, pcov = curve_fit(objective, C_terms, y_data, p0=initial_guess, maxfev=100000,
                           bounds=bou)

    return popt, pcov


def weightedNormHeuristicMTSIG(C_terms, Lp, Lnp, W, L0, b, condition, conditionPorNP):
    result = np.zeros_like(C_terms)  # Initialize result array of same shape as C_terms
    for i, C in enumerate(C_terms):
        if conditionPorNP[i]:
            if condition[i]:
                result[i] = ((C * Lp * W) + (L0*b)) / ((C * W) + L0)
            else:
                result[i] = ((C * Lp) + (L0*b)) / (C + L0)
        else:
            if condition[i]:
                result[i] = ((C * Lnp * W) + (L0*b)) / ((C * W) + L0)
            else:
                result[i] = ((C * Lnp) + (L0*b)) / (C + L0)
    return result


def fit_dataWeightedNormHeuristic(C_terms, y_data, condition, conditionPorNP):
    # Define the objective function for curve fitting
    def objective(C_terms, Lp, Lnp, W, L0, b):
        return weightedNormHeuristicMTSIG(C_terms, Lp, Lnp, W, L0, b, condition, conditionPorNP)

    # Initial guess for L, w, and sigma
    initial_guess = [100, 50, 0.5, 0.05, y_data[0]]
    bou = [[0, 0, 0, 0, 0],
           [np.inf, np.inf, np.inf, np.inf, np.inf]]

    # Perform the curve fitting
    popt, pcov = curve_fit(objective, C_terms, y_data, p0=initial_guess, maxfev=100000,
                           bounds=bou)

    return popt, pcov


resp = np.concatenate((prefCenterMean, prefPeriMean, nonprefCenterMean, nonprefPeriMean), axis=0) * 100
cont = np.concatenate((contrasts, contrasts, contrasts, contrasts), axis=0)
conditionArray = np.concatenate((np.full(7, False), np.full(7, True),
                                 np.full(7, False), np.full(7, True)),
                                axis=0)
condition2Array = np.concatenate((np.full(14, True), np.full(14, False)), axis=0)

####### Population
# fit weighted norm
popt, pcov = fit_data(cont, resp, conditionArray, condition2Array)
yPred = weightedNormMTSIG(cont, *popt, conditionArray, condition2Array)
print(r2_score(resp, yPred))
print(popt)

# fit weighted norm heuristic
popt, pcov = fit_dataWeightedNormHeuristic(cont, resp, conditionArray, condition2Array)
yPred = weightedNormHeuristicMTSIG(cont, *popt, conditionArray, condition2Array)
print(r2_score(resp, yPred))
print(popt)

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
    plt.scatter(contrasts, prefCenterMean*100, color='green')

    # pref peri
    yFit = weightedNormMTSIG(xFit, *popt, np.full(100, True), np.full(100, True))
    plt.plot(xFit, yFit, color='green', alpha=0.5, linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[7:], *popt, conditionArray[7:])
    # plt.plot(cont[7:], yFit, color='green', alpha=0.5, linestyle='dashed')
    plt.scatter(contrasts, prefPeriMean*100, color='green', alpha=0.5)

    # nonpref center
    yFit = weightedNormMTSIG(xFit, *popt, np.full(100, False), np.full(100, False))
    plt.plot(xFit, yFit, color='red', linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[:7], *popt, conditionArray[:7])
    # plt.plot(cont[:7], yFit, color='green', linestyle='dashed')
    plt.scatter(contrasts, nonprefCenterMean*100, color='red')

    # nonpref peri
    yFit = weightedNormMTSIG(xFit, *popt, np.full(100, True), np.full(100, False))
    plt.plot(xFit, yFit, color='red', alpha=0.5, linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[7:], *popt, conditionArray[7:])
    # plt.plot(cont[7:], yFit, color='green', alpha=0.5, linestyle='dashed')
    plt.scatter(contrasts, nonprefPeriMean*100, color='red', alpha=0.5)

    plt.xscale('symlog', linthresh=0.1)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(loc='lower right')
    plt.show()


# plot population RF weighted normalization (heuristic)
for filler in range(1):
    popt, pcov = fit_dataWeightedNormHeuristic(cont, resp, conditionArray, condition2Array)
    yPred = weightedNormHeuristicMTSIG(cont, *popt, conditionArray, condition2Array)
    print(r2_score(resp, yPred))

    # plotting
    xFit = np.logspace(-1, 2, 100)

    # pref center
    yFit = weightedNormHeuristicMTSIG(xFit, *popt, np.full(100, False), np.full(100, True))
    plt.plot(xFit, yFit, color='green', linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[:7], *popt, conditionArray[:7])
    # plt.plot(cont[:7], yFit, color='green', linestyle='dashed')
    plt.scatter(contrasts, prefCenterMean*100, color='green')

    # pref peri
    yFit = weightedNormHeuristicMTSIG(xFit, *popt, np.full(100, True), np.full(100, True))
    plt.plot(xFit, yFit, color='green', alpha=0.5, linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[7:], *popt, conditionArray[7:])
    # plt.plot(cont[7:], yFit, color='green', alpha=0.5, linestyle='dashed')
    plt.scatter(contrasts, prefPeriMean*100, color='green', alpha=0.5)

    # nonpref center
    yFit = weightedNormHeuristicMTSIG(xFit, *popt, np.full(100, False), np.full(100, False))
    plt.plot(xFit, yFit, color='red', linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[:7], *popt, conditionArray[:7])
    # plt.plot(cont[:7], yFit, color='green', linestyle='dashed')
    plt.scatter(contrasts, nonprefCenterMean*100, color='red')

    # nonpref peri
    yFit = weightedNormHeuristicMTSIG(xFit, *popt, np.full(100, True), np.full(100, False))
    plt.plot(xFit, yFit, color='red', alpha=0.5, linestyle='dashed')
    # yFit = weightedNormMTSIG(cont[7:], *popt, conditionArray[7:])
    # plt.plot(cont[7:], yFit, color='green', alpha=0.5, linestyle='dashed')
    plt.scatter(contrasts, nonprefPeriMean*100, color='red', alpha=0.5)

    plt.xscale('symlog', linthresh=0.1)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(loc='lower right')
    plt.show()


######## individual neuron fits
# fit weighted norm
r2ScoresWeighted = []
r2ScoresHeuristic = []
for unit in range(len(unitList)):
    resp = np.concatenate((prefCenterPopResp[unit], prefPeriPopResp[unit],
                           nonprefCenterPopResp[unit], nonprefPeriPopResp[unit]),
                          axis=0) * 100
    popt, pcov = fit_data(cont, resp, conditionArray, condition2Array)
    yPred = weightedNormMTSIG(cont, *popt, conditionArray, condition2Array)
    r2ScoresWeighted.append((r2_score(resp, yPred)))

    popt, pcov = fit_data(cont, resp, conditionArray, condition2Array)
    yPred = weightedNormHeuristicMTSIG(cont, *popt, conditionArray, condition2Array)
    r2ScoresHeuristic.append((r2_score(resp, yPred)))

r2ScoresWeighted = np.array(r2ScoresWeighted)
r2ScoresHeuristic = np.array(r2ScoresHeuristic)

# plot scatter of r2 scores
fig, ax1 = plt.subplots()

ax1.scatter(r2ScoresWeighted, r2ScoresHeuristic)
ax1.set_xlabel('r2 Scores Weighted')
ax1.set_ylabel('r2 Scores Heuristic')

x_min, x_max = ax1.get_xlim()
y_min, y_max = ax1.get_ylim()
min_val = min(x_min, y_min)
max_val = max(x_max, y_max)
ax1.plot([min_val, max_val], [min_val, max_val], 'k--', label='Unity Line', zorder=0)

# Reset limits to ensure the unity line fits perfectly
ax1.set_xlim([min_val, max_val])
ax1.set_ylim([min_val, max_val])
stat, p_value = wilcoxon(r2ScoresWeighted, r2ScoresHeuristic, alternative='two-sided')
print(p_value)

plt.show()



### BACKUP
# normalization fit.
def weightedNormMTSIG(C_terms, Lp, Lnp, W, sigma, b, condition, conditionPorNP):
    result = np.zeros_like(C_terms)  # Initialize result array of same shape as C_terms
    for i, C in enumerate(C_terms):
        if conditionPorNP[i]:
            if condition[i]:
                result[i] = ((C * Lp * W) / ((C * W) + sigma)) + b
            else:
                result[i] = ((C * Lp) / (C + sigma)) + b
        else:
            if condition[i]:
                result[i] = ((C * Lnp * W) / ((C * W) + sigma)) + b
            else:
                result[i] = ((C * Lnp) / (C + sigma)) + b
    return result


def fit_data(C_terms, y_data, condition, conditionPorNP):
    # Define the objective function for curve fitting
    def objective(C_terms, Lp, Lnp, W, sigma, b):
        return weightedNormMTSIG(C_terms, Lp, Lnp, W, sigma, b, condition, conditionPorNP)

    # Initial guess for L, w, and sigma
    initial_guess = [100, 50, 0.5, 0.05, y_data[0]]
    bou = [[0, 0, 0, 0, 0],
           [np.inf, np.inf, np.inf, np.inf, np.inf]]

    # Perform the curve fitting
    popt, pcov = curve_fit(objective, C_terms, y_data, p0=initial_guess, maxfev=100000,
                           bounds=bou)

    return popt, pcov


def weightedNormHeuristicMTSIG(C_terms, Lp, Lnp, W, L0, b, condition, conditionPorNP):
    result = np.zeros_like(C_terms)  # Initialize result array of same shape as C_terms
    for i, C in enumerate(C_terms):
        if conditionPorNP[i]:
            if condition[i]:
                result[i] = ((C * Lp * W) + (L0*b)) / ((C * W) + L0)
            else:
                result[i] = ((C * Lp) + (L0*b)) / (C + L0)
        else:
            if condition[i]:
                result[i] = ((C * Lnp * W) + (L0*b)) / ((C * W) + L0)
            else:
                result[i] = ((C * Lnp) + (L0*b)) / (C + L0)
    return result


def fit_dataWeightedNormHeuristic(C_terms, y_data, condition, conditionPorNP):
    # Define the objective function for curve fitting
    def objective(C_terms, Lp, Lnp, W, L0, b):
        return weightedNormHeuristicMTSIG(C_terms, Lp, Lnp, W, L0, b, condition, conditionPorNP)

    # Initial guess for L, w, and sigma
    initial_guess = [100, 50, 0.5, 0.05, y_data[0]]
    bou = [[0, 0, 0, 0, 0],
           [np.inf, np.inf, np.inf, np.inf, np.inf]]

    # Perform the curve fitting
    popt, pcov = curve_fit(objective, C_terms, y_data, p0=initial_guess, maxfev=100000,
                           bounds=bou)

    return popt, pcov


