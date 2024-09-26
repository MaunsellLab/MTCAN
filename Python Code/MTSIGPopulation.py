"""
MTSigma population main analysis script. This script will extract spike counts during
valid trials from units active in a session and put them in the appropriate
matrix. Then we will add this to a population matrix to measure how c50 changes
with rMax at a population level.

Chery - July 2024

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

# # non-curated file List Akshan (C50 at center and peri for pref unmatched)
fileList = ['Akshan_240603', 'Akshan_240606', 'Akshan_240610', 'Akshan_240701',
            'Akshan_240903', 'Akshan_240903_2', 'Akshan_240905', 'Akshan_240906',
            'Akshan_240909', 'Akshan_240911', 'Akshan_240913', 'Akshan_240916',
            'Akshan_240917', 'Akshan_240922', 'Akshan_240924']
unitList = ['240603_167', '240606_176', '240610_169', '240701_1', '240903_3',
            '240903_8_2', '240905_4', '240906_25', '240909_10', '240911_17',
            '240913_1', '240913_27', '240913_36', '240913_47', '240913_48',
            '240916_7', '240917_44', '240922_6', '240924_22', '240924_35',
            '240924_21']

# # curated file List Akshan (C50 at center and peri for pref unmatched)
# fileList = ['Akshan_240603', 'Akshan_240606', 'Akshan_240610', 'Akshan_240903_2',
#             'Akshan_240905', 'Akshan_240906', 'Akshan_240913']
# unitList = ['240603_167', '240606_176', '240610_169', '240903_8_2', '240905_4',
#             '240906_25', '240913_27', '240913_48', '240913_36']


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

    # fit CRF

    for unit in sessionGoodUnits:
        count = np.where(units == unit)[0][0]
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
        prefCenterResp = meanSpikeReshaped[count][0][0]
        prefCenterResp = np.insert(prefCenterResp, 0, baselineResp)
        prefPeriResp = meanSpikeReshaped[count][1][0]
        prefPeriResp = np.insert(prefPeriResp, 0, baselineResp)
        prefCenterPopResp.append(prefCenterResp/np.max(meanSpikeReshaped[count]))
        prefPeriPopResp.append(prefPeriResp/np.max(meanSpikeReshaped[count]))

        # create population avg CRF curve for nonpref resp at center and peri
        nonprefCenterResp = meanSpikeReshaped[count][0][1]
        nonprefCenterResp = np.insert(nonprefCenterResp, 0, baselineResp)
        nonprefPeriResp = meanSpikeReshaped[count][1][1]
        nonprefPeriResp = np.insert(nonprefPeriResp, 0, baselineResp)
        nonprefCenterPopResp.append(nonprefCenterResp/np.max(meanSpikeReshaped[count]))
        nonprefPeriPopResp.append(nonprefPeriResp/np.max(meanSpikeReshaped[count]))

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

        rMaxPerChange.append(prefRMaxChange)
        c50PerChange.append(prefC50Change)
        # rMaxPerChange.append(npRMaxChange)
        # c50PerChange.append(npC50Change)

    # to open another file in the loop
    os.chdir('../../Python Code')


#################################################### Plotting ####################################################
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


# plotting figure
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
    print(pOpt)

    # plot pref peri
    initialGuess = [prefPeriMean[0], max(prefPeriMean), np.median(contrasts), 2.0]
    # pOpt, pCov = curve_fit(contrastFn, contrasts, prefPeriMean, bounds=([prefPeriMean[0], 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
    pOpt, pCov = curve_fit(contrastFn, contrasts, prefPeriMean,
                           bounds=([prefPeriMean[0], 0.8812, 0, 0], [np.inf, 0.8812+0.001, np.inf, np.inf]))
    print(pOpt)
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
    ax[1].set_ylim(0, 30)
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
            if condition[i]:  # Check condition for each C term
                result[i] = (C * Lp * w**1) / ((C * w) + sigma) + b
            else:
                result[i] = (C * Lp) / (C + sigma) + b
        else:
            if condition[i]:  # Check condition for each C term
                result[i] = (C * Lnp * w**1) / ((C * w) + sigma) + b
            else:
                result[i] = (C * Lnp) / (C + sigma) + b
    return result


# Example of the fitting procedure using curve_fit
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



###### WORKING VERSION

# Define the conditional model
def weightedNormMTSIG(C_terms, L, w, sigma, b, condition):
    result = np.zeros_like(C_terms)  # Initialize result array of same shape as C_terms
    for i, C in enumerate(C_terms):
        if condition[i]:  # Check condition for each C term
            result[i] = (C * L * w**1) / ((C * w) + sigma) + b
        else:
            result[i] = (C * L) / (C + sigma) + b
    return result


# Example of the fitting procedure using curve_fit
def fit_data(C_terms, y_data, condition):
    # Define the objective function for curve fitting
    def objective(C_terms, L, w, sigma, b):
        return weightedNormMTSIG(C_terms, L, w, sigma, b, condition)

    # Initial guess for L, w, and sigma
    initial_guess = [.5, 0.5, 0.05, resp[0]]

    # Perform the curve fitting
    popt, pcov = curve_fit(objective, C_terms, y_data, p0=initial_guess)

    return popt, pcov


resp = np.concatenate((nonprefCenterMean, nonprefPeriMean), axis=0)
cont = np.concatenate((contrasts, contrasts), axis=0)
conditionArray = np.concatenate((np.full(7, False), np.full(7, True)), axis=0)

popt, pcov = fit_data(cont, resp, conditionArray)


# plotting
xFit = np.logspace(-1, 2, 100)

yFit = weightedNormMTSIG(xFit, *popt, np.full(100, False))
plt.plot(xFit, yFit, color='green', linestyle='dashed')
# yFit = weightedNormMTSIG(cont[:7], *popt, conditionArray[:7])
# plt.plot(cont[:7], yFit, color='green', linestyle='dashed')
plt.plot(contrasts, prefCenterMean, color='green')

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


yFit = weightedNormMTSIG(xFit, *popt, np.full(100, True))
plt.plot(xFit, yFit, color='green', alpha=0.5, linestyle='dashed')
# yFit = weightedNormMTSIG(cont[7:], *popt, conditionArray[7:])
# plt.plot(cont[7:], yFit, color='green', alpha=0.5, linestyle='dashed')
plt.plot(contrasts, prefPeriMean, color='green', alpha=0.5)

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

# plt.xscale('symlog', linthresh=0.1)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend(loc='lower right')
plt.show()





