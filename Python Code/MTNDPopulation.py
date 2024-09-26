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

# # pref + nonpref (Meetz) (n=18)
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
#             'Akshan_240704', 'Akshan_240705', 'Akshan_240709']
#
# unitList = ['240530_55', '240603_167', '240606_176', '240607_137', '240607_170',
#             '240607_179', '240611_136', '240613_138', '240628_19', '240703_19',
#             '240703_91', '240704_86', '240705_74', '240705_162', '240709_53',
#             '240709_108']


# # both monkeys (Meetz + Akshan)
fileList = ['Meetz_230620', 'Meetz_230626', 'Meetz_230627', 'Meetz_230628',
            'Meetz_230630', 'Meetz_230707', 'Meetz_230710', 'Meetz_230711',
            'Meetz_230718', 'Meetz_230719', 'Meetz_230720', 'Akshan_240530',
            'Akshan_240603', 'Akshan_240606', 'Akshan_240607', 'Akshan_240611',
            'Akshan_240613', 'Akshan_240628', 'Akshan_240703', 'Akshan_240704',
            'Akshan_240705', 'Akshan_240709']

unitList = ['230620_58', '230626_131', '230627_147', '230627_141', '230627_106',
            '230627_71', '230627_70', '230627_14', '230628_143', '230630_50',
            '230707_140', '230710_137', '230711_45', '230711_50', '230718_178',
            '230719_147', '230719_156', '230720_149', '240530_55', '240603_167',
            '240606_176', '240607_137', '240607_170', '240607_179', '240611_136',
            '240613_138', '240628_19', '240703_19', '240703_91', '240704_86',
            '240705_74', '240705_162', '240709_53', '240709_108']


# population arrays
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

                    # only units that are considered good for population analysis
                    for unit in sessionGoodUnits:
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
    for unit in sessionGoodUnits:
        count = np.where(units == unit)[0][0]

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

    # to open another file in the loop
    os.chdir('../../Python Code')

print(time.time()-p0)

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

    # pCentPred = (((prefMean[0]) + (nonprefMean[1:])) / (
    #             weights[0] + weights[1:]))
    #
    # npCentPred = (((nonprefMean[0]) + (prefMean[1:])) / (
    #              weights[0] + weights[1:]))

    pCentPred = (((prefMean[0]) + (nonprefMean[0] * (weights[1:]**2))) / (
                weights[0] + weights[1:]) + 0.08)

    npCentPred = (((nonprefMean[0]) + (prefMean[0] * (weights[1:]**2))) / (
                 weights[0] + weights[1:]) + 0.08)

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

    ppCentPred = ((prefMean[0] * weights[0]) + (prefMean[1:] * weights[1:]) + (meanSpon * (weights[0] + weights[1:] + 0.07) )) / (
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
plt.hist(data3, bins=20, alpha=1, label='pnNMIPop', color='blue')
plt.hist(data4, bins=20, alpha=1, label='npNMIPop', color='blue')
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
