"""
MTND main analysis script. This script will extract spike counts during
valid trials from units active in a session and put them in the appropriate
matrix. From there, this script plots how normalization changes with distance.
We also plot the preferred response across the diameter of the receptive field

Chery - June 2023
"""

# Imports
from usefulFns import *

# subpar sessions: 230606, 230706, 230712

p0 = time.time()

# fileList = ['230720']

# # master file list
# fileList = ['230607', '230608', '230612', '230613', '230615',
#             '230620', '230626', '230627', '230628', '230630',
#             '230707', '230710', '230711', '230719', '230720]

# # master unit list
# unitList = ['230607_147', '230607_146', '230607_144', '230607_140', '230607_126',
#             '230607_125', '230608_72', '230608_137', '230608_143', '230608_145',
#             '230608_146', '230608_149', '230612_213', '230613_187', '230615_83',
#             '230615_164', '230615_166', '230620_58', '230626_131', '230627_147'
#             '230627_141', '230627_106', '230627_71', '230627_70', '230627_14',
#             '230628_143', '230630_50', '230707_140', '230710_137', '230711_45',
#             '230711_50', '230718_178', '230719_147', '230719_156','230720_149']

# # pref + null
# fileList = ['230607', '230608', '230612', '230613', '230615']

# unitList = ['230607_147', '230607_146', '230607_144', '230607_140', '230607_126',
#             '230607_125', '230608_72', '230608_137', '230608_143', '230608_145',
#             '230608_146', '230608_149', '230612_213', '230613_187', '230615_83',
#             '230615_164', '230615_166']

# pref + nonpref
fileList = ['230620', '230626', '230627', '230628', '230630',
            '230707', '230710', '230711', '230718', '230719',
            '230720']

unitList = ['230620_58', '230626_131', '230627_147', '230627_141', '230627_106',
            '230627_71', '230627_70', '230627_14', '230628_143', '230630_50',
            '230707_140', '230710_137', '230711_45', '230711_50', '230718_178',
            '230719_147', '230719_156', '230720_149']

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
offsetDegSepNormPop = []
adaptationMat = np.zeros((100, 4, 100))
adaptationMat[:] = np.nan
adaptC = 0
allBlocksDone = []

for file in fileList:
    # Load relevant file here with pyMat reader
    monkeyName = 'Meetz'
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
        ax.set_rticks(np.linspace(0, np.max(r)*1.1, 5))

        # RF Heatmap
        ax1 = fig.add_subplot(gs01[0, 0])

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

        # normalization vs distance (PN, NP)
        ax2 = fig.add_subplot(gs02[0, 0])
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
        pCentPred = (prefOnly[0] + nullOnly[1:]) / (
                1 + prefOnly[1:] / np.max(prefOnly))
        npCentPred = (nullOnly[0] + prefOnly[1:]) / (
                1 + prefOnly[1:] / np.max(prefOnly))
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

        # pref response across diameter of RF
        ax3 = fig.add_subplot(gs03[0, 0])
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
        ax3.plot(x, prefTransect, color='black')
        ax3.errorbar(x, prefTransect, yerr=prefTransectSEM, fmt='o', ecolor='black',
                     color='black', markersize=2)
        ax3.set_xlabel('stimulus offset positions')
        ax3.set_ylabel('Response (spikes/sec)')
        ax3.set_ylim(bottom=0)
        ax3.axhline(y=spon, linestyle='--', color='blue', alpha=0.8)

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
            # PN
            centerP = meanSpikeReshaped[count, 2, 0]
            offsetN = meanSpikeReshaped[count, 0, i*2-1]
            PN = meanSpikeReshaped[count, 2, i*2-1]
            pnNMI = ((centerP + offsetN) - PN) / ((centerP + offsetN) + PN)
            pCentNMI.append(pnNMI)

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

            # NN
            NN = meanSpikeReshaped[count, 1, i*2-1]
            tempNMI = ((centerN + offsetN) - NN) / ((centerN + offsetN) + NN)
            nnNMI.append(tempNMI)

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
            spon = meanSpikeReshaped[count, 0, 0] / np.max(meanSpikeReshaped[count])
            nullOnly = (np.concatenate(([meanSpikeReshaped[count, 1, 0]],
                                       meanSpikeReshaped[count, 0, 1::2]), axis=0) /
                        np.max(meanSpikeReshaped[count]))
            prefOnly = (np.concatenate(([meanSpikeReshaped[count, 2, 0]],
                                       meanSpikeReshaped[count, 0, 2::2]), axis=0) /
                        np.max(meanSpikeReshaped[count]))
            NP = meanSpikeReshaped[count, 1, 2::2] / np.max(meanSpikeReshaped[count])
            PN = meanSpikeReshaped[count, 2, 1::2] / np.max(meanSpikeReshaped[count])
            PP = meanSpikeReshaped[count, 2, 2::2] / np.max(meanSpikeReshaped[count])
            NN = meanSpikeReshaped[count, 1, 1::2] / np.max(meanSpikeReshaped[count])
            prefTransect = prefTransect / np.max(meanSpikeReshaped[count])

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

numSteps = 4
numUnits = len(sponNormalized)
offsetDegSepNormPop = offsetDegSepNormPop.reshape((numUnits, numSteps*2+1))
meanSpon = np.mean(sponNormalized)
semSpon = np.std(sponNormalized) / np.sqrt(numUnits)

# figure
fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4)
fig.set_size_inches(16, 8)

# bin size = n, exponent for prediction = exp
n = 20
exp = 1
fig.suptitle(f'population responses pref+null neurons, exp={exp}')

# transect response
transX = np.concatenate((np.arange(-numSteps, 0),
                         np.arange(0, numSteps + 1)),
                        axis=0)
transectMean = np.mean(transectNormalized, axis=0)
transectSEM = np.std(transectNormalized, axis=0) / np.sqrt(numUnits)
ax1.plot(transX, transectMean, color='black', label='Pref')
ax1.errorbar(transX, transectMean, yerr=transectSEM, fmt='o', ecolor='black',
             color='black', markersize=2)
ax1.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
ax1.fill_between(x=transX, y1=meanSpon+semSpon, y2=meanSpon-semSpon,
                 color='blue', alpha=0.2)
ax1.set_ylabel('Normalized Response')

ax1.set_xlabel('Offset Positions')
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
# pCentPred = ((prefMean[0] + nonprefMean[1:])**1) / (
#             1 + prefMean[1:] / np.max(prefMean))
# npCentPred = ((nonprefMean[0] + prefMean[1:])**1) / (
#             1 + prefMean[1:] / np.max(prefMean))
pCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / (
            ((prefMean[0] + prefMean[1:])**1) / np.max(transectMean))
npCentPred = ((nonprefMean[0] + prefMean[1:])**exp) / (
             ((prefMean[0] + prefMean[1:])**1) / np.max(transectMean))
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

# NMI pop plot
pnNMIPopMean = np.mean(pnNMIPop, axis=0)
pnNMISEM = np.std(pnNMIPop, axis=0) / np.sqrt(numUnits)
npNMIPopMean = np.mean(npNMIPop, axis=0)
npNMISEM = np.std(npNMIPop, axis=0) / np.sqrt(numUnits)
ppNMIPopMean = np.mean(ppNMIPop, axis=0)
ppNMISEM = np.std(ppNMIPop, axis=0) / np.sqrt(numUnits)
nnNMIPopMean = np.mean(nnNMIPop, axis=0)
nnNMISEM = np.std(nnNMIPop, axis=0) / np.sqrt(numUnits)

ax4.plot(xNorm, pnNMIPopMean, color='green', label='PN')
ax4.errorbar(xNorm, pnNMIPopMean, yerr=pnNMISEM, fmt='o', ecolor='green',
             color='green', markersize=2)
ax4.scatter(xNorm, pnNMIPopMean, color='green')
ax4.plot(xNorm, npNMIPopMean, color='red', label='NP')
ax4.errorbar(xNorm, npNMIPopMean, yerr=npNMISEM, fmt='o', ecolor='red',
             color='red', markersize=2)
ax4.scatter(xNorm, npNMIPopMean, color='red')
ax4.plot(xNorm, ppNMIPopMean, color='black', label='PP', linestyle='--')
ax4.errorbar(xNorm, ppNMIPopMean, yerr=ppNMISEM, fmt='o', ecolor='black',
             color='black', markersize=2)
ax4.scatter(xNorm, ppNMIPopMean, color='black')
ax4.plot(xNorm, nnNMIPopMean, color='grey', label='NN', linestyle='--')
ax4.errorbar(xNorm, nnNMIPopMean, yerr=nnNMISEM, fmt='o', ecolor='grey',
             color='grey', markersize=2)
ax4.scatter(xNorm, nnNMIPopMean, color='grey')
ax4.set_xlabel('stimulus offset position')
ax4.set_ylabel('NMI')
ax4.set_ylim(bottom=0)
ax4.set_title('NMI, Ni and Maunsell, 2017')

# raw vs predicted norm PN, NP, PP, NN
pCentPred = ((np.reshape(prefNormalized[:, 0], (numUnits, 1)) + nonprefNormalized[:, 1:])**exp) / (
            (np.reshape(prefNormalized[:, 0], (numUnits, 1)) + prefNormalized[:, 1:]) / (
             (np.reshape(np.max(transectNormalized, axis=1), (numUnits, 1)))))
meanPCentPred = ((prefMean[0] + nonprefMean[1:])**exp) / (
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
ax5.set_ylim([0, 1.5])
ax5.set_xlim([0, 1.5])
ax5.set_xlabel('Real Response')
ax5.set_ylabel('Predicted Response')
ax5.set_title('real vs predicated normalization response')
line = lines.Line2D([0, 1], [0, 1], color='black')
transform = ax5.transAxes
line.set_transform(transform)
ax5.add_line(line)
ax5.legend()

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
binSEMPrefResp = np.array([np.std(i) for i in equalBinsPrefResp]) / n
binMeanNonprefResp = np.array([np.mean(i) for i in equalBinsNonprefResp])
binSEMNonpreResp = np.array([np.std(i) for i in equalBinsNonprefResp]) / n

# pn/np response
sep = offsetDegSepNormPop[:, 5:].reshape(numUnits*numSteps)
pnResp = pnNormalized.reshape(numUnits*numSteps)
npResp = npNormalized.reshape(numUnits*numSteps)
sortIndex = np.argsort(sep)
sortedSep = sep[sortIndex]
sortedPNResp = pnResp[sortIndex]
sortedNPResp = npResp[sortIndex]
# manual bins - equally populated bins
equalBinsSep = [sortedSep[i:i + n] for i in range(0, len(sortedSep), n)]
equalBinsPNResp = [sortedPNResp[i:i + n] for i in range(0, len(sortedPNResp), n)]
equalBinsNPResp = [sortedNPResp[i:i + n] for i in range(0, len(sortedNPResp), n)]
binMeanSep = np.array([np.mean(i) for i in equalBinsSep])
binMeanPNResp = np.array([np.mean(i) for i in equalBinsPNResp])
binSEMPNResp = np.array([np.std(i) for i in equalBinsPNResp]) / n
binMeanNPResp = np.array([np.mean(i) for i in equalBinsNPResp])
binSEMNPResp = np.array([np.std(i) for i in equalBinsNPResp]) / n

# pn/np pred
pnPred = ((binMeanPrefResp[0] + binMeanNonprefResp[1:])**exp) / (
    (binMeanPrefResp[0] + binMeanPrefResp[1:]) / np.max(binMeanTransectResp))
npPred = ((binMeanNonprefResp[0] + binMeanPrefResp[1:])**exp) / (
    (binMeanPrefResp[0] + binMeanPrefResp[1:]) / np.max(binMeanTransectResp))

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
ax7.set_xlabel('Normalized Separation')
ax7.set_ylabel('Normalized Binned Resp')
ax7.set_ylim(bottom=0)
ax7.set_title(f'Equal sized bins ({n}) of normalized separation vs Response')
ax7.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
ax7.fill_between(x=np.linspace(0, 1, 10), y1=meanSpon+semSpon,
                 y2=meanSpon-semSpon, color='blue', alpha=0.2)

ax8.set_visible(False)

plt.tight_layout()
plt.show()

# adaptation mat
adaptationMean = np.nanmean(adaptationMat[:numUnits, :, :50], axis=0)
plt.plot(adaptationMean[0]); plt.show()


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

"""

# weighting by max response of neuron to direction tuning
count = 6

x = np.concatenate((np.arange(-numSteps, 0),
                    np.arange(0, numSteps + 1)),
                   axis=0)
xNorm = np.arange(1, numSteps+1)
transReverse = transectMeanSpike[count, ::-1]
prefTransect = np.concatenate((transReverse,
                               [meanSpikeReshaped[count, 2, 0]],
                               meanSpikeReshaped[count, 0, 2::2]),
                              axis=0)


pCentPred = (meanSpikeReshaped[count, 2, 0] + meanSpikeReshaped[count, 0, 1::2]) / (
            1 + meanSpikeReshaped[count, 0, 2::2] / 32.41106719)

npCentPred = (meanSpikeReshaped[count, 1, 0] + meanSpikeReshaped[count, 0, 2::2]) / (
            1 + meanSpikeReshaped[count, 0, 2::2] / 32.41106719)

# pCentPred = (meanSpikeReshaped[count, 2, 0] + meanSpikeReshaped[count, 0, 2::2]) / (
#             1 + meanSpikeReshaped[count, 0, 2::2] / 32.41106719)
#
# npCentPred = (meanSpikeReshaped[count, 1, 0] + meanSpikeReshaped[count, 0, 1::2]) / (
#             1 + meanSpikeReshaped[count, 0, 2::2] / 32.41106719)

plt.plot(x, prefTransect, color='black')
plt.plot(xNorm, pCentPred, color='green', linestyle='--')
plt.plot(xNorm, npCentPred, color='red', linestyle='--')

plt.show()


plt.plot(transectOneOverTime, color='green')
plt.plot(testOneOverTime, color='gold')
plt.plot(centerPrefOverTime, color='black')
plt.show()
