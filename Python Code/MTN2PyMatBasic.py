
# imports
import seaborn as sns
import numpy as np
import numpy.ma as ma
from usefulFns import *
from itertools import combinations
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time

# for loop to run through all files
t0 = time.time()
totUnits = []
totR2 = []
fileList = ['221010', '221108', '221110', '221115', '221117', '221124',
            '221128', '221010', '221208']

for fileIterator in fileList:

    # Load relevant file here with pyMat reader
    monkeyName = 'Meetz'
    seshDate = fileIterator
    fileName = f'{monkeyName}_{seshDate}_MTNC_Spikes.mat'

    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

    if not os.path.exists('Normalization'):
        os.makedirs('Normalization')
    os.chdir('Normalization/')

    # list of indices of correctTrials (non-instruct, valid trialCertify)
    corrTrials = correctTrialsMTX(allTrials)

    # generate list of unique active units, and their channel
    units = activeUnits('spikeData', allTrials)
    unitCluster = allTrials[corrTrials[0]]['spikeTempInfo']['cgs']
    unitsChannel = unitsInfo(units, corrTrials, allTrials)

    # change stimDesc to be list of dictionaries
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        nStim = len(currTrial['stimDesc']['data']['listType'])
        currTrial['stimDesc']['data'] = [{k:v[i] for k,v in currTrial['stimDesc']['data'].items()}
                                        for i in range(nStim)]

    # assert: stim sequence list is frozen
    seqList = []
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        stimDesc = currTrial['stimDesc']['data']
        for stim in stimDesc:
            if stim['stimLoc'] == 0 and stim['listType'] == 1:
                if len(seqList) < 49:
                    seqList.append(stim['stimIndex'])
                    seqArr = np.array(seqList)
                    lastIndex = stim['stimIndex']
                else:
                    posLastIndex = np.where(seqArr==lastIndex)[0][0]
                    if posLastIndex == len(seqArr)-1:
                        if stim['stimIndex'] != seqArr[0]:
                            print('out of sequence')
                        else:
                            lastIndex = stim['stimIndex']
                    else:
                        if stim['stimIndex'] != seqArr[posLastIndex+1]:
                            print('out of sequence')
                        else:
                            lastIndex = stim['stimIndex']

    # assert: are there correct trials without spikeData
    noSpikeData = []
    for trialCount, currTrial in enumerate(allTrials):
        trial = currTrial['trial']['data']
        extendedEOT = currTrial['extendedEOT']['data']
        if extendedEOT == 0 and trial['instructTrial'] != 1:
            if 'spikeData' not in currTrial:
                noSpikeData.append(trialCount)

    # assert: frame consistency during stimlus duration
    frameRateHz = header['frameRateHz']['data']
    stimDurFrame = []
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
        trueStimDurMS = np.int32(np.around(1000/frameRateHz * stimDurFrame[0]))

    # generates a dictionary, numpy array, and Pandas Dataframe of stim Index
    # and corresponding directions/contrasts
    stimIndexDict = {}
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        stimDesc = currTrial['stimDesc']['data']
        for stim in stimDesc:
            if stim['stimLoc'] != 2:
                index = stim['stimIndex']
                if index not in stimIndexDict:
                    stimIndexDict[index] = {}
                    if stim['stimLoc'] not in stimIndexDict[index]:
                        stimIndexDict[index][stim['stimLoc']] = \
                        {'direction': stim['directionDeg'],
                            'contrast': np.around(stim['contrast'],2)}
                else:
                    if stim['stimLoc'] not in stimIndexDict[index]:
                        stimIndexDict[index][stim['stimLoc']] = \
                        {'direction': stim['directionDeg'],
                            'contrast': np.around(stim['contrast'],2)}
    stimIndexArray = np.zeros((49,4))
    for i in range(len(stimIndexDict)):
        stimIndexArray[i][0] = stimIndexDict[i][0]['direction']
        stimIndexArray[i][1] = stimIndexDict[i][0]['contrast']
        stimIndexArray[i][2] = stimIndexDict[i][1]['direction']
        stimIndexArray[i][3] = stimIndexDict[i][1]['contrast']
    stimIndexDF = pd.DataFrame(stimIndexArray, columns=['loc0 Direction', 'loc0 Contrast',
                                                'loc1 Direction', 'loc1 Contrast'])

    # initialize lists/arrays/dataframes for counting spikeCounts and for analysis
    blocksDone = allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone']
    highContrast, zeroContrast = max(stimIndexDF['loc0 Contrast'].unique()), \
                                 min(stimIndexDF['loc0 Contrast'].unique())
    zeroDir = 0
    dirArray = np.array([0,60,120,180,240,300])
    angleMat = np.arange(180,900,60)
    spikeCountMat = np.zeros((len(units),blocksDone+1,49))
    spikeCountLong = []
    sponSpikeCountLong = []
    onLatency = 25/1000  # time in MS for counting window latency after stim on
    offLatency = 100/1000  # time in MS for counting window latency after stim off
    histPrePostMS = 100  # 100ms window pre/post stimulus on/off
    sponWindowMS = 50  # 50ms window before stimulus onset
    spikeHists = np.zeros((len(units),49, trueStimDurMS + (2*histPrePostMS+1)))
    stimIndexCount = np.zeros(49)

    # insert spike counts into matrix of unique stimulus sets
    for corrTrial in corrTrials:
        currTrial = allTrials[corrTrial]
        if 'spikeData' in currTrial:
            stimDesc = currTrial['stimDesc']['data']
            stim1TimeS = currTrial['taskEvents']['stimulusOn']['time'][0]
            fixateTimeS = currTrial['taskEvents']['fixate']['time']
            for stim in stimDesc:
                if stim['stimLoc'] == 0 and stim['listType'] == 1:
                    stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'])
                                    / 1000) + stim1TimeS
                    stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'])
                                    / 1000) + stim1TimeS
                    stimIndex = np.int32(stim['stimIndex'])
                    stCount = int(stimIndexCount[stimIndex])
                    stimIndexCount[stimIndex] += 1
                    for unitCount, unit in enumerate(units):
                        if unit in currTrial['spikeData']['unit']:
                            unitIndex = np.where(currTrial['spikeData']['unit'] == unit)[0]
                            # added 50ms onset latency for spike counts (100 for offset)
                            unitTimeStamps = currTrial['spikeData']['timeStamp'][unitIndex]
                            stimSpikes = np.where((unitTimeStamps >= (stimOnTimeS + onLatency)) &
                                         (unitTimeStamps <= (stimOffTimeS + offLatency)))[0]
                            spikeCountMat[unitCount][stCount][stimIndex] \
                                = len(stimSpikes)
                            spikeCountLong.append([unit,
                                                   stimIndex,
                                                   stimIndexCount[stimIndex],
                                                   len(stimSpikes)])

                            # Spontaneous Spikes
                            sponSpikes = np.where((unitTimeStamps >= (stimOnTimeS-(sponWindowMS/1000)))
                                                  & (unitTimeStamps <= stimOnTimeS))[0]
                            sponSpikeCountLong.append([unit, len(sponSpikes)])

                            # PSTHs
                            stimOnPreSNEV = stimOnTimeS - (histPrePostMS/1000)
                            stimOffPostSNEV = stimOffTimeS + (histPrePostMS/1000)
                            histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)
                                                            & (unitTimeStamps <= stimOffPostSNEV)
                                                             )] - stimOnPreSNEV
                            histStimSpikes = np.int32(histStimSpikes*1000)
                            # spikeHists[unitCount, stimIndex, histStimSpikes] += 1
                            spikeHists[unitCount][stimIndex][histStimSpikes] += 1

    # mean, SEM, and reshaping of spikeCount matrices
    # create pandas dataframe of spikeCount with corresponding unit, stimIndex
    spikeCountDF = pd.DataFrame(spikeCountLong, columns=['unit', 'stimIndex',
                                                         'stimCount', 'stimSpikes'])
    sponSpikeCountDF = pd.DataFrame(sponSpikeCountLong, columns=['unit', 'sponSpikes'])
    sponSpikesMean = np.zeros(len(units)) # spikes in 50ms window
    sponSpikesSEM = np.zeros(len(units))
    for unitCount, unit in enumerate(units):
        sponSpikesMean[unitCount] = sponSpikeCountDF.loc[sponSpikeCountDF['unit'] == unit].mean()[1]
        sponSpikesSEM[unitCount] = sponSpikeCountDF.loc[sponSpikeCountDF['unit'] == unit].sem()[1]
    meanSpike = np.mean(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSD = np.std(spikeCountMat[:, :blocksDone, :], axis=1)
    spikeCountSEM = spikeCountSD/np.sqrt(blocksDone)
    meanSpikeReshaped = np.zeros((len(units), 1, 49))
    for count, i in enumerate(meanSpikeReshaped):
        # row1 of 7x7 grid when reshaped
        i[:, 0:6] = meanSpike[count][0:6]
        i[:, 6] = meanSpike[count][42]
        # row2
        i[:, 7:13] = meanSpike[count][6:12]
        i[:, 13] = meanSpike[count][43]
        # row3
        i[:, 14:20] = meanSpike[count][12:18]
        i[:, 20] = meanSpike[count][44]
        # row4
        i[:, 21:27] = meanSpike[count][18:24]
        i[:, 27] = meanSpike[count][45]
        # row5
        i[:, 28:34] = meanSpike[count][24:30]
        i[:, 34] = meanSpike[count][46]
        # row6
        i[:, 35:41] = meanSpike[count][30:36]
        i[:, 41] = meanSpike[count][47]
        # row7
        i[:, 42:48] = meanSpike[count][36:42]
        i[:, 48] = meanSpike[count][48]

    # Unit's preferred direction based off of
    # gaussian fit of MTNC stimuli:
    angleMat = np.arange(180, 900, 60)
    unitGaussMean = np.zeros(len(units))
    unitGaussSig = np.zeros(len(units))
    for unitCount, unit in enumerate(units):
        loc0Resp = meanSpike[unitCount][36:42]
        loc1Resp = meanSpike[unitCount][42:48]
        combResp = (loc0Resp + loc1Resp) / 2
        extRespMat = np.concatenate((combResp[3:], combResp, combResp[:3]), axis=0)
        maxIndex = np.where(combResp == np.max(combResp))[0][0] + 3
        x = angleMat[maxIndex-3:maxIndex+4]
        y = extRespMat[maxIndex-3:maxIndex+4]
        params = gaussFit(x, y)
        unitGaussMean[unitCount] = params[2] % 360
        unitGaussSig[unitCount] = params[3]

    # scipy curvefit Normalization parameters
    unitAlphaIndex = np.zeros(len(units))
    unitGaussFit = []
    unitPairedNormFit = []
    unitGaussR2 = np.zeros(len(units))
    unitPairedEV = np.zeros(len(units))
    unitPairedNormR2 = np.zeros(len(units))
    unitNormFitEstimate = [[] for i in range(len(units))]
    scaledLoc = np.zeros(len(units))
    scalar = np.zeros(len(units))
    for unitCount, unit in enumerate(units):
        b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS
        bReIndex = np.zeros((7, 7))

        # find direction tested that is closest to the pref dir
        # and reindex around this so that it is in the middle of the grid
        prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
        nullDirIndex = np.where(dirArray == nullDir)[0][0]
        reIndex = (np.array([0,1,2,3,4,5])+nullDirIndex) % 6

        tempMain = b[:6, :6][:, reIndex]
        tempMain = tempMain[:6, :6][reIndex, :]
        temp0Blank = b[:6, 6][reIndex]
        temp1Blank = b[6, :6][reIndex]
        bReIndex[:6, :6] = tempMain
        bReIndex[:6, 6] = temp0Blank
        bReIndex[6, :6] = temp1Blank
        bReIndex[6, 6] = b[6, 6]

        # fixed (independent) variables - matrix of corresponding stim Indexes
        stimMat = np.zeros((7, 7))
        stimMat[:6, :6] = np.arange(36).reshape(6, 6)
        stimMat[6, :6] = np.arange(36, 42)
        stimMat[:, 6] = np.arange(42, 49)

        # reshape fixed variables to match dependent variables
        stimMatReIndex = np.zeros((7, 7))
        tempMain = stimMat[:6, :6][:, reIndex]
        tempMain = np.squeeze(tempMain[:6, :6][reIndex, :])
        temp0Blank = stimMat[:6, 6][reIndex]
        temp1Blank = stimMat[6, :6][reIndex]
        stimMatReIndex[:6, :6] = tempMain
        stimMatReIndex[:6, 6] = temp0Blank
        stimMatReIndex[6, :6] = temp1Blank
        stimMatReIndex[6, 6] = stimMat[6, 6]

        # Generic Normalization (L1+L2)/(al1+al2+sig) w.o scalar
        # fits L0-L6 for both locations separately loc0 and loc1 will
        # have different L0-L6 values
        resp = b.reshape(49)
        fixedVals = fixedValsForGenericNorm(stimMat.reshape(49), stimIndexDict)
        pOpt, pCov = curve_fit(genericNormNoScalar, fixedVals, resp.squeeze(),
                               bounds=((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                       (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
                                        np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
                                        3, 3, .10)))
        print(unit, pOpt)
        y_pred = genericNormNoScalar(fixedVals, *pOpt)
        r2 = r2_score(resp.squeeze(), y_pred)
        print(r2)

        totR2.append(r2)
        # unitPairedEV[unitCount] = EV
        unitPairedNormR2[unitCount] = r2
        # scaledLoc[unitCount] = sLoc
        unitPairedNormFit.append(pOpt)
        unitNormFitEstimate[unitCount] = pOpt

    # generate paired stimulus correlation, selectivity index,
    # suppression, and NMI index for each gabor pair for each pair
    # of neurons (similar to Bram fig 2c)
    pairPairedCorr = []
    pairSelectivityIndex = []
    pairNonPrefSuppIndex = []
    pairNMIIndex = []
    pairSimpleNMIIndex = []

    pairSingleCorr = []
    pairSingleGaborNMI = []
    pairSingleGaborSelectivity = []

    tempMat = stats.zscore(spikeCountMat[:, :blocksDone, :], axis=1, nan_policy='omit')
    stimIndexMat = np.arange(36).reshape(6, 6)
    upperTriangle = upperTriMasking(stimIndexMat)

    for pairCount, pair in enumerate(combs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]

        neuronPairAvgSelectivity = []
        neuronPairAvgNMIloc0 = []
        neuronPairAvgNMIloc1 = []
        # indices and correlations for paired Gabor stimuli
        for i in range(36):
            # correlation for that Gabor pair b/w 2 units
            # pairStimCorr = stats.pearsonr(tempMat[n1, :, i], tempMat[n2, :, i])
            pairStimCorr = stats.pearsonr(spikeCountMat[n1, :blocksDone, i],
                                          spikeCountMat[n2, :blocksDone, i])

            # extract directions of the gabor pair
            loc0Dir = stimIndexDict[i][0]['direction']
            loc1Dir = stimIndexDict[i][1]['direction']

            # n1 selectivity and suppression index
            loc0Resp = unitNormFitEstimate[n1][np.where(dirArray == loc0Dir)[0][0]]
            loc1Resp = unitNormFitEstimate[n1][np.where(dirArray == loc1Dir)[0][0] + 6]
            if loc0Resp > loc1Resp:
                n1Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
                n1NonPrefSupp = unitNormFitEstimate[n1][13] / (
                    unitNormFitEstimate[n1][13] + unitNormFitEstimate[n1][12])
                n1PrefResp = 0
            else:
                n1Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
                n1NonPrefSupp = unitNormFitEstimate[n1][12] / (
                    unitNormFitEstimate[n1][13] + unitNormFitEstimate[n1][12])
                n1PrefResp = 1

            # n2 selectivity and suppression index
            loc0Resp = unitNormFitEstimate[n2][np.where(dirArray == loc0Dir)[0][0]]
            loc1Resp = unitNormFitEstimate[n2][np.where(dirArray == loc1Dir)[0][0] + 6]
            if loc0Resp > loc1Resp:
                n2Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
                n2NonPrefSupp = unitNormFitEstimate[n2][13] / (
                    unitNormFitEstimate[n2][13] + unitNormFitEstimate[n2][12])
                n2PrefResp = 0
            else:
                n2Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
                n2NonPrefSupp = unitNormFitEstimate[n2][12] / (
                    unitNormFitEstimate[n2][13] + unitNormFitEstimate[n2][12])
                n2PrefResp = 1

            # n1 NMI
            b = meanSpikeReshaped[n1].reshape(7, 7) * 1000/trueStimDurMS
            pairedMat = b[:6, :6].reshape(36)
            loc0Single = b[6, :6]
            loc1Single = b[:6, 6]
            loc0SingleResp = loc0Single[np.where(dirArray == loc0Dir)[0]][0]
            loc1SingleResp = loc1Single[np.where(dirArray == loc1Dir)[0]][0]
            n1NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[i]) / (
                  (loc0SingleResp + loc1SingleResp) + pairedMat[i])
            n1SimpleNMI = pairedMat[i] / (loc0SingleResp + loc1SingleResp)

            # n2 NMI
            b = meanSpikeReshaped[n2].reshape(7, 7) * 1000/trueStimDurMS
            pairedMat = b[:6, :6].reshape(36)
            loc0Single = b[6, :6]
            loc1Single = b[:6, 6]
            loc0SingleResp = loc0Single[np.where(dirArray == loc0Dir)[0]][0]
            loc1SingleResp = loc1Single[np.where(dirArray == loc1Dir)[0]][0]
            n2NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[i]) / (
                  (loc0SingleResp + loc1SingleResp) + pairedMat[i])
            n2simpleNMI = pairedMat[i] / (loc0SingleResp + loc1SingleResp)

            # pair selectivity, suppression, and NMI index
            pairSelectivity = np.sqrt(n1Selectivity * n2Selectivity)
            if n1PrefResp != n2PrefResp:
                pairSelectivity = -pairSelectivity
            pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)
            pairNMI = (n1NMI + n2NMI) / 2
            pairSimpleNMI = np.sqrt(n1SimpleNMI * n2simpleNMI)

            pairPairedCorr.append(pairStimCorr[0])
            pairSelectivityIndex.append(pairSelectivity)
            pairNonPrefSuppIndex.append(pairSuppression)
            pairNMIIndex.append(pairNMI)
            pairSimpleNMIIndex.append(pairSimpleNMI)

            # appending pair Selectivity and pair NMI to get average
            # for the single Gabor figure
            neuronPairAvgSelectivity.append(pairSelectivity)
            if i in upperTriangle:
                neuronPairAvgNMIloc0.append(pairNMI)
            else:
                neuronPairAvgNMIloc1.append(pairNMI)

        # indices and correlations for single Gabor b/w pair of neuron
        for j in range(36, 48):
            pairSingCorr = stats.pearsonr(spikeCountMat[n1, :blocksDone, j],
                                          spikeCountMat[n2, :blocksDone, j])

            # extract location and direction of tested gabor
            loc0Dir = stimIndexDict[j][0]['direction']
            loc0Con = stimIndexDict[j][0]['contrast']
            loc1Dir = stimIndexDict[j][1]['direction']
            loc1Con = stimIndexDict[j][1]['contrast']
            if loc0Con == 1:
                avgNMI = np.mean(neuronPairAvgNMIloc1)
            else:
                avgNMI = np.mean(neuronPairAvgNMIloc0)

            pairSingleGaborNMI.append(avgNMI)
            pairSingleGaborSelectivity.append(np.mean(neuronPairAvgSelectivity))
            pairSingleCorr.append(pairSingCorr[0])

    pairPairedCorr = np.array(pairPairedCorr)
    pairSelectivityIndex = np.array(pairSelectivityIndex)
    pairNonPrefSuppIndex = np.array(pairNonPrefSuppIndex)
    pairNMIIndex = np.array(pairNMIIndex)
    pairSimpleNMIIndex = np.array(pairSimpleNMIIndex)
    pairSingleCorr = np.array(pairSingleCorr)
    pairSingleGaborNMI = np.array(pairSingleGaborNMI)
    pairSingleGaborSelectivity = np.array(pairSingleGaborSelectivity)

    ## compile similarity scores and correlation matrix into one matrix
    gaborPairCompiledArr = np.array([pairPairedCorr,
                                     pairSelectivityIndex,
                                     pairNonPrefSuppIndex,
                                     pairNMIIndex,
                                     pairSimpleNMIIndex])
    gaborSingleCompiledArr = np.array([pairSingleCorr,
                                       pairSingleGaborNMI,
                                       pairSingleGaborSelectivity])
    np.save(f'../../gaborSingleCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborSingleCompiledArr)
    np.save(f'../../gaborPairCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborPairCompiledArr)
    os.chdir('../../../Python Code')

    totUnits.append(len(units))

print(time.time()-t0)