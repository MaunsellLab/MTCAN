'''
MTN2 Analysis Script

to do:
convert heatmap to spikes/sec, it's at spikes/stimDurMS
'''

import seaborn as sns
import numpy as np
import numpy.ma as ma
from usefulFns import *
from itertools import combinations
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time

####### START HERE ######
# for loop to run through all files
t0 = time.time()
totUnits = []
totR2 = []
fileList = ['221010', '221108', '221110', '221115', '221117', '221124',
            '221128', '221010', '221208', '221229']
# fileList = ['221115', '221117', '221124',
#             '221128', '221010', '221208']

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


    '''
    ## Filter Units
    # inclusion criteria:
    # 1. pref response > 2x null response
    # 2. pref response > baseline (need to incorporate this)
    filterUnits = []
    for unitCount, unit in enumerate(units):
        b = meanSpikeReshaped[unitCount].reshape(7,7)
        bSmooth = gaussian_filter(b, sigma=1)
        prefDir, nullDir = unitPrefNullDir(bSmooth)
    
        orientList = [(prefDir,'blank'),(nullDir,'blank'),
                      ('blank',prefDir),('blank',nullDir)] 
        respList = []
        for i in orientList:
            loc0Dir, loc1Dir = i
            loc0Con = highContrast
            loc1Con = highContrast
            if loc0Dir == 'blank':
                loc0Dir = 0
                loc0Con = 0
            if loc1Dir == 'blank':
                loc1Dir = 0
                loc1Con = 0
            sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                       (stimIndexDF['loc0 Contrast'] == loc0Con) &
                                       (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                       (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
            unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit) 
                        & (spikeCountDF['stimIndex'] == sIndex)]
            unitDF = unitDF.iloc[:blocksDone].mean()
            meanResp = unitDF[3] * 1000/trueStimDurMS
            respList.append(meanResp)
        
        if respList[0] > 2*respList[1] or respList[2] > 2*respList[3]:
            filterUnits.append(unit)
    '''

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

    # Z-scored Correlations
    # filtered units test
    # combs = [i for i in combinations(filterUnits, 2)]

    # unfiltered units
    combs = [i for i in combinations(units, 2)]
    corrMat = np.zeros(len(combs))

    # z-scored spikeCountMat
    zSpikeCountMat = stats.zscore(spikeCountMat[:, :blocksDone, :], axis=1, nan_policy='omit')
    zSpikeCountMat = np.nan_to_num(zSpikeCountMat)
    zSpikeCountMat = np.reshape(zSpikeCountMat, (len(units), blocksDone*49))
    for count, i in enumerate(combs):
        n1 = np.where(units == i[0])[0][0]
        n2 = np.where(units == i[1])[0][0]
        pairCorr = stats.pearsonr(zSpikeCountMat[n1], zSpikeCountMat[n2])
        corrMat[count] = pairCorr[0]

    popCorr = np.mean(corrMat)

    # distance of neuron pair (based on channel number)
    pairChannelDistance = []
    for count, i in enumerate(combs):
        n1Chan = int(unitsChannel[np.where(units == i[0])[0][0]])
        n2Chan = int(unitsChannel[np.where(units == i[1])[0][0]])
        channelDiff = abs(n1Chan-n2Chan)
        pairChannelDistance.append(channelDiff)

    pairChannelDistance = np.array(pairChannelDistance)

    ## RF location tuning similarity b/w neurons Bhattacharyya Distance 2D
    RFLocMat = np.load('../RFLoc Tuning/unitsRFLocMat.npy')
    pairLocSimScore = np.zeros((len(combs),1))
    for pairCount, pair in enumerate(combs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]

        # a = np.flip(RFLocMat[n1], axis=0)
        # b = np.flip(RFLocMat[n2], axis=0)
        a = RFLocMat[n1]
        b = RFLocMat[n2]
        m1, cov1, p = gauss2dParams(a)
        m2, cov2, p2 = gauss2dParams(b)
        BC = bhattCoef2D(m1,m2,cov1,cov2)
        pairLocSimScore[pairCount] = BC

    # direction tuning similarity b/w neurons using MTNC stim: Bhattacharyya Distance
    pairSimPrefDir = np.zeros((len(combs), 1))
    pairSimScore = np.zeros((len(combs), 1))
    for pairCount, pair in enumerate(combs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]
        m1 = unitGaussMean[n1]
        v1 = unitGaussSig[n1]**2
        m2 = unitGaussMean[n2]
        v2 = unitGaussSig[n2]**2

        # similarity of pref dirs only
        if abs(m1 - m2) > 180:
            if m1 > m2:
                m1 = m2 - (360 - (m1 - m2))
            else:
                m2 = m1 - (360 - (m2 - m1))
        pairSimPrefDir[pairCount] = 1 - (abs(m1 - m2) / 180)
        # bhattacharyya similarity score
        BC = bhattCoef(m1, m2, v1, v2)
        pairSimScore[pairCount] = BC

    '''
    ## direction tuning similarity b/w neurons Bhattacharyya Distance
    dirTuningMat = np.load('../Direction Tuning/unitsDirTuningMat.npy') #load from directions folder
    # unitsBaselineMat = np.load('../Direction Tuning/unitsBaselineMat.npy')
    extTunMat = np.concatenate((dirTuningMat[:,3:], dirTuningMat[:,:], 
                                dirTuningMat[:,:3]), axis=1)
    angleMat = np.arange(180,900,60)
    # unitsPrefDirMat = np.zeros(len(units))
    combs = [i for i in combinations(units, 2)]
    pairSimPrefDir = np.zeros((len(combs),1))
    pairSimScore = np.zeros((len(combs),1))
    for pairCount, pair in enumerate(combs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]
        n1Max = int(np.where(dirTuningMat[n1] == np.max(dirTuningMat[n1]))[0][0] + 3)
        n1X = angleMat[n1Max-3:n1Max+4]
        # n1Y = extTunMat[n1][n1Max-3:n1Max+4]
        n1Y = extTunMat[n1][n1Max-3:n1Max+4]/max(extTunMat[n1][n1Max-3:n1Max+4])
        n1XFull = np.linspace(n1X[0],n1X[-1],1000)
        params = gaussFit(n1X, n1Y)
        # n1YFull = gauss(n1XFull, *params)
        m1 = params[2] # mean neuron 1
        v1 = params[3]**2 # var neuron 1
        n1TrueMean = m1 % 360
        
        n2Max = int(np.where(dirTuningMat[n2] == np.max(dirTuningMat[n2]))[0][0] + 3)
        n2X = angleMat[n2Max-3:n2Max+4]
        # n2Y = extTunMat[n2][n2Max-3:n2Max+4]
        n2Y = extTunMat[n2][n2Max-3:n2Max+4]/max(extTunMat[n2][n2Max-3:n2Max+4])
        n2XFull = np.linspace(n2X[0], n2X[-1],1000)
        params = gaussFit(n2X, n2Y)
        # n2YFull = gauss(n2XFull, *params)
        m2 = params[2]
        v2 = params[3]**2
        n2TrueMean = m2 % 360
    
        m1 = m1%360
        m2 = m2%360
        if abs(m1-m2) > 180:
            if m1 > m2:
                m1 = m2-(360-(m1-m2))
            else:
                m2 = m1-(360-(m2-m1))
    
        # similarity of pref dirs only
        pairSimPrefDir[pairCount] = 1 - (abs(m1-m2)/180)
        # bhattacharyya similarity score 
        BC = bhattCoef(m1, m2, v1, v2)
        pairSimScore[pairCount] = BC
        # add unit's pref dir to unitPrefDir mat
        # unitsPrefDirMat[n1] = n1TrueMean
        # unitsPrefDirMat[n2] = n2TrueMean
    '''

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
        resp = b.reshape(49)[:-1]
        fixedVals = fixedValsForGenericNorm(stimMat.reshape(49)[:-1], stimIndexDict)
        pOpt, pCov = curve_fit(genericNormNoScalar, fixedVals, resp.squeeze(),
                               bounds=((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                       (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
                                        np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
                                        1, 1)))
        print(unit, pOpt)
        y_pred = genericNormNoScalar(fixedVals, *pOpt)
        r2 = r2_score(resp.squeeze(), y_pred)
        print(r2)

        # Generic Normalization (L1+L2)/(al1+al0+sig) w Scalar

        # resp = b.reshape(49)
        # fixedVals = fixedValsForGenericNorm(stimMat.reshape(49), stimIndexDict)
        # if max(bReIndex[:6, 6]) > max(bReIndex[6, :6]):
        #     pOpt, pCov = curve_fit(genericNorm1, fixedVals, resp.squeeze(), bounds=(
        #         (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        #         (np.inf, np.inf, np.inf, np.inf, np.inf,
        #          np.inf, 1, 5, 5, np.inf)))
        #     print(unit, pOpt)
        #     y_pred = genericNorm1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     sLoc = 0
        # else:
        #     print('using norm func 0')
        #     pOpt, pCov = curve_fit(genericNorm0, fixedVals, resp.squeeze(), bounds=(
        #         (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        #         (np.inf, np.inf, np.inf, np.inf, np.inf,
        #          np.inf, 1, 5, 5, np.inf)))
        #     print(unit, pOpt)
        #     y_pred = genericNorm0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     sLoc = 1

        # for fitting EMS normalization generic (without Gauss tuning)

        # resp = bReIndex.reshape(49)[:-1]
        # fixedVals = fixedValsForEMSGen(stimMatReIndex.reshape(49)[:-1], stimIndexDict)
        #
        # if max(bReIndex[:6, 6])-min(bReIndex[:6, 6]) > max(bReIndex[6, :6])-min(bReIndex[6, :6]):
        #     pOpt, pCov = curve_fit(emsNormGen1, fixedVals, resp, bounds=(
        #         (0, 0, 0, 0, 0, 0, 0, 0, 0),
        #         (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
        #     print(unit, pOpt)
        #     y_pred = emsNormGen1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     sLoc = 0
        # else:
        #     print('using norm func 0')
        #     pOpt, pCov = curve_fit(emsNormGen0, fixedVals, resp, bounds=(
        #         (0, 0, 0, 0, 0, 0, 0, 0, 0),
        #         (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
        #     print(unit, pOpt)
        #     y_pred = emsNormGen0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     sLoc = 1

        # # EMS normalization two steps with Gauss tuning curve

        # # step 1: fitting gaussian function to single presentations
        # # to extract gaussian fit and scaling factor (scalar on loc0)
        # # if scalar >1 then loc0 has stronger resp
        # resp = np.concatenate((bReIndex[6, :6], bReIndex[:6, 6]), axis=0)
        # singleStim = np.concatenate((stimMatReIndex[6, :6], stimMatReIndex[:6, 6]),
        #                             axis=0)
        # fixedVals = fixedValsForCurveFit(prefDir, singleStim, stimIndexDict)
        # if max(bReIndex[6, :6]) > max(bReIndex[:6, 6]):
        #     pOptGauss, pCov = curve_fit(gaussNormFunc0, fixedVals, resp.squeeze(),
        #                                 bounds=((0, 0, 0, 0, 0),
        #                                 (np.inf, max(resp)*1.5, 360, 360, 1)))
        #     y_pred = gaussNormFunc0(fixedVals, *pOptGauss)
        #     print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')
        #     sLoc = 1
        # else:
        #     pOptGauss, pCov = curve_fit(gaussNormFunc1, fixedVals, resp.squeeze(),
        #                                 bounds=((0, 0, 0, 0, 0),
        #                                 (np.inf, max(resp)*1.5, 360, 360, 1)))
        #     y_pred = gaussNormFunc1(fixedVals, *pOptGauss)
        #     print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')
        #     print('scaled loc=0')
        #     sLoc = 0
        # unitGaussR2[unitCount] = r2_score(resp.squeeze(), y_pred)
        # unitGaussFit.append(pOptGauss)
        # scalar[unitCount] = pOptGauss[4]
        #
        # # step 2: fitting  EMS normalization eqn with fixed gauss tuning curve
        # # to paired stimuli
        # # only two free parameters (a0, a1)
        # resp = bReIndex[:6, :6].reshape(36)
        # pairedStim = stimMatReIndex[:6, :6].reshape(36)
        # pairedStim = pairedStim.astype(int)
        # fixedVals = fixedValsCurveFitForPairedStim(prefDir, pairedStim, stimIndexDict, pOptGauss)
        # if sLoc == 1:
        #     pOpt, pCov = curve_fit(EMSGaussNormStepFunc0, fixedVals, resp.squeeze(),
        #                            bounds=((-1, -1), (5, 5)))
        #     y_pred = EMSGaussNormStepFunc0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)
        #     print('using EMS Func0, sLoc==1')
        # else:
        #     pOpt, pCov = curve_fit(EMSGaussNormStepFunc1, fixedVals, resp.squeeze(),
        #                            bounds=((-1, -1), (5, 5)))
        #     y_pred = EMSGaussNormStepFunc1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)


        # # Fitting two step EMS Normalization Equation with L0-L6
        # # instead of guassian

        # # step 1, fitting single stimulus with L0-L6 for direction tuning
        # resp = np.concatenate((bReIndex[6, :6], bReIndex[:6, 6]), axis=0)
        # singleStim = np.concatenate((stimMatReIndex[6, :6], stimMatReIndex[:6, 6]),
        #                             axis=0)
        # fixedVals = fixedValsForEMSGen(singleStim, stimIndexDict)
        # if max(bReIndex[6, :6]) > max(bReIndex[:6, 6]):
        #     pOpt, pCov = curve_fit(emsSingleStim0, fixedVals, resp.squeeze(),
        #                            bounds=((0, 0, 0, 0, 0, 0, 0),
        #                                    (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
        #     y_pred = emsSingleStim0(fixedVals, *pOpt)
        #     print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')
        #     sLoc = 1
        # else:
        #     pOpt, pCov = curve_fit(emsSingleStim1, fixedVals, resp.squeeze(),
        #                            bounds=((0, 0, 0, 0, 0, 0, 0),
        #                                    (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
        #     y_pred = emsSingleStim1(fixedVals, *pOpt)
        #     print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')
        #     print('scaled loc=0')
        #     sLoc = 0
        #
        # # step 2, fitting paired stimulus to EMS normalization eqn with
        # # fixed L0-L6 and S values. Only two free parameters (a0,a1)
        # resp = bReIndex[:6, :6].reshape(36)
        # pairedStim = stimMatReIndex[:6, :6].reshape(36)
        # pairedStim = pairedStim.astype(int)
        # fixedVals = fixedValsForPairedStimL0L6(pairedStim, stimIndexDict, pOpt)
        # if sLoc == 1:
        #     pOpt, pCov = curve_fit(EMSGenNormPaired0, fixedVals, resp.squeeze())
        #     y_pred = EMSGenNormPaired0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)
        #     print('using EMS Func0, sLoc==1')
        # else:
        #     pOpt, pCov = curve_fit(EMSGenNormPaired1, fixedVals, resp.squeeze())
        #     y_pred = EMSGenNormPaired1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)

        totR2.append(r2)
        # unitPairedEV[unitCount] = EV
        unitPairedNormR2[unitCount] = r2
        # scaledLoc[unitCount] = sLoc
        unitPairedNormFit.append(pOpt)
        unitNormFitEstimate[unitCount] = pOpt
        unitAlphaIndex[unitCount] = (pOpt[0] + pOpt[1]) / 2

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

    # cm = plt.cm.get_cmap('seismic')
    # y = pairSelectivityIndex # pair selectivity
    # x = pairNMIIndex # pair non-pref alpha index
    # z = pairCorr # noise correlations
    # sc = plt.scatter(x, y, c=z, s=35, cmap=cm, vmin=-0.6, vmax=0.6)
    # # vmin=-0.6, vmax=0.6
    # plt.colorbar(sc)
    # plt.xlabel('non-preferred suppression')
    # plt.ylabel('pair selectivity')
    # plt.show()

    # generate pair's combined Norm Value
    pairAlphaMulti = np.zeros((len(combs)))
    pairR2 = np.zeros(len(combs))
    for pairCount, pair in enumerate(combs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]

        pairR2[pairCount] = np.sqrt(unitPairedNormR2[n1] * unitPairedNormR2[n2])

        pairAlphaMulti[pairCount] = (unitAlphaIndex[n1] + unitAlphaIndex[n2]) / 2

    pairCombinedSimScore = np.squeeze((pairLocSimScore + pairSimScore) / 2)
    pairCombinedSimScoreMulti = np.squeeze((pairLocSimScore * pairSimScore))
    pairCombinedSimScorePrefDir = np.squeeze((pairLocSimScore + pairSimPrefDir) / 2)
    pairCombinedPrefDirAlpha = (pairSimPrefDir.flatten() + pairAlphaMulti) / 2

    ## Create subset of Pandas DF for P, N, P+N conditions for all units
    pnContrastPlotDF = pd.DataFrame(columns=['unit', 'stimIndex', 'stimCount',
                                    'stimSpikes', 'contrast', 'prefNullStr'])
    for unitCount, unit in enumerate(units):
        # find direction tested that is closest to unit's preferred and null direction
        prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
        orientCount = 0
        orientList = ['pref+pref', 'pref+null', 'null+pref', 'null+null',
                     'pref+blank', 'null+blank','blank+pref', 'blank+null']
        for j in [(prefDir,prefDir),(prefDir,nullDir),(nullDir,prefDir),(nullDir,nullDir),]:
            loc0Dir, loc1Dir = j
            sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                       (stimIndexDF['loc0 Contrast'] == highContrast) &
                                       (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                       (stimIndexDF['loc1 Contrast'] == highContrast)][0]
            unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                        & (spikeCountDF['stimIndex'] == sIndex)]
            unitDF = unitDF.iloc[:blocksDone]
            unitDF['contrast'] = highContrast
            unitDF['prefNullStr'] = orientList[orientCount]
            pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
            orientCount += 1

        #Loc 0 Pref/Null only
        for i in [(prefDir, zeroDir), (nullDir, zeroDir)]:
            loc0Dir, loc1Dir = i
            sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                       (stimIndexDF['loc0 Contrast'] == highContrast) &
                                       (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                       (stimIndexDF['loc1 Contrast'] == zeroContrast)][0]
            unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                        & (spikeCountDF['stimIndex'] == sIndex)]
            unitDF = unitDF.iloc[:blocksDone]
            unitDF['contrast'] = zeroContrast
            unitDF['prefNullStr'] = orientList[orientCount]
            pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
            orientCount += 1

        #loc 1 Pref/Null only
        for x in [(zeroDir,prefDir), (zeroDir,nullDir)]:
            loc0Dir, loc1Dir = x
            sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                       (stimIndexDF['loc0 Contrast'] == zeroContrast) &
                                       (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                       (stimIndexDF['loc1 Contrast'] == highContrast)][0]
            unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                        & (spikeCountDF['stimIndex'] == sIndex)]
            unitDF = unitDF.iloc[:blocksDone]
            unitDF['contrast'] = zeroContrast
            unitDF['prefNullStr'] = orientList[orientCount]
            pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
            orientCount += 1

    ## generate within unit NMI score for P+N and N+P condition only
    unitWithinNMI = np.zeros(len(units))
    unitNMI0 = np.zeros(len(units))
    unitNMI1 = np.zeros(len(units))
    for unitCount, unit in enumerate(units):
        unitDF = pnContrastPlotDF.loc[pnContrastPlotDF['unit'] == unit]
        #spikeCounts in spikes/sec
        unitDF['stimSpikes'] = unitDF['stimSpikes'] * 1000/trueStimDurMS
        singleResp = np.zeros(4)
        doubleResp = np.zeros(4)
        for count, indx in enumerate(['pref+blank', 'null+blank', 'blank+pref', 'blank+null']):
            mean = unitDF.loc[(unitDF['prefNullStr'] == indx)].mean()[3]
            singleResp[count] = mean
        for count, indx in enumerate(['pref+pref', 'pref+null', 'null+pref', 'null+null']):
            mean = unitDF.loc[(unitDF['prefNullStr'] == indx)].mean()[3]
            doubleResp[count] = mean

        # NMI
        # avgP = (singleResp[0] + singleResp[2]) / 2
        # avgN = (singleResp[1] + singleResp[3]) / 2
        # avgPN = (doubleResp[1] + doubleResp[2]) / 2
        # unitWithinNMI[unitCount] = avgPN / avgP + avgN + avgPN

        #NMI loc 0
        NMI0 = (doubleResp[0] + doubleResp[3] - doubleResp[1]) / \
               (doubleResp[0] + doubleResp[3] + doubleResp[1])
        #NMI loc 1
        NMI1 = (doubleResp[0] + doubleResp[3] - doubleResp[2]) / \
               (doubleResp[0] + doubleResp[3] + doubleResp[2])
        unitNMI0[unitCount] = NMI0
        unitNMI1[unitCount] = NMI1
        unitWithinNMI[unitCount] = NMI0 + NMI1

        # NMI loc 0
        # NMI0 = (singleResp[0] - singleResp[3] - doubleResp[1] - singleResp[3]) / (
        #         singleResp[0] - singleResp[3] + doubleResp[1] - singleResp[3])
        # NMI1 = (singleResp[1] - singleResp[2] - doubleResp[2] - singleResp[2]) / (
        #         singleResp[1] - singleResp[2] + doubleResp[2] - singleResp[2])
        # unitWithinNMI[unitCount] = (NMI0 + NMI1) / 2

        # # #NMI loc 0
        # NMI0 = (doubleResp[0] + doubleResp[3] - doubleResp[1]) / (
        #         doubleResp[0] + doubleResp[3] + doubleResp[1])
        # NMI1 = (doubleResp[0] + doubleResp[3] - doubleResp[2]) / (
        #         doubleResp[0] + doubleResp[3] + doubleResp[2])
        # unitWithinNMI[unitCount] = (NMI0 + NMI1) / 2

        # #NMI loc 0
        # NMI0 = doubleResp[1] / doubleResp[0]
        # NMI1 = doubleResp[2] / doubleResp[0]
        # unitWithinNMI[unitCount] = min([NMI0, NMI1])

    # NMI P+N/N+P across pairs
    pairNMI = np.zeros(len(combs))
    for pairCount, pair in enumerate(combs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]

        pairNMI[pairCount] = (unitNMI0[n1] + unitNMI0[n2]) - (unitNMI1[n1] + unitNMI1[n2]) / (
                              unitNMI0[n1] + unitNMI0[n2]) + (unitNMI1[n2] + unitNMI1[n2])


    # generate within unit NMI for every paired stimulus
    unitAllNMI = np.zeros((len(units), 36))
    for unitCount, unit in enumerate(units):
        b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS
        pairedMat = b[:6,:6].reshape(36)
        loc0Single = b[6,:6]
        loc1Single = b[:6,6]
        directionSet = np.array([0, 60, 120, 180, 240, 300])

        for i in range(36):
            loc0Dir = stimIndexDict[i][0]['direction']
            loc1Dir = stimIndexDict[i][1]['direction']
            loc0SingleResp = loc0Single[np.where(directionSet == loc0Dir)[0]][0]
            loc1SingleResp = loc1Single[np.where(directionSet == loc1Dir)[0]][0]
            # NMI = ((loc0SingleResp + loc1SingleResp) - pairedMat[i]) / (
            #       (loc0SingleResp + loc1SingleResp) + pairedMat[i])
            NMI = pairedMat[i] / (loc0SingleResp + loc1SingleResp)
            unitAllNMI[unitCount][i] = NMI

    pairAllNMI = np.zeros(len(combs))
    for pairCount, pair in enumerate(combs):
        n1 = np.where(units == pair[0])[0][0]
        n2 = np.where(units == pair[1])[0][0]

        m1 = np.median(unitAllNMI[n1])
        v1 = np.std(unitAllNMI[n1]) ** 2
        m2 = np.median(unitAllNMI[n2])
        v2 = np.std(unitAllNMI[n2]) ** 2
        BC = bhattCoef(m1, m2, v1, v2)
        # pairAllNMI[pairCount] = BC
        pairAllNMI[pairCount] = (m1 * m2)

    # cm = plt.cm.get_cmap('seismic')
    # x = np.squeeze((pairLocSimScore + pairSimPrefDir) / 2) # RFoverlap and same dir pref
    # y = pairAllNMI # pair NMI mean
    # z = corrMat # noise correlations
    # sc = plt.scatter(x, y, c=z, s=35, cmap=cm, vmin=-0.5, vmax=0.5)
    # plt.colorbar(sc)
    # plt.show()
    # print(stats.pearsonr(pairAllNMI, corrMat))
    # print(stats.pearsonr(pairSimPrefDir.flatten(), corrMat))
    # print(stats.pearsonr((pairSimPrefDir.flatten() + pairLocSimScore.flatten())/2, corrMat))
    # print(stats.pearsonr((pairSimPrefDir.flatten() + pairAllNMI + pairLocSimScore.flatten())/3, corrMat))
    # plt.scatter((pairSimPrefDir.flatten() + pairAllNMI + pairLocSimScore.flatten())/3, corrMat); plt.show()


    ## compile similarity scores and correlation matrix into one matrix
    compiledArr = np.array([corrMat,
                            np.squeeze(pairSimScore),
                            np.squeeze(pairSimPrefDir),
                            np.squeeze(pairLocSimScore),
                            pairCombinedSimScore,
                            pairCombinedSimScoreMulti,
                            pairCombinedSimScorePrefDir,
                            pairAlphaMulti,
                            pairR2,
                            pairNMI,
                            pairCombinedPrefDirAlpha,
                            pairAllNMI])
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
    np.save(f'../../corrMaster/pairCorrelationsAndSimScores{seshDate}',
            compiledArr)
    np.save(f'../../gaborPairCorrMaster/pairCorrelationsAndIndices{seshDate}',
            gaborPairCompiledArr)
    os.chdir('../../../Python Code')

    totUnits.append(len(units))

print(time.time()-t0)

## SUPERPLOT OF PSTH, Normalization (pref+null+blank) heatmap, and bar plot
## Create subset of Pandas DF for P, N, P+N conditions for all units
pnContrastPlotDF = pd.DataFrame(columns=['unit', 'stimIndex', 'stimCount',
                                'stimSpikes', 'contrast', 'prefNullStr'])
for unitCount, unit in enumerate(units):
    # find direction tested that is closest to unit's preferred and null direction
    prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
    orientCount = 0
    orientList = ['pref+pref', 'pref+null', 'null+pref', 'null+null',
                 'pref+blank', 'null+blank','blank+pref', 'blank+null']
    for j in [(prefDir,prefDir),(prefDir,nullDir),(nullDir,prefDir),(nullDir,nullDir),]:
        loc0Dir, loc1Dir = j
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == highContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == highContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                    & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = highContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1

    #Loc 0 Pref/Null only
    for i in [(prefDir, zeroDir), (nullDir, zeroDir)]:
        loc0Dir, loc1Dir = i
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == highContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == zeroContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                    & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = zeroContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1

    #loc 1 Pref/Null only
    for x in [(zeroDir,prefDir), (zeroDir,nullDir)]:
        loc0Dir, loc1Dir = x
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == zeroContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == highContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                    & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = zeroContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1

for unitCount, unit in enumerate(units):
    b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS
    bSmooth = gaussian_filter(b, sigma=1)
    prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])

    ## figure 
    # fig = plt.figure(constrained_layout=True)
    fig = plt.figure()
    fig.set_size_inches(17, 7)

    gs0 = gridspec.GridSpec(1, 5)

    heights = [1, 1, 1]
    widths = [4, 4, 4]
    gs00 = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs0[0],
                                            width_ratios=widths,
                                            height_ratios=heights)
    gs01 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[1:3])
    gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[3])
    gs03 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[4])

    # PSTHs Plot
    locDirList = [(nullDir, nullDir), (prefDir, nullDir), ('blank', nullDir),
                  (nullDir, prefDir), (prefDir, prefDir), ('blank', prefDir),
                  (nullDir, 'blank'), (prefDir, 'blank'), ('blank', 'blank')]
    yMax = 0
    plotCount = 0
    for row in range(3):
        for col in range(3):
            ax = fig.add_subplot(gs00[row,col])
            locDir = locDirList[plotCount]
            loc0Con = highContrast
            loc1Con = highContrast
            loc0Dir, loc1Dir = locDir
            loc0Title = loc0Dir
            loc1Title = loc1Dir
            if loc0Dir == 'blank':
                loc0Dir = zeroDir
                loc0Con = zeroContrast
            if loc1Dir == 'blank':
                loc1Dir = zeroDir
                loc1Con = zeroContrast
            histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) & 
                                          (stimIndexDF['loc0 Contrast'] == loc0Con) &
                                          (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                          (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
            dirPlot = spikeHists[unitCount, histIndex, :] * 1000/stimIndexCount[histIndex]
            smoothPlot = gaussian_filter1d(dirPlot, 5)
            if max(smoothPlot) > yMax:
                yMax = max(smoothPlot)
            ax.plot(smoothPlot)
            ax.set_title(f'loc0: {loc0Title}, loc1: {loc1Title}', fontsize=5)
            ax.set_xticks([0,
                           histPrePostMS,
                           histPrePostMS+trueStimDurMS,
                           2*histPrePostMS+trueStimDurMS])
            ax.set_xticklabels([])
            ax.set_xlim([0, trueStimDurMS+(2*histPrePostMS+1)])
            ax.axvspan(histPrePostMS, histPrePostMS+trueStimDurMS,
                       color='grey', alpha=0.1)
            ax.axhline(y=meanSpike[unitCount][48]*1000/trueStimDurMS,
                       linestyle='--', color='grey')
            if plotCount == 6:
                ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax.set_xlabel('Stimulus Duration (ms)', fontsize=7)
                ax.set_xticklabels([-histPrePostMS, 0, 0+trueStimDurMS,
                                    trueStimDurMS+histPrePostMS],
                                   fontsize=7)
            plotCount += 1

        axes = fig.get_axes()
        for ax in axes:
            ax.set_ylim([0, yMax*1.1])


    # Normalization Plot
    # ReIndexing to have pref direction in the middle Gauss Smooth
    nullIndex = np.where(dirArray == nullDir)[0][0]
    reIndex = (np.array([0, 1, 2, 3, 4, 5])+nullIndex) % 6
    tickLabels = np.array(['0', '60', '120', '180', '240', '300'])[reIndex]
    tickLabels = np.append(tickLabels, ['blank'])

    # smoothed data
    bSmoothReIndex = np.zeros((7, 7))
    tempMain = bSmooth[:6, :6][:, reIndex]
    tempMain = tempMain[:6, :6][reIndex, :]
    temp0Blank = bSmooth[:6, 6][reIndex]
    temp1Blank = bSmooth[6, :6][reIndex]
    bSmoothReIndex[:6, :6] = tempMain
    bSmoothReIndex[:6, 6] = temp0Blank
    bSmoothReIndex[6, :6] = temp1Blank
    bSmoothReIndex[6, 6] = bSmooth[6, 6]

    vMax = np.max(bSmoothReIndex)

    # raw data
    bReIndex = np.zeros((7, 7))
    tempMain = b[:6, :6][:, reIndex]
    tempMain = tempMain[:6, :6][reIndex, :]
    temp0Blank = b[:6, 6][reIndex]
    temp1Blank = b[6, :6][reIndex]
    bReIndex[:6, :6] = tempMain
    bReIndex[:6, 6] = temp0Blank
    bReIndex[6, :6] = temp1Blank
    bReIndex[6, 6] = b[6, 6]

    if np.max(bReIndex) > vMax:
        vMax = np.max(bReIndex)

    # fitted data
    stimMat = np.zeros((7, 7))
    stimMat[:6, :6] = np.arange(36).reshape(6, 6)
    stimMat[6, :6] = np.arange(36, 42)
    stimMat[:, 6] = np.arange(42, 49)

    stimMatReIndex = np.zeros((7, 7))
    tempMain = stimMat[:6, :6][:, reIndex]
    tempMain = np.squeeze(tempMain[:6, :6][reIndex, :])
    temp0Blank = stimMat[:6, 6][reIndex]
    temp1Blank = stimMat[6, :6][reIndex]
    stimMatReIndex[:6, :6] = tempMain
    stimMatReIndex[:6, 6] = temp0Blank
    stimMatReIndex[6, :6] = temp1Blank
    stimMatReIndex[6, 6] = stimMat[6, 6]

    y_pred = np.zeros(49)
    # fitting Gaussian EMS to paired stimuli
    fixedVals = fixedValsForGenericNorm(stimMatReIndex.reshape(49)[:-1],
                                        stimIndexDict)

    y_pred[:-1] = genericNormNoScalar(fixedVals, *unitNormFitEstimate[unitCount])
    y_pred[-1] = b[6, 6]
    y_pred = y_pred.reshape((7, 7))

    if np.max(y_pred) > vMax:
        vMax = np.max(y_pred)

    # residual (diff between real vs predicted resp)
    residual = abs(y_pred - bReIndex)
    if np.max(residual) > vMax:
        vMax = np.max(residual)

    # scatter plot of fit responses vs real responses
    respReal = bReIndex.reshape(49)
    respFit = y_pred.reshape(49)
    ax2 = fig.add_subplot(gs01[0,1])
    ax2.scatter(respFit, respReal)
    ax2.set_ylabel('Real Responses (spikes/sec)', fontsize=7)
    ax2.set_xlabel('Fit Responses (spikes/sec)', fontsize=7)
    ax2.xaxis.set_label_position('top')
    ax2.set_ylim([-2,np.max(respReal)*1.10])
    ax2.set_xlim([-2,np.max(respReal)*1.10])
    line = lines.Line2D([0, 1], [0, 1], color='red')
    transform = ax2.transAxes
    line.set_transform(transform)
    ax2.add_line(line)
    # ax2 = plt.axis('equal')

    ax3 = fig.add_subplot(gs01[0,0])
    ax3 = sns.heatmap(bReIndex, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
    ax3.set_xticks(np.arange(7)+0.5)
    ax3.set_title(f'Raw Data', y=-0.1, fontsize=7)
    ax3.set_xlabel('Location 0 Stimulus Direction', fontsize=7)
    ax3.xaxis.set_label_position('top')
    ax3.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax3.set_ylabel('Location 1 Stimulus Direction', fontsize=7)
    ax3.xaxis.set_ticks_position("top")
    ax3.set_yticks(np.arange(7)+0.5)
    ax3.set_yticklabels(tickLabels, rotation=0, fontsize=7)

    ax6 = fig.add_subplot(gs01[1, 0])
    ax6 = sns.heatmap(y_pred, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
    ax6.set_xticks(np.arange(7) + 0.5)
    ax6.set_title(f'Predicted', y=-0.1, fontsize=7)
    ax6.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax6.set_ylabel('Location 1 Stimulus Direction', fontsize=7)
    ax6.xaxis.set_ticks_position("top")
    ax6.set_yticks(np.arange(7) + 0.5)
    ax6.set_yticklabels(tickLabels, rotation=0, fontsize=7)

    ax7 = fig.add_subplot(gs01[1, 1])
    ax7 = sns.heatmap(residual, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
                      # cbar=True, cbar_kws={'orientation': 'horizontal'})
    ax7.set_xticks(np.arange(7) + 0.5)
    ax7.set_title(f'Residual (Predicted-Raw)', y=-0.1, fontsize=7)
    ax7.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax7.xaxis.set_ticks_position("top")
    ax7.set_yticks(np.arange(7) + 0.5)
    ax7.set_yticklabels(tickLabels, rotation=0, fontsize=7)


    # bar plot
    ax4 = fig.add_subplot(gs02[0,0])
    unitDF = pnContrastPlotDF.loc[pnContrastPlotDF['unit'] == unit]

    #spikeCounts in spikes/sec
    unitDF['stimSpikes'] = unitDF['stimSpikes'] * 1000/trueStimDurMS 

    offset = lambda p: transforms.ScaledTranslation(p/72., 0, plt.gcf().dpi_scale_trans)
    trans = plt.gca().transData
    offsetCount = 0

    colorList = ['b', 'g', 'r', 'c']
    # plot category by category
    yMax = 0
    for pos in [0.35, 0.65]:
        if pos == 0.35:
            for count, indx in enumerate(['pref+blank', 'null+blank',
                                          'blank+pref', 'blank+null']):
                mean = unitDF.loc[(unitDF['prefNullStr'] == indx)].mean()[3]
                sem = unitDF.loc[(unitDF['prefNullStr'] == indx)].sem()[3]
                ax4 = plt.scatter(pos, mean, transform=trans+offset(offsetCount),
                                  label=indx, color=colorList[count])
                ax4 = plt.errorbar(pos, mean, yerr=sem, fmt='o',
                                   transform=trans+offset(offsetCount),
                                   color=colorList[count])
                offsetCount += 5
                stimIndex = unitDF.loc[(unitDF['prefNullStr'] == indx), 'stimIndex'].iloc[0]

                # # using EMS generic for L0-L6
                # fixedVals = fixedValsForEMSGen(stimIndex, stimIndexDict)
                # if scaledLoc[unitCount] == 0:
                #     y_pred = emsNormGen1(fixedVals, *unitNormFitEstimate[unitCount])
                # else:
                #     y_pred = emsNormGen0(fixedVals, *unitNormFitEstimate[unitCount])

                # using gaussian EMS for single stim
                fixedVals = fixedValsForGenericNorm(stimIndex, stimIndexDict)
                y_pred = genericNormNoScalar(fixedVals, *unitNormFitEstimate[unitCount])
                ax4 = plt.scatter(pos, y_pred, transform=trans+offset(offsetCount),
                                  label=indx, color=lightenColor(colorList[count], 0.5))
                offsetCount += 5
                if mean > yMax:
                    yMax = mean
        else:
            offsetCount = 0
            for count, indx in enumerate(['pref+pref', 'pref+null',
                                          'null+pref', 'null+null']):
                mean = unitDF.loc[(unitDF['prefNullStr'] == indx)].mean()[3]
                sem = unitDF.loc[(unitDF['prefNullStr'] == indx)].sem()[3]
                ax4 = plt.scatter(pos,mean, transform=trans+offset(offsetCount),
                                  label=indx, color=colorList[count])
                ax4 = plt.errorbar(pos, mean, yerr=sem, fmt="o",
                                   transform=trans+offset(offsetCount),
                                   color=colorList[count])
                offsetCount += 5
                stimIndex = unitDF.loc[(unitDF['prefNullStr'] == indx), 'stimIndex'].iloc[0]
                fixedVals = fixedValsForGenericNorm(stimIndex, stimIndexDict)
                y_pred = genericNormNoScalar(fixedVals, *unitNormFitEstimate[unitCount])
                ax4 = plt.scatter(pos, y_pred, transform=trans+offset(offsetCount),
                                  label=indx, color=lightenColor(colorList[count], 0.5))
                offsetCount += 5
                if mean > yMax:
                    yMax = mean
    ax4 = plt.axhline(y=meanSpike[unitCount][48]*1000/trueStimDurMS, linestyle='--', color='grey')
    ax4 = plt.ylim([0,yMax*1.5])
    ax4 = plt.xticks([0.35,0.65], ['Singular Stimulus', 'Dual Stimulus'])
    ax4 = plt.xlim(left=0.2,right=0.8)
    ax4 = plt.ylabel('Firing Rate spikes/sec')
    ax4 = plt.legend(loc='upper right', prop={'size': 6}, bbox_to_anchor=(1.25, 1.0))

    # norm fit parameters with EMS generic (l1-l6)
    ax5 = fig.add_subplot(gs03[0,0])
    ax5.text(0.5,0.5, f'the scaled location is: {int(scaledLoc[unitCount])}\n\
    L0_0: {unitNormFitEstimate[unitCount][0]:.2f}\n\
    L0_60: {unitNormFitEstimate[unitCount][1]:.2f}\n\
    L0_120: {unitNormFitEstimate[unitCount][2]:.2f}\n\
    L0_180: {unitNormFitEstimate[unitCount][3]:.2f}\n\
    L0_240: {unitNormFitEstimate[unitCount][4]:.2f}\n\
    L0_300: {unitNormFitEstimate[unitCount][5]:.2f}\n\
    L1_0: {unitNormFitEstimate[unitCount][6]:.2f}\n\
    L1_60: {unitNormFitEstimate[unitCount][7]:.2f}\n\
    L1_120: {unitNormFitEstimate[unitCount][8]:.2f}\n\
    L1_180: {unitNormFitEstimate[unitCount][9]:.2f}\n\
    L1_240: {unitNormFitEstimate[unitCount][10]:.2f}\n\
    L1_300: {unitNormFitEstimate[unitCount][11]:.2f}\n\
    alpha loc0: {unitNormFitEstimate[unitCount][12]:.2f}\n\
    alpha loc1: {unitNormFitEstimate[unitCount][13]:.2f}\n\
    fit R2: {unitPairedNormR2[unitCount]:.2f}', size=10, ha='center', transform=ax5.transAxes)
    ax5.axis('off')

    ax8 = fig.add_subplot(gs03[1, 0], polar='True')
    theta = np.radians(np.arange(0, 420, 60))
    sponTheta = np.radians(np.arange(0, 360, 360/100))
    sponTheta = np.append(sponTheta, sponTheta[0])

    # loc 0
    r0 = (np.append(meanSpike[unitCount][36:42], meanSpike[unitCount][36])) \
          * 1000/trueStimDurMS
    er0 = (np.append(spikeCountSEM[unitCount][36:42], spikeCountSEM[unitCount][36])) \
          * 1000/trueStimDurMS
    ax8.plot(theta, r0, markersize=2, color='green', label='location 0')
    ax8.errorbar(theta, r0, yerr=er0, fmt='o', ecolor='green',
                 color='green', markersize=2)

    # loc 1
    r1 = (np.append(meanSpike[unitCount][42:48], meanSpike[unitCount][42])) \
        * 1000/trueStimDurMS
    er1 = (np.append(spikeCountSEM[unitCount][42:48], spikeCountSEM[unitCount][42])) \
        * 1000/trueStimDurMS
    ax8.plot(theta, r1, markersize=2, color='red', label='location 1')
    ax8.errorbar(theta, r1, yerr=er1, fmt='x', ecolor='red',
                 color='red', markersize=2)
    ax8.set_theta_zero_location('W')
    ax8.set_title('Direction tuning at both locations')
    spon = np.array([meanSpike[unitCount][48] * 1000/trueStimDurMS] * len(sponTheta))
    ax8.plot(sponTheta, spon, linestyle='--', color='blue',
             label='spontaneous rate')
    ax8.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0,
               prop={'size': 6})

    plt.tight_layout()

    plt.savefig(f'{unit}EMS2PartsuperPlot.pdf')
    plt.close('all')


## 7x7 plot of all hists

for unitCount, unit in enumerate(units):

    fig, ax = plt.subplots(ncols=7, nrows=7, constrained_layout=True)

    count = 0
    yMax = 0
    for i in range(6):
        for j in range(6):
            dirPlot = spikeHists[unitCount, count, :] * 1000/stimIndexCount[count]
            smoothPlot = gaussian_filter1d(dirPlot, 5)
            if max(smoothPlot) > yMax:
                yMax = max(smoothPlot)
            ax[i,j].plot(smoothPlot)
            ax[i,j].set_xticks([0, histPrePostMS, histPrePostMS+trueStimDurMS,
                                2*histPrePostMS+trueStimDurMS])
            ax[i,j].set_xticklabels([])
            ax[i,j].set_xlim([0,trueStimDurMS+(2*histPrePostMS+1)])
            ax[i,j].axvspan(histPrePostMS, histPrePostMS+trueStimDurMS,
                            color='grey', alpha=0.1)
            ax[i,j].axhline(y=meanSpike[unitCount][48]*1000/trueStimDurMS,
                            linestyle='--', color='grey')
            count += 1

    for subCount, x in enumerate(range(36,42)):
        dirPlot = spikeHists[unitCount, x, :] * 1000 / stimIndexCount[x]
        smoothPlot = gaussian_filter1d(dirPlot, 5)
        if max(smoothPlot) > yMax:
            yMax = max(smoothPlot)
        ax[6, subCount].plot(smoothPlot)
        ax[6, subCount].set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                                    2 * histPrePostMS + trueStimDurMS])
        ax[6, subCount].set_xticklabels([])
        ax[6, subCount].set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
        ax[6, subCount].axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                                color='grey', alpha=0.1)
        ax[6, subCount].axhline(y=meanSpike[unitCount][48] * 1000 / trueStimDurMS,
                                linestyle='--', color='grey')
        if subCount == 0:
            ax[6, subCount].set_xlabel('Time (ms)')
            ax[6, subCount].set_ylabel('Firing Rate (spikes/sec)')
            ax[6, subCount].set_xticklabels([-(histPrePostMS), 0, 0 + trueStimDurMS,
                                            trueStimDurMS + histPrePostMS], fontsize=7)

    for subCount, x in enumerate(range(42,49)):
        dirPlot = spikeHists[unitCount, x, :] * 1000 / stimIndexCount[x]
        smoothPlot = gaussian_filter1d(dirPlot, 5)
        if max(smoothPlot) > yMax:
            yMax = max(smoothPlot)
        ax[subCount, 6].plot(smoothPlot)
        ax[subCount, 6].set_xticks([0, histPrePostMS, histPrePostMS + trueStimDurMS,
                                    2 * histPrePostMS + trueStimDurMS])
        ax[subCount, 6].set_xticklabels([])
        ax[subCount, 6].set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
        ax[subCount, 6].axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                                color='grey', alpha=0.1)
        ax[subCount, 6].axhline(y=meanSpike[unitCount][48] * 1000 / trueStimDurMS,
                                linestyle='--', color='grey')

    axes = fig.get_axes()
    for ax in axes:
        ax.set_ylim([0,yMax*1.1])

    plt.show()


###################################### STOP ########################################
###################################### HERE ########################################
'''


EXTRA CODE 


'''

# correlations for each paired stimulus condition across pairs:
tempMat = stats.zscore(spikeCountMat[:,:blocksDone,:], axis=1, nan_policy='omit')
individualCorrMat = np.zeros((len(combs), 36))
for count, i in enumerate(combs):
    n1 = np.where(units == i[0])[0][0]
    n2 = np.where(units == i[1])[0][0]
    for x in range(36):
        pairStimCorr = stats.pearsonr(tempMat[n1, :, x], tempMat[n2, :, x])
        individualCorrMat[count][x] = pairStimCorr[0]
individualCorrMat = individualCorrMat.reshape(len(combs) * 36)

## Plot Polar plot of dir tuning for each unit in both locations
for unitCount, unit in enumerate(units):

    ## figure
    fig, ax1 = plt.subplots(subplot_kw={'projection': 'polar'})
    theta = np.radians(np.arange(0,420,60))
    sponTheta = np.radians(np.arange(0,360,360/100))
    sponTheta = np.append(sponTheta, sponTheta[0])

    #loc 0
    r0 = (np.append(meanSpike[unitCount][36:42], meanSpike[unitCount][36])) \
         * 1000/trueStimDurMS
    er0 = (np.append(spikeCountSEM[unitCount][36:42], spikeCountSEM[unitCount][36])) \
           * 1000/trueStimDurMS
    ax1.plot(theta,r0, markersize=2, color='green', label='location 0')
    ax1.errorbar(theta, r0, yerr=er0, fmt='o', ecolor='green',
                 color='green', markersize=2)

    # loc 1
    r1 = (np.append(meanSpike[unitCount][42:48],meanSpike[unitCount][42])) \
        * 1000/trueStimDurMS
    er1 = (np.append(spikeCountSEM[unitCount][42:48], spikeCountSEM[unitCount][42])) \
        * 1000/trueStimDurMS
    ax1.plot(theta,r1, markersize=2, color='red', label='location 1')
    ax1.errorbar(theta, r1, yerr=er1, fmt='x', ecolor='red',
                 color='red', markersize=2)
    ax1.set_theta_zero_location('W')
    ax1.set_title('Direction tuning at both locations')
    spon = np.array([meanSpike[unitCount][48]] * len(sponTheta))
    ax1.plot(sponTheta,spon*1000/sponWindowMS, linestyle='--', color='blue',
             label='spontaneous spikes/sec')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

    plt.show()


# compiledPD = pd.DataFrame(compiledArr, columns=['correlations','pairSimScore',
#                           'pairSimPrefDir','pairLocSimScore','pairCombinedSimScore',
#                           'pairCombinedSimScoreMulti','pairCombinedSimScorePrefDir',])


# Plot Noise Correlations as a function of combined location similarity
pairCombinedSimScore = np.squeeze((pairLocSimScore + pairSimScore) / 2)
plt.scatter(pairCombinedSimScore, corrMat)

pairCombinedSimScoreMulti = np.squeeze((pairLocSimScore * pairSimScore))
plt.scatter(pairCombinedSimScoreMulti, corrMat)

pairCombinedSimScorePrefDir = np.squeeze((pairLocSimScore + pairSimPrefDir) / 2)
plt.scatter(pairCombinedSimScorePrefDir, corrMat)

m, b = np.polyfit(pairCombinedSimScore, corrMat, 1)
plt.plot(pairCombinedSimScore, m*pairCombinedSimScore+b)
pearsonR, pValue = stats.pearsonr(pairCombinedSimScore,corrMat)

plt.axis('equal')
plt.xlabel('Averaged Sim Score (dir tuning and RF location')
plt.ylabel("Pearson's Correlation")
plt.show()

# unit Direction Selectivity
unitsBaselineMat = np.load('../Direction Tuning/unitsBaselineMat.npy')
unitSelectivity = np.zeros(len(units))
for unit in filterUnits:
    unitIndex = np.where(units == unit)[0][0]
    nMax = np.where(dirTuningMat[unitIndex] == np.max(dirTuningMat[unitIndex]))[0].squeeze()
    prefDir = dirArray[nMax]
    nullDir = (prefDir + 180) % 360
    nullResp = dirTuningMat[unitIndex][np.where(dirArray==nullDir)[0][0]]
    prefResp = dirTuningMat[unitIndex][nMax]
    baseResp = unitsBaselineMat[unitIndex]
    DS = 1 - ((nullResp-baseResp)/(prefResp-baseResp))
    unitSelectivity[unitIndex] = DS

# pair Selectivity
combs = [i for i in combinations(filterUnits, 2)]
pairSelScore = np.zeros((len(combs),1))
for pairCount, pair in enumerate(combs):
    n1 = np.where(units == pair[0])[0][0]
    n2 = np.where(units == pair[1])[0][0]
    geoMean = np.sqrt((unitSelectivity[n1]*unitSelectivity[n2]))
    pairSelScore[pairCount] = geoMean


# generate paired stimulus correlation, selectivity index,
# and suppression index for each gabor pair for each pair
# of neurons (similar to Bram fig 2c)
pairCorr = []
pairSelectivityIndex = []
pairNonPrefSuppIndex = []
tempMat = stats.zscore(spikeCountMat[:, :blocksDone, :], axis=1, nan_policy='omit')
for pairCount, pair in enumerate(combs):
    n1 = np.where(units == pair[0])[0][0]
    n2 = np.where(units == pair[1])[0][0]
    for i in range(36):
        # correlation for that pair b/w 2 units
        # pairStimCorr = stats.pearsonr(tempMat[n1, :, i], tempMat[n2, :, i])
        pairStimCorr = stats.pearsonr(spikeCountMat[n1, :blocksDone, i],
                                      spikeCountMat[n2, :blocksDone, i])

        # extract directions of the gabor pair
        loc0Dir = stimIndexDict[i][0]['direction']
        loc1Dir = stimIndexDict[i][1]['direction']

        # n1 selectivity and suppression index
        n1ScaledLoc = int(scaledLoc[n1])
        n1Scalar = unitNormFitEstimate[n1][6]
        loc0Resp = unitNormFitEstimate[n1][np.where(dirArray == loc0Dir)[0][0]]
        loc1Resp = unitNormFitEstimate[n1][np.where(dirArray == loc1Dir)[0][0]]
        if n1ScaledLoc == 0:
            loc0Resp = loc0Resp * n1Scalar
        else:
            loc1Resp = loc1Resp * n1Scalar
        if loc0Resp > loc1Resp:
            n1Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
            n1NonPrefSupp = unitNormFitEstimate[n1][8] / (
                unitNormFitEstimate[n1][8] + unitNormFitEstimate[n1][7])
            n1PrefResp = 0
        else:
            n1Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
            n1NonPrefSupp = unitNormFitEstimate[n1][7] / (
                unitNormFitEstimate[n1][8] + unitNormFitEstimate[n1][7])
            n1PrefResp = 1

        # n2 selectivity and suppression index
        n2ScaledLoc = int(scaledLoc[n2])
        n2Scalar = unitNormFitEstimate[n2][6]
        loc0Resp = unitNormFitEstimate[n2][np.where(dirArray == loc0Dir)[0][0]]
        loc1Resp = unitNormFitEstimate[n2][np.where(dirArray == loc1Dir)[0][0]]
        if n2ScaledLoc == 0:
            loc0Resp = loc0Resp * n2Scalar
        else:
            loc1Resp = loc1Resp * n2Scalar
        if loc0Resp > loc1Resp:
            n2Selectivity = (loc0Resp - loc1Resp) / (loc0Resp + loc1Resp)
            n2NonPrefSupp = unitNormFitEstimate[n2][8] / (
                unitNormFitEstimate[n2][8] + unitNormFitEstimate[n2][7])
            n2PrefResp = 0
        else:
            n2Selectivity = (loc1Resp - loc0Resp) / (loc0Resp + loc1Resp)
            n2NonPrefSupp = unitNormFitEstimate[n2][7] / (
                unitNormFitEstimate[n2][8] + unitNormFitEstimate[n2][7])
            n2PrefResp = 1

        # pair selectivity and suppression index
        pairSelectivity = np.sqrt(n1Selectivity * n2Selectivity)
        if n1PrefResp != n2PrefResp:
            pairSelectivity = -pairSelectivity
        pairSuppression = np.sqrt(n1NonPrefSupp * n2NonPrefSupp)

        pairCorr.append(pairStimCorr[0])
        pairSelectivityIndex.append(pairSelectivity)
        pairNonPrefSuppIndex.append(pairSuppression)


## superPlot Original with EMS gauss and two step fitting
## SUPERPLOT OF PSTH, Normalization (pref+null+blank) heatmap, and bar plot
## Create subset of Pandas DF for P, N, P+N conditions for all units
pnContrastPlotDF = pd.DataFrame(columns=['unit', 'stimIndex', 'stimCount',
                                         'stimSpikes', 'contrast', 'prefNullStr'])
for unitCount, unit in enumerate(units):
    # find direction tested that is closest to unit's preferred and null direction
    prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
    orientCount = 0
    orientList = ['pref+pref', 'pref+null', 'null+pref', 'null+null',
                  'pref+blank', 'null+blank', 'blank+pref', 'blank+null']
    for j in [(prefDir, prefDir), (prefDir, nullDir), (nullDir, prefDir), (nullDir, nullDir), ]:
        loc0Dir, loc1Dir = j
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == highContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == highContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                                  & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = highContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1

    # Loc 0 Pref/Null only
    for i in [(prefDir, zeroDir), (nullDir, zeroDir)]:
        loc0Dir, loc1Dir = i
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == highContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == zeroContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                                  & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = zeroContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1

    # loc 1 Pref/Null only
    for x in [(zeroDir, prefDir), (zeroDir, nullDir)]:
        loc0Dir, loc1Dir = x
        sIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                   (stimIndexDF['loc0 Contrast'] == zeroContrast) &
                                   (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                   (stimIndexDF['loc1 Contrast'] == highContrast)][0]
        unitDF = spikeCountDF.loc[(spikeCountDF['unit'] == unit)
                                  & (spikeCountDF['stimIndex'] == sIndex)]
        unitDF = unitDF.iloc[:blocksDone]
        unitDF['contrast'] = zeroContrast
        unitDF['prefNullStr'] = orientList[orientCount]
        pnContrastPlotDF = pd.concat([pnContrastPlotDF, unitDF])
        orientCount += 1

for unitCount, unit in enumerate(units):
    b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000 / trueStimDurMS
    bSmooth = gaussian_filter(b, sigma=1)
    prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])

    ## figure
    # fig = plt.figure(constrained_layout=True)
    fig = plt.figure()
    fig.set_size_inches(17, 7)

    gs0 = gridspec.GridSpec(1, 5)

    heights = [1, 1, 1]
    widths = [4, 4, 4]
    gs00 = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs0[0],
                                            width_ratios=widths,
                                            height_ratios=heights)
    gs01 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[1:3])
    gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[3])
    gs03 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[4])

    # PSTHs Plot
    locDirList = [(nullDir, nullDir), (prefDir, nullDir), ('blank', nullDir),
                  (nullDir, prefDir), (prefDir, prefDir), ('blank', prefDir),
                  (nullDir, 'blank'), (prefDir, 'blank'), ('blank', 'blank')]
    yMax = 0
    plotCount = 0
    for row in range(3):
        for col in range(3):
            ax = fig.add_subplot(gs00[row, col])
            locDir = locDirList[plotCount]
            loc0Con = highContrast
            loc1Con = highContrast
            loc0Dir, loc1Dir = locDir
            loc0Title = loc0Dir
            loc1Title = loc1Dir
            if loc0Dir == 'blank':
                loc0Dir = zeroDir
                loc0Con = zeroContrast
            if loc1Dir == 'blank':
                loc1Dir = zeroDir
                loc1Con = zeroContrast
            histIndex = stimIndexDF.index[(stimIndexDF['loc0 Direction'] == loc0Dir) &
                                          (stimIndexDF['loc0 Contrast'] == loc0Con) &
                                          (stimIndexDF['loc1 Direction'] == loc1Dir) &
                                          (stimIndexDF['loc1 Contrast'] == loc1Con)][0]
            dirPlot = spikeHists[unitCount, histIndex, :] * 1000 / stimIndexCount[histIndex]
            smoothPlot = gaussian_filter1d(dirPlot, 5)
            if max(smoothPlot) > yMax:
                yMax = max(smoothPlot)
            ax.plot(smoothPlot)
            ax.set_title(f'loc0: {loc0Title}, loc1: {loc1Title}', fontsize=5)
            ax.set_xticks([0,
                           histPrePostMS,
                           histPrePostMS + trueStimDurMS,
                           2 * histPrePostMS + trueStimDurMS])
            ax.set_xticklabels([])
            ax.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
            ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                       color='grey', alpha=0.1)
            ax.axhline(y=meanSpike[unitCount][48] * 1000 / trueStimDurMS,
                       linestyle='--', color='grey')
            if plotCount == 6:
                ax.set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax.set_xlabel('Stimulus Duration (ms)', fontsize=7)
                ax.set_xticklabels([-histPrePostMS, 0, 0 + trueStimDurMS,
                                    trueStimDurMS + histPrePostMS],
                                   fontsize=7)
            plotCount += 1

        axes = fig.get_axes()
        for ax in axes:
            ax.set_ylim([0, yMax * 1.1])

    # Normalization Plot
    # ReIndexing to have pref direction in the middle Gauss Smooth
    nullIndex = np.where(dirArray == nullDir)[0][0]
    reIndex = (np.array([0, 1, 2, 3, 4, 5]) + nullIndex) % 6
    tickLabels = np.array(['0', '60', '120', '180', '240', '300'])[reIndex]
    tickLabels = np.append(tickLabels, ['blank'])

    # smoothed data
    bSmoothReIndex = np.zeros((7, 7))
    tempMain = bSmooth[:6, :6][:, reIndex]
    tempMain = tempMain[:6, :6][reIndex, :]
    temp0Blank = bSmooth[:6, 6][reIndex]
    temp1Blank = bSmooth[6, :6][reIndex]
    bSmoothReIndex[:6, :6] = tempMain
    bSmoothReIndex[:6, 6] = temp0Blank
    bSmoothReIndex[6, :6] = temp1Blank
    bSmoothReIndex[6, 6] = bSmooth[6, 6]

    vMax = np.max(bSmoothReIndex)

    # raw data
    bReIndex = np.zeros((7, 7))
    tempMain = b[:6, :6][:, reIndex]
    tempMain = tempMain[:6, :6][reIndex, :]
    temp0Blank = b[:6, 6][reIndex]
    temp1Blank = b[6, :6][reIndex]
    bReIndex[:6, :6] = tempMain
    bReIndex[:6, 6] = temp0Blank
    bReIndex[6, :6] = temp1Blank
    bReIndex[6, 6] = b[6, 6]

    if np.max(bReIndex) > vMax:
        vMax = np.max(bReIndex)

    # fitted data
    stimMat = np.zeros((7, 7))
    stimMat[:6, :6] = np.arange(36).reshape(6, 6)
    stimMat[6, :6] = np.arange(36, 42)
    stimMat[:, 6] = np.arange(42, 49)

    stimMatReIndex = np.zeros((7, 7))
    tempMain = stimMat[:6, :6][:, reIndex]
    tempMain = np.squeeze(tempMain[:6, :6][reIndex, :])
    temp0Blank = stimMat[:6, 6][reIndex]
    temp1Blank = stimMat[6, :6][reIndex]
    stimMatReIndex[:6, :6] = tempMain
    stimMatReIndex[:6, 6] = temp0Blank
    stimMatReIndex[6, :6] = temp1Blank
    stimMatReIndex[6, 6] = stimMat[6, 6]

    y_pred = np.zeros((7, 7))
    # fitting Gaussian EMS to paired stimuli
    pairedStim = stimMatReIndex[:6, :6].reshape(36)
    fixedVals = fixedValsCurveFitForPairedStim(prefDir,
                                               pairedStim,
                                               stimIndexDict,
                                               unitGaussFit[unitCount])
    if scaledLoc[unitCount] == 0:
        pred = EMSGaussNormStepFunc1(fixedVals, *unitPairedNormFit[unitCount]
                                     ).reshape((6, 6))
    else:
        pred = EMSGaussNormStepFunc0(fixedVals, *unitPairedNormFit[unitCount]
                                     ).reshape((6, 6))
    y_pred[:6, :6] = pred

    # fitting Gaussian to single stimuli
    singleStim = np.concatenate((stimMatReIndex[6, :6], stimMatReIndex[:6, 6]), axis=0)
    fixedVals = fixedValsForCurveFit(prefDir, singleStim, stimIndexDict)

    if scaledLoc[unitCount] == 1:
        pred = gaussNormFunc0(fixedVals, *unitGaussFit[unitCount])
    else:
        pred = gaussNormFunc1(fixedVals, *unitGaussFit[unitCount])
    y_pred[6, :6] = pred[:6]
    y_pred[:6, 6] = pred[6:]

    # # using EMS generic, L0-L6
    # fixedVals = fixedValsForEMSGen(stimMatReIndex.reshape(49), stimIndexDict)
    # if scaledLoc[unitCount] == 0:
    #     y_pred = emsNormGen1(fixedVals, *unitNormFitEstimate[unitCount]).reshape((7,7))
    # else:
    #     y_pred = emsNormGen0(fixedVals, *unitNormFitEstimate[unitCount]).reshape((7,7))

    if np.max(y_pred) > vMax:
        vMax = np.max(y_pred)

    # residual (diff between real vs predicted resp)
    residual = abs(y_pred - bReIndex)
    if np.max(residual) > vMax:
        vMax = np.max(residual)

    # ax2 = fig.add_subplot(gs01[0,1])
    # ax2 = sns.heatmap(bSmoothReIndex, square=True, linewidths=0.2, vmin=0,
    #                   vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
    # ax2.set_xticks(np.arange(7)+0.5)
    # ax2.set_title(f'Gaussian Smoothed', y=-0.1, fontsize=7)
    # ax2.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    # ax2.xaxis.set_ticks_position('top')
    # ax2.set_xlabel('Location 0 Stimulus Direction', fontsize=7)
    # ax2.xaxis.set_label_position('top')
    # ax2.set_yticks(np.arange(7)+0.5)
    # ax2.set_yticklabels(tickLabels, rotation=0, fontsize=7)

    respReal = bReIndex.reshape(49)
    respFit = y_pred.reshape(49)
    ax2 = fig.add_subplot(gs01[0, 1])
    ax2.scatter(respFit, respReal)
    ax2.set_ylabel('Real Responses (spikes/sec)', fontsize=7)
    ax2.set_xlabel('Fit Responses (spikes/sec)', fontsize=7)
    ax2.xaxis.set_label_position('top')
    ax2.set_ylim([-2, np.max(respReal) * 1.10])
    ax2.set_xlim([-2, np.max(respReal) * 1.10])
    line = lines.Line2D([0, 1], [0, 1], color='red')
    transform = ax2.transAxes
    line.set_transform(transform)
    ax2.add_line(line)
    # ax2 = plt.axis('equal')

    ax3 = fig.add_subplot(gs01[0, 0])
    ax3 = sns.heatmap(bReIndex, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
    ax3.set_xticks(np.arange(7) + 0.5)
    ax3.set_title(f'Raw Data', y=-0.1, fontsize=7)
    ax3.set_xlabel('Location 0 Stimulus Direction', fontsize=7)
    ax3.xaxis.set_label_position('top')
    ax3.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax3.set_ylabel('Location 1 Stimulus Direction', fontsize=7)
    ax3.xaxis.set_ticks_position("top")
    ax3.set_yticks(np.arange(7) + 0.5)
    ax3.set_yticklabels(tickLabels, rotation=0, fontsize=7)

    ax6 = fig.add_subplot(gs01[1, 0])
    ax6 = sns.heatmap(y_pred, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
    ax6.set_xticks(np.arange(7) + 0.5)
    ax6.set_title(f'Predicted', y=-0.1, fontsize=7)
    ax6.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax6.set_ylabel('Location 1 Stimulus Direction', fontsize=7)
    ax6.xaxis.set_ticks_position("top")
    ax6.set_yticks(np.arange(7) + 0.5)
    ax6.set_yticklabels(tickLabels, rotation=0, fontsize=7)

    ax7 = fig.add_subplot(gs01[1, 1])
    ax7 = sns.heatmap(residual, square=True, linewidths=0.2, vmin=0,
                      vmax=vMax, annot=True, annot_kws={'fontsize': 7}, cbar=False)
    # cbar=True, cbar_kws={'orientation': 'horizontal'})
    ax7.set_xticks(np.arange(7) + 0.5)
    ax7.set_title(f'Residual (Predicted-Raw)', y=-0.1, fontsize=7)
    ax7.set_xticklabels(tickLabels, rotation=45, fontsize=7)
    ax7.xaxis.set_ticks_position("top")
    ax7.set_yticks(np.arange(7) + 0.5)
    ax7.set_yticklabels(tickLabels, rotation=0, fontsize=7)

    # bar plot
    ax4 = fig.add_subplot(gs02[0, 0])
    unitDF = pnContrastPlotDF.loc[pnContrastPlotDF['unit'] == unit]

    # spikeCounts in spikes/sec
    unitDF['stimSpikes'] = unitDF['stimSpikes'] * 1000 / trueStimDurMS

    offset = lambda p: transforms.ScaledTranslation(p / 72., 0, plt.gcf().dpi_scale_trans)
    trans = plt.gca().transData
    offsetCount = 0

    colorList = ['b', 'g', 'r', 'c']
    # plot category by category
    yMax = 0
    for pos in [0.35, 0.65]:
        if pos == 0.35:
            for count, indx in enumerate(['pref+blank', 'null+blank', 'blank+pref', 'blank+null']):
                mean = unitDF.loc[(unitDF['prefNullStr'] == indx)].mean()[3]
                sem = unitDF.loc[(unitDF['prefNullStr'] == indx)].sem()[3]
                ax4 = plt.scatter(pos, mean, transform=trans + offset(offsetCount),
                                  label=indx, color=colorList[count])
                ax4 = plt.errorbar(pos, mean, yerr=sem, fmt='o',
                                   transform=trans + offset(offsetCount),
                                   color=colorList[count])
                offsetCount += 5
                stimIndex = unitDF.loc[(unitDF['prefNullStr'] == indx), 'stimIndex'].iloc[0]

                # # using EMS generic for L0-L6
                # fixedVals = fixedValsForEMSGen(stimIndex, stimIndexDict)
                # if scaledLoc[unitCount] == 0:
                #     y_pred = emsNormGen1(fixedVals, *unitNormFitEstimate[unitCount])
                # else:
                #     y_pred = emsNormGen0(fixedVals, *unitNormFitEstimate[unitCount])

                # using gaussian EMS for single stim
                fixedVals = fixedValsForCurveFit(prefDir, stimIndex, stimIndexDict)
                if scaledLoc[unitCount] == 1:
                    y_pred = gaussNormFunc0(fixedVals, *unitGaussFit[unitCount])
                else:
                    y_pred = gaussNormFunc1(fixedVals, *unitGaussFit[unitCount])
                ax4 = plt.scatter(pos, y_pred, transform=trans + offset(offsetCount),
                                  label=indx, color=lightenColor(colorList[count], 0.5))
                offsetCount += 5
                if mean > yMax:
                    yMax = mean
        else:
            offsetCount = 0
            for count, indx in enumerate(['pref+pref', 'pref+null', 'null+pref', 'null+null']):
                mean = unitDF.loc[(unitDF['prefNullStr'] == indx)].mean()[3]
                sem = unitDF.loc[(unitDF['prefNullStr'] == indx)].sem()[3]
                ax4 = plt.scatter(pos, mean, transform=trans + offset(offsetCount),
                                  label=indx, color=colorList[count])
                ax4 = plt.errorbar(pos, mean, yerr=sem, fmt="o",
                                   transform=trans + offset(offsetCount),
                                   color=colorList[count])
                offsetCount += 5
                stimIndex = unitDF.loc[(unitDF['prefNullStr'] == indx), 'stimIndex'].iloc[0]
                # fixedVals = fixedValsForEMSGen(stimIndex, stimIndexDict)
                # if scaledLoc[unitCount] == 0:
                #     y_pred = emsNormGen1(fixedVals, *unitNormFitEstimate[unitCount])
                # else:
                #     y_pred = emsNormGen0(fixedVals, *unitNormFitEstimate[unitCount])
                fixedVals = fixedValsCurveFitForPairedStim(prefDir,
                                                           stimIndex,
                                                           stimIndexDict,
                                                           unitGaussFit[unitCount])
                if scaledLoc[unitCount] == 0:
                    y_pred = EMSGaussNormStepFunc1(fixedVals,
                                                   *unitPairedNormFit[unitCount])
                else:
                    y_pred = EMSGaussNormStepFunc0(fixedVals,
                                                   *unitPairedNormFit[unitCount])
                ax4 = plt.scatter(pos, y_pred, transform=trans + offset(offsetCount),
                                  label=indx, color=lightenColor(colorList[count], 0.5))
                offsetCount += 5
                if mean > yMax:
                    yMax = mean
    ax4 = plt.axhline(y=meanSpike[unitCount][48] * 1000 / trueStimDurMS, linestyle='--', color='grey')
    ax4 = plt.ylim([0, yMax * 1.5])
    ax4 = plt.xticks([0.35, 0.65], ['Singular Stimulus', 'Dual Stimulus'])
    ax4 = plt.xlim(left=0.2, right=0.8)
    ax4 = plt.ylabel('Firing Rate spikes/sec')
    ax4 = plt.legend(loc='upper right', prop={'size': 6}, bbox_to_anchor=(1.25, 1.0))

    # Normalization Fit Parameters for EMS with Gaussian
    ax5 = fig.add_subplot(gs03[0, 0])
    ax5.text(0.5, 0.5, f'the scaled location is: {int(scaledLoc[unitCount])}\n\
    gaussian offset: {unitGaussFit[unitCount][0]:.2f}\n\
    gaussian amplitude: {unitGaussFit[unitCount][1]:.2f}\n\
    gaussian mean: {unitGaussFit[unitCount][2]:.2f}\n\
    gaussian sigma: {unitGaussFit[unitCount][3]:.2f}\n\
    location scalar: {unitGaussFit[unitCount][4]:.2f}\n\
    normalization alpha 0: {unitPairedNormFit[unitCount][0]:.2f}\n\
    normalization alpha 1: {unitPairedNormFit[unitCount][1]:.2f}\n\
    gauss fit R2: {unitGaussR2[unitCount]:.2f}\n\
    EMS norm fit R2: {unitPairedNormR2[unitCount]:.2f}', size=10,
             ha='center', transform=ax5.transAxes)
    ax5.axis('off')
    # baseline: {unitNormFitEstimate[unitCount][8]:.2f}\n\

    # norm fit parameters with EMS generic (l1-l6)
    # ax5 = fig.add_subplot(gs03[0,0])
    # ax5.text(0.5,0.5, f'the scaled location is: {int(scaledLoc[unitCount])}\n\
    # L0: {unitNormFitEstimate[unitCount][0]:.2f}\n\
    # L60: {unitNormFitEstimate[unitCount][1]:.2f}\n\
    # L120: {unitNormFitEstimate[unitCount][2]:.2f}\n\
    # L180: {unitNormFitEstimate[unitCount][3]:.2f}\n\
    # L240: {unitNormFitEstimate[unitCount][4]:.2f}\n\
    # L300: {unitNormFitEstimate[unitCount][5]:.2f}\n\
    # Scalar: {unitNormFitEstimate[unitCount][6]:.2f}\n\
    # alpha 1: {unitNormFitEstimate[unitCount][7]:.2f}\n\
    # alpha 2: {unitNormFitEstimate[unitCount][8]:.2f}\n\
    # sigma: {unitNormFitEstimate[unitCount][9]:.2f}\n\
    # baseline: {unitNormFitEstimate[unitCount][10]:.2f}\n\
    # fit R2: {unitR2[unitCount]:.2f}', size=10, ha='center', transform=ax5.transAxes)
    # ax5.axis('off')

    ax8 = fig.add_subplot(gs03[1, 0], polar='True')
    theta = np.radians(np.arange(0, 420, 60))
    sponTheta = np.radians(np.arange(0, 360, 360 / 100))
    sponTheta = np.append(sponTheta, sponTheta[0])
    # loc 0
    r0 = (np.append(meanSpike[unitCount][36:42], meanSpike[unitCount][36])) \
         * 1000 / trueStimDurMS
    er0 = (np.append(spikeCountSEM[unitCount][36:42], spikeCountSEM[unitCount][36])) \
          * 1000 / trueStimDurMS
    ax8.plot(theta, r0, markersize=2, color='green', label='location 0')
    ax8.errorbar(theta, r0, yerr=er0, fmt='o', ecolor='green',
                 color='green', markersize=2)

    # loc 1
    r1 = (np.append(meanSpike[unitCount][42:48], meanSpike[unitCount][42])) \
         * 1000 / trueStimDurMS
    er1 = (np.append(spikeCountSEM[unitCount][42:48], spikeCountSEM[unitCount][42])) \
          * 1000 / trueStimDurMS
    ax8.plot(theta, r1, markersize=2, color='red', label='location 1')
    ax8.errorbar(theta, r1, yerr=er1, fmt='x', ecolor='red',
                 color='red', markersize=2)
    ax8.set_theta_zero_location('W')
    ax8.set_title('Direction tuning at both locations')
    spon = np.array([meanSpike[unitCount][48] * 1000 / trueStimDurMS] * len(sponTheta))
    ax8.plot(sponTheta, spon, linestyle='--', color='blue',
             label='spontaneous rate')
    ax8.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0,
               prop={'size': 6})

    plt.tight_layout()

    plt.savefig(f'{unit}EMS2PartsuperPlot.pdf')
    plt.close('all')

