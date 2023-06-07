# Imports
from usefulFns import *

# Load relevant file here with pyMat reader
monkeyName = 'Meetz'
seshDate = '230606'
fileName = f'{monkeyName}_{seshDate}_MTND_Spikes.mat'
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

# get mapping stimulus positions and index
stepIndex = []
loc1AziEle = []
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    stimDesc = currTrial['stimDesc']['data']
    for stim in stimDesc:
        if stim['stimLoc'] == 1 and stim['stepIndex'] not in stepIndex:
            stepIndex.append(stim['stepIndex'])
            loc1AziEle.append((stim['azimuthDeg'], stim['elevationDeg']))

stepIndex = np.array(stepIndex)

# initialize lists/arrays/dataframes for counting spikeCounts and for analysis
blocksDone = allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone'][0]
numOffsets = header['blockStatus']['data']['numOffsets'][0]
numSteps = int((numOffsets - 1) / 2)
spikeCountMat = np.zeros((len(units), blocksDone+1, numSteps*7+3))
stimCount = np.zeros((3, numOffsets), dtype=int)
stimCountIndex = np.arange(numSteps*7-2)
stimCountIndex = stimCountIndex.reshape(3, int(len(stimCountIndex)/3))
transectStimCount = np.zeros(numSteps, dtype=int)
transectCountIndex = np.arange(numSteps*7-2, numSteps*7+3)
onLatency = 50 / 1000  # time in MS for counting window latency after stim on
offLatency = 150 / 1000  # time in MS for counting window latency after stim off
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
meanSpikeReshaped = np.zeros((len(units), 3, int((numSteps*7-2)/3)))
for count, i in enumerate(meanSpikeReshaped):
    meanSpikeReshaped[count] = (meanSpike[count, :numSteps*7-2].reshape(3, int((numSteps*7-2)/3)) *
                                1000/trueStimDurMS)

# plot line graph of responses across distances
for count in range(len(units)):
    x = np.arange(0, numSteps+1)
    xNorm = np.arange(1, numSteps+1)
    nullOnly = np.concatenate(([meanSpikeReshaped[count, 1, 0]],
                              meanSpikeReshaped[count, 0, 1::2]), axis=0)
    prefOnly = np.concatenate(([meanSpikeReshaped[count, 2, 0]],
                              meanSpikeReshaped[count, 0, 2::2]), axis=0)
    NP = meanSpikeReshaped[count, 1, 1::2]
    PN = meanSpikeReshaped[count, 2, 2::2]
    plt.plot(x, prefOnly, color='black', label='Preferred Only')
    plt.plot(x, nullOnly, color='grey', label='Null Only')
    plt.plot(xNorm, NP, color='red', label='Null Center, Pref Offset')
    plt.plot(xNorm, PN, color='green', label='Pref Center, Null Offset')
    plt.xlabel('stimulus offset positions')
    plt.ylabel('Response (spikes/sec)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.show()

for count in range(len(units)):




# basic version
# initialize lists/arrays/dataframes for counting spikeCounts and for analysis
blocksDone = allTrials[corrTrials[-2]]['blockStatus']['data']['blocksDone']
numOffsets = header['blockStatus']['data']['numOffsets']
numSteps = int((numOffsets - 1) / 2)
stimCount = np.zeros((3, numOffsets))
transectStimCount = np.zeros(numSteps)
for corrTrial in corrTrials:
    currTrial = allTrials[corrTrial]
    stimDesc = currTrial['stimDesc']['data']
    for stim in stimDesc:
        centIndex = stim['centerIndex']
        offIdx = stim['offsetIndex']
        step = stim['stepIndex']
        if stim['listType'] == 1:
            if step > numSteps:
                if stim['stimLoc'] == 1:
                    indx = step - numSteps - 1
                    transectStimCount[indx] += 1
            else:
                if stim['stimLoc'] == 0:
                    col = (int(offIdx != 0) * step * 2
                           - int(offIdx == 1))
                    # ^ basically indexing into stimCount stepIndex * 2 - 1
                    # for every column after the first, which is the blank
                    # condition
                    stimCount[centIndex, col] += 1

                    if stim['centerIndex'] == 0 and stim['offsetIndex'] == 0:
                        stimCount[0, 0] += 1
                    if stim['centerIndex'] == 1 and stim['offsetIndex'] == 0:
                        stimCount[1, 0] += 1
                    if stim['centerIndex'] == 2 and stim['offsetIndex'] == 0:
                        stimCount[2, 0] += 1
                    if stim['centerIndex'] == 0 and stim['offsetIndex'] == 1:
                        stimCount[0, stim['stepIndex'] * 2 - 1] += 1
                    if stim['centerIndex'] == 0 and stim['offsetIndex'] == 2:
                        stimCount[0, stim['stepIndex'] * 2] += 1
                    if stim['centerIndex'] == 1 and stim['offsetIndex'] == 1:
                        stimCount[1, stim['stepIndex'] * 2 - 1] += 1
                    if stim['centerIndex'] == 1 and stim['offsetIndex'] == 2:
                        stimCount[1, stim['stepIndex'] * 2] += 1
                    if stim['centerIndex'] == 2 and stim['offsetIndex'] == 1:
                        stimCount[2, stim['stepIndex'] * 2 - 1] += 1
                    if stim['centerIndex'] == 2 and stim['offsetIndex'] == 2:
                        stimCount[2, stim['stepIndex'] * 2] += 1

