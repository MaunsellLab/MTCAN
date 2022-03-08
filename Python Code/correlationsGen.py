# Pairwise Correlations
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

corr = []
trialsSpikes = np.zeros((20,2))
trialsSpikes[:,:] = np.nan
numNeurons = 2
fRate = 20 #spikes/sec
stimDur = 493

for j in range(100):
    for i in range(20):
        spikes = np.random.poisson(fRate/(1000/stimDur), numNeurons)
        trialsSpikes[i] = spikes

    R = np.zeros((numNeurons,numNeurons))
    for neuronI in range(numNeurons):
        for neuronJ in range(numNeurons):
            if neuronI == neuronJ:
                R[neuronI,neuronJ] = 1
            else:
                R[neuronI, neuronJ] = 0.1

    L = np.linalg.cholesky(R)
    trialsSpikes = np.matmul(trialsSpikes, L)

    stimCorr = np.corrcoef(trialsSpikes[:,0], trialsSpikes[:,1])
    corr.append(stimCorr[0][1])
print(np.mean(corr))



# neurons with diff tuning curves (diff firing rates)
corr = []
trialsSpikes = np.zeros((10,2))
trialsSpikes[:,:] = np.nan
numNeurons = 2
fRate1 = 25
fRate2 = 10
stimDur = 493
for j in range(100):

    for i in range(10):
        spike1 = np.random.poisson(fRate1/(1000/stimDur))
        spike2 = np.random.poisson(fRate2/(1000/stimDur))
        spikes = [spike1, spike2]
        spikes = np.array(spikes)

        R = np.zeros((numNeurons,numNeurons))
        for neuronI in range(numNeurons):
            for neuronJ in range(numNeurons):
                if neuronI == neuronJ:
                    R[neuronI,neuronJ] = 1
                else:
                    R[neuronI, neuronJ] = 0.1

        L = np.linalg.cholesky(R)
        trialsSpikes[i] = np.matmul(spikes,L)

    stimCorr = np.corrcoef(trialsSpikes[:,0], trialsSpikes[:,1])
    corr.append(stimCorr[0][1])

# alt strat for neurons with diff tuning curves (diff firing rates)
corr = []
trialsSpikes = np.zeros((10,2))
trialsSpikes[:,:] = np.nan
numNeurons = 2
fRate1 = 35
fRate2 = 25
sharedRate = 20
stimDur = 493
for j in range(100):
    for i in range(trialSpikes.shape[0]):
        if np.random.uniform(0,1) < 0.10:
            shareSpike = np.random.poisson(sharedRate/(1000/stimDur))
            spikes = [shareSpike, shareSpike]
            spikes = np.array(spikes)
            trialsSpikes[i] = spikes
        else:
            spike1 = np.random.poisson(fRate1/(1000/stimDur))
            spike2 = np.random.poisson(fRate2/(1000/stimDur))
            spikes = [spike1, spike2]
            
            trialsSpikes[i] = spikes
    
    stimCorr = np.corrcoef(trialsSpikes[:,0], trialsSpikes[:,1])
    corr.append(stimCorr[0][1])

print(np.mean(corr))
plt.hist(corr)
plt.show()




def insertStimSpikeData(units, index, stimOnTimeSNEV):
    '''
    this function will generate a poisson spike train and return the normalized
    response for the RF for a given stimulus configuration
    Inputs:
        x: neuron number
        index: the index of the corresponding stimulus configuration
        sigma: semisaturation constant from neuron's contrast response function
    Outputs:
        RFSpikes: normalized response 

    spikes_4rows = np.tile(spikes, (4,1))
    '''


    numNeurons = len(units)

    C0 = stimIndexDict[index][0]['contrast']
    C1 = stimIndexDict[index][1]['contrast']
    L0 = tcDict[x+1][stimIndexDict[index][0]['direction']]
    L1 = tcDict[x+1][stimIndexDict[index][1]['direction']]
    sigma = 0.1
    expectedNormSpikeRate = int(((C0*L0) + (C1*L1))/(C0 + C1 + sigma))
    stimDur = 493
    popMean = expectedNormSpikeRate/(1000/stimDur)

    spikes = popMean + np.random.rand(1,numNeurons)           
    R = np.zeros((numNeurons,numNeurons))
    for neuronI in range(numNeurons):
        for neuronJ in range(numNeurons):
            if neuronI == neuronJ:
                R[neuronI,neuronJ] = 1
            else:
                R[neuronI, neuronJ] = 0.1

    L = np.linalg.cholesky(R)
    spikes = np.matmul(spikes,L)   
    spikes = np.around(spikes) 


    for count, i in enumerate(spikes[0]):
        if count == 0 and i != 0
            unit = units[0]
            channelIdentity = int(unit[0:unit.find('_')])
            channel = np.array([channelIdentity] * stimDur)
            spikeTimeMS = (np.sort(np.random.choice(np.arange(stimDur), int(i),
            replace = False)))/1000
            currTrial['spikeData']['timeStamp'] = np.append(currTrial['spikeData'] \
            ['timeStamp'], stimOnTimeSNEV + spikeTimeMS, 0)
            currTrial['spikeData']['unit'] = np.append(currTrial['spikeData'] \
            ['unit'], [unit] * len(spikeTimeMS), 0)
            currTrial['spikeData']['channel'] = np.append(currTrial['spikeData'] \
            ['channel'], [channelIdentity] * len(spikeTimeMS), 0)
        elif count == 1 and i!= 0:
            unit = units[1]
            unit = units[0]
            channelIdentity = int(unit[0:unit.find('_')])
            channel = np.array([channelIdentity] * stimDur)
            spikeTimeMS = (np.sort(np.random.choice(np.arange(stimDur), int(i), 
            replace = False)))/1000
            currTrial['spikeData']['timeStamp'] = np.append(currTrial['spikeData'] \
            ['timeStamp'], stimOnTimeSNEV + spikeTimeMS, 0)
            currTrial['spikeData']['unit'] = np.append(currTrial['spikeData'] \
            ['unit'], [unit] * len(spikeTimeMS), 0)
            currTrial['spikeData']['channel'] = np.append(currTrial['spikeData'] \
            ['channel'], [channelIdentity] * len(spikeTimeMS), 0)
        elif count == 2 and i != 0:
            unit = units[2] 
            unit = units[0]
            channelIdentity = int(unit[0:unit.find('_')])
            channel = np.array([channelIdentity] * stimDur)
            spikeTimeMS = (np.sort(np.random.choice(np.arange(stimDur), int(i), 
            replace = False)))/1000
            currTrial['spikeData']['timeStamp'] = np.append(currTrial['spikeData'] \
            ['timeStamp'], stimOnTimeSNEV + spikeTimeMS, 0)
            currTrial['spikeData']['unit'] = np.append(currTrial['spikeData'] \
            ['unit'], [unit] * len(spikeTimeMS), 0)
            currTrial['spikeData']['channel'] = np.append(currTrial['spikeData'] \
            ['channel'], [channelIdentity] * len(spikeTimeMS), 0)

