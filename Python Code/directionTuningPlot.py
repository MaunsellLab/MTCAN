import scipy.io as sp
import numpy as np

allTrials = sp.loadmat('testa00.mat', squeeze_me = True)

mapFile = allTrials['file']
trials = allTrials['trials']

numDir = mapFile['mapSettings'].item()['data'].item()['directionDeg'].item()['n'].tolist()
stimDurMS = mapFile['mapStimDurationMS'].item()['data'].tolist()
maxNumStim = 1000
hisPrePostMS = 50
numStim = zeros(1,numDir)
spikeCounts = zeros(numDir, maxNumStim)
spikeHists = zeros(numDir, stimDurMS + 2 * histPrePostMS)

for trial in trials:
    


