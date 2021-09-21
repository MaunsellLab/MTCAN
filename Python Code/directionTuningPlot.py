import scipy.io as sp
import numpy as np

allTrials = sp.loadmat('testa00.mat', squeeze_me = True)

file = allTrials['file']
mapSettings = file['mapSettings'].item()['data'].item()
directionDeg = mapSettings['directionDeg'].item()

numDir = file['mapSettings'].item()['data'].item()['directionDeg'].item()['n'].tolist()
stimDurMS = file['mapStimDurationMS'].item()['data'].tolist()
maxNumStim = 1000
hisPrePostMS = 50
numStim = zeros(1,numDir)
spikeCounts = zeros(numDir, maxNumStim)
spikeHists = zeros(numDir, stimDurMS + 2 * histPrePostMS)

