import scipy.io as sp
import numpy as np

allTrials = sp.loadmat('testa00.mat', squeeze_me = True)

file = allTrials['file']
mapSettings = file['mapSettings'].item()['data'].item()
directionDeg = mapSettings['directionDeg'].item()


