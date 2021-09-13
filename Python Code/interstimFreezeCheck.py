import scipy.io as sp
import numpy as np

# when importing .mat file, make sure it is the combined file with both the 
# header and all trials saved 

allTrials = sp.loadmat('alltrials_testing_2021_0912.mat', squeeze_me = True)
allTrialsData = allTrials['trials']
header = allTrials['header']

frameRateHz = header['frameRateHz'].item()['data'].tolist()
stimDuration = header['stimSetting'].item()['data'].item()['stimDurationMS'].tolist()


for i in range(0,len(allTrialsData.item()[0])-1):
    currTrial = allTrialsData.item()[0][i]

    # does every trial have a stimDesc field/ sufficient to check 'frozen' 
    # interstim time only on valid trials (i.e. exclude instruction trials?)
