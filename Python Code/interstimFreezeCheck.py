import scipy.io as sp
import numpy as np

def interstimFreezeCheck(currTrial):
    '''
    takes the current trial and sees if the stimulus presentations have frozen
    interstim time / compares it to the expected interstim time
    '''

    frameRateHz = header['frameRateHz'].item()['data'].tolist()
    stimDur = header['stimSetting'].item()['data'].item()['stimDurationMS'].tolist()

    stimDesc = currTrial['stimDesc']
    for d in stimDesc.item #### cont. statement for frame on and frame off

    # does every trial have a stimDesc field/ sufficient to check 'frozen' 
    # interstim time only on valid trials (i.e. exclude instruction trials?)
