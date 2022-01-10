import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import scipy.io as sp
import math
import os
from collections import defaultdict
from usefulFns import *



def randTuningCurve(numNeurons):
    '''
    functon will generate random tuning cruves for x number of neurons

    Inputs:
        numNueurons (int): number of neurons 
    Outputs:
        tuningMat (2D array): matrix of tuning curve values for each
                              neuron
    '''
    
    tuningMat = np.zeros((numNeurons + 1, 6))
    tuningMat[0] = np.arange(0,360,60)

    for i in range(1, tuningMat.shape[0]):
        y_translate = np.random.randint(10,50)
        x_translate = np.random.randint(60,120)
        tuningMat[i,:] = np.sin((tuningMat[0,:] * np.pi / 180) + x_translate) + y_translate
    return tuningMat


tuningCurves = randTuningCurve(2)
