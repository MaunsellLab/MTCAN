'''
MTN Analysis Script
 - heatmap of Norm responses
 - potentially some correlations
'''

import numpy as np
from usefulFns import *

# load my file 

def activeUnits(unitData):

    '''
    function returns the active units across trials for a session as a list

    Inputs: unitData (str): spikeData
    Outputs: units (list): active units for a sessioon

    '''
    units = []
    for currTrial in allTrials:
        if unitData in currTrial:
            uniqueUnits = np.unique(currTrial[unitData]['unit'])
            for unique in uniqueUnits:
                if unique not in units:
                    units.append(unique)
    
    return units


# generate list of unique active units
units = activeUnits('spikeData')

spikeCountMat = np.zeros((numNeurons,30,169))
spikeCountMat[:,0,:] = np.arange(0,169)
spikeCountMat[:,1:,:] = np.nan


for currTrial in allTrials: 
    extendedEOT = currTrial['extendedEOT']['data']
    trial = currTrial['trial']['data']
    if trial['instructTrial'] != 1 and 'spikeData' in currTrial:


