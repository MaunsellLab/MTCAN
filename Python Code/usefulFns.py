def loadMatFile(fileName):
    '''
    Loads the given matfile and assigns variables to access trial data

    Inputs: matfile name, str
    Outputs: variables, nd. array
    '''
    allTrials = sp.loadmat(fileName, squeeze_me = True)
    allTrialsData = allTrials['trials']
    header = allTrials['header']


def fieldInTrial(trial, fieldList):
    '''
    Function will check whether a field or all fields in the list is in the trial
    
    Inputs: trial (data struct from MATLAB)
            list of fields (list)
    Outputs: bool
    '''

    for field in fieldList:
        if field not in trial.dtype.names:
            return False
        else:
            return True

def targetOnsetIndexStimDesc(stimDesc):
    '''
    fn will identify the index of the target onset stimulus
    fn will put out ListTypes subdictionary within stimDesc 
    
    Inputs: stimDesc (variable)
    Outputs: the index of the target onset stim (int)
    '''

    for count, d in enumerate(stimDesc['listTypes']):
                if 2 in d:
                    targetOnsetStim = count
                    break
    
    return count