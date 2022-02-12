
import mat73

def loadMatFile(fileName):
    '''
    Loads the given matfile and assigns variables to access trial data

    Inputs: matfile name, (str)
    Outputs: variables, (nd.array)
    '''
    
    allTrials = mat73.loadmat(f'../Matlab Data/{fileName}', use_attrdict = True)
    allTrialsData = allTrials.trials
    header = allTrials.header

    return allTrialsData, header