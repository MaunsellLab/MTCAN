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

            