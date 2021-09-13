import scipy.io as sp
import numpy as np
from collections import defaultdict ## All edits to incorporate this marked w/ comments
import argparse

parser = argparse.ArgumentParser('')
parser.add_argument('data_fn', type=str)
parser.add_argument('n_stimuli', type=int, default=9)

# Read in command line args
args = parser.parse_args()

def stimSeqCheck(attendLoc, stimDesc, stimSequence, prevStim):
    '''
    Function to create a stimulus seq based on attended location.
    Function will also ensure that subsequence trials follow the stimulus
    sequence that was generated. 

    params: attendLoc    (int)        -- RF ID to attend
            stimDesc     (np.ndarray) -- stimulus description
            stimSequence (list)       -- list that indicates stimulus ordering
            prevStim     (dict)       -- dict of prev stim by location    
    returns: locDict, stimSequence, prevStim
    '''
    for d in stimDesc['stimTypes'][1:targetOnsetStim]:

        stimulus = d[2:4].tolist()

        # If this is the first time (e.g. all 9 haven't been seen), 
        # add it to the sequence
        if len(stimSequence) < args.n_stimuli:
            stimSequence.append(stimulus)
        else:
            # Assert that the currrent stimulus at this location
            # is in the correct order
            prev = stimSequence.index(prevStim[attendLoc])
            cur = stimSequence.index(stimulus)
            assert (cur % args.n_stimuli) == ((prev + 1) % 9), \
                f"Violation of stimulus sequence: prev.= {prev}, cur.={cur}"

            lastStim[attendLoc] = currStim[attendLoc]

        # Increment the count for this (location, stimulus) pair by 1
        locDict[(attendLoc, stimulus)] += 1

        # Record that this was the last seen stimulus at the current location
        prevStim[attendLoc] = stimulus

    return locDict, stimSequence, prevStim

if __name__ == "__main__":

    # Import mat file with multiple trials
    allTrials = sp.loadmat(args.data_fn, squeeze_me = True)
    allTrialsData = allTrials['trials'].item()[0]

    # Make structures for storing stim sequence + counts 
    locDict = defaultdict(int)
    stimSequence = []
    prevStim = {}

    # Run through diff trials 
    for trial in allTrialsData:
        trial_id  = trial['trial'].item()['data']
        trial_end = trial['trialEnd'].item()['data']

        if trial_end == 0: 
            if trial_id['instructTrial'] == 0:
                stimDesc = trial['stimDesc'].item()['data'].item()
                for count, d in enumerate(stimDesc['listTypes']):
                    if 2 in d:
                        targetOnsetStim = count
                        break # maybe continue? if we want to keep looping through the remaining trials
                attendLoc = trial_id['attendLoc'].tolist()
                locDict, stimSequence, prevStim = stimSeqCheck(attendLoc, stimDesc, stimSequence, prevStim)
