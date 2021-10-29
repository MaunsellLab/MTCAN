import scipy.io as sp
import numpy as np
import matplotlib.pyplot as plt
import math

allTrials = sp.loadmat('Meetz_2021_1028_1.mat', squeeze_me = True)
allTrialsData = allTrials['trials']

# code to check whether valid RT are being rewarded or not
# check EOT to make sure that trialEnd is not including 'wrong/distracted' 
correcttNoRwrd = []
for n,currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data']
    trial = currTrial['trial'].item()['data'].item()
    if extendedEOT == 6 and 'reactTimeMS' in currTrial.dtype.names:
        reactTimeMS = currTrial['reactTimeMS'].item()['data'].item()
        if reactTimeMS > 0 and reactTimeMS < 500:
            correcttNoRwrd.append(n)

# List of trials with early's but Target appears on screen
RTs = []
for n,currTrial in enumerate(allTrialsData.item()[0]):
    trial = currTrial['trial'].item()['data'].item()
    extendedEOT = currTrial['extendedEOT'].item()['data'].item()
    if extendedEOT == 6 and 'targetOn' in currTrial.dtype.names:
        RTs.append(n)

# trial RTs with early's but Target appears on screen
earlyWithRT = []
for n, currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data']
    if extendedEOT == 6 and 'targetOn' in currTrial.dtype.names and 'reactTimeMS' in currTrial.dtype.names:
        earlyWithRT.append((n, currTrial['reactTimeMS'].item()['data'].item()))

# Short trials with short RT
correctShortRT = []
for n, currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data']
    trial = currTrial['trial'].item()['data'].item()
    if extendedEOT == 0 and trial['catchTrial'] != 1:
        reactTimeMS = currTrial['reactTimeMS'].item()['data'].item()
        if reactTimeMS <100:
            correctShortRT.append((n, reactTimeMS))

# convert Left eyeXY to eyePos
eyePosDurTarget = []
for i in correctShortRT:
    trialNumber, reactTime = i
    currTrial = allTrialsData.item()[0][trialNumber]
    leftXDeg = []
    leftYDeg = []
    rightXDeg = []
    rightYDeg = []
    eyeLX = currTrial['eyeLXData'].item()['data'].item() 
    eyeLY = currTrial['eyeLYData'].item()['data'].item()
    eyeRX = currTrial['eyeRXData'].item()['data'].item()
    eyeRY = currTrial['eyeRYData'].item()['data'].item()
    eyeLeftCal = currTrial['eyeLeftCalibrationData'].item()['data'].item()['cal'].item()
    eyeRightCal = currTrial['eyeRightCalibrationData'].item()['data'].item()['cal'].item()

    count = min([len(eyeLX), len(eyeLY), len(eyeRX), len(eyeRY)])
    for l in range(0, count):
        xDegConvert = (eyeLX[l] * eyeLeftCal['m11'].item()) + (eyeLY[l] * eyeLeftCal['m21'].item()) + eyeLeftCal['tX'].item()
        leftXDeg.append(xDegConvert)
        yDegConvert = (eyeLX[l] * eyeLeftCal['m12'].item()) + (eyeLY[l] * eyeLeftCal['m22'].item()) + eyeLeftCal['tY'].item()
        leftYDeg.append(yDegConvert)
    # for r in range(0, count):
        xDegConvert = (eyeRX[l] * eyeRightCal['m11'].item()) + (eyeRY[l] * eyeRightCal['m21'].item()) + eyeRightCal['tX'].item()
        rightXDeg.append(xDegConvert)
        yDegConvert = (eyeRX[l] * eyeRightCal['m12'].item()) + (eyeRY[l] * eyeRightCal['m22'].item()) + eyeRightCal['tY'].item()
        rightYDeg.append(yDegConvert)

    fixate = currTrial['fixate'].item()['timeMS'].item()
    trialStart = currTrial['trialStart'].item()['timeMS'].item()
    fixateTimeWRTEyePos = math.floor((fixate - trialStart)/2)
    fixateTimeForIndex = fixateTimeWRTEyePos
    windowWidth = currTrial['fixWindowData'].item()['data'].item()['windowDeg'].item()['size'].item()['width'].item()
    windowHeight = currTrial['fixWindowData'].item()['data'].item()['windowDeg'].item()['size'].item()['height'].item()
    for index, xDeg in enumerate(leftXDeg[fixateTimeWRTEyePos:]):
        if -(windowWidth/2) > xDeg or xDeg > (windowWidth/2) or -(windowHeight/2) > leftYDeg[fixateTimeForIndex] or leftYDeg[fixateTimeForIndex] > (windowHeight/2):
            eyePosSac = index + fixateTimeWRTEyePos
            break
        else:
            fixateTimeForIndex += 1
    eyePosSacTimeMS = eyePosSac * 2
    saccadeTime = currTrial['saccade'].item()['timeMS'].item()
    diff = (saccadeTime - trialStart) - eyePosSacTimeMS
    
    targetOn = currTrial['targetOn'].item()['timeMS'].item()
    eyePosReactTime = eyePosSacTimeMS - (targetOn - trialStart)
    eyePosDurTarget.append((trialNumber, eyePosReactTime))
    knotRTSaccTarget = saccadeTime - targetOn
    targetOnTimeEyePos = (targetOn - trialStart)/2

    plt.figure()
    plt.plot(leftXDeg, color = 'olive')
    plt.plot(leftYDeg, color = 'green')
    plt.plot(rightXDeg, color = 'blue')
    plt.plot(rightYDeg, color = 'red')
    plt.axvline(x = targetOnTimeEyePos)
    plt.title(eyePosReactTime)

# scatter plot RTs aligned with targetOnset
for y,x in eyePosDurTarget:
    plt.scatter(x,y) 
plt.axvline(x=0)

# scatter plot of trialOutcomes 
for n,outcome in enumerate(trialOutcomes):
    if outcome == 0:
        color = 'green'
    elif outcome == 6:
        color = 'orange'
    else:
        color = 'blue'
    plt.scatter(n, outcome, c = color, s = 0.5)
for x_axis in catchTrials:
    plt.axvline(x = x_axis, color = 'pink')
for x_axis in invalidTrials:
    plt.axvline(x = x_axis, color = 'black')
plt.show()

# diff in the len of the stimOn and stimOff list
stimDiff = []
for currTrial in allTrialsData.item()[0]:
    trial = currTrial['trial'].item()['data'].item()
    if trial['instructTrial'] != 1:
        if 'stimulusOn' in currTrial.dtype.names and 'stimulusOff' in currTrial.dtype.names:
            stimulusOn = currTrial['stimulusOn'].item()['timeMS'].item()
            stimulusOff = currTrial['stimulusOff'].item()['timeMS'].item()
            if type(stimulusOn) == int:
                stimDiff.append(0)
            elif type(stimulusOff) == int:
                diff = len(stimulusOn) - 1
                stimDiff.append(diff)
            else:
                diff = len(stimulusOn) - len(stimulusOff)
                stimDiff.append(diff)

dist = {}
for i in stimDiff:
    dist[i] = dist.get(i,0)+1

#hist
plt.hist(np.array(distRT).flatten(), bins = 20)

#weird Trials
weirdTrials = []
for n, currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data']
    if extendedEOT == 6 and 'reactTimeMS' in currTrial.dtype.names:
        weirdTrials.append(n)

#stim diff for trials preceding the weird trials
stimDiff = []
for i in weirdTrials:
    currTrial = allTrialsData.item()[0][i-1]
    trial = currTrial['trial'].item()['data'].item()
    if trial['instructTrial'] != 1:
        if 'stimulusOn' in currTrial.dtype.names and 'stimulusOff' in currTrial.dtype.names:
            stimulusOn = currTrial['stimulusOn'].item()['timeMS'].item()
            stimulusOff = currTrial['stimulusOff'].item()['timeMS'].item()
            if type(stimulusOn) == int:
                stimDiff.append(0)
            elif type(stimulusOff) == int:
                diff = len(stimulusOn) - 1
                stimDiff.append(diff)
            else:
                diff = len(stimulusOn) - len(stimulusOff)
                stimDiff.append(diff)
