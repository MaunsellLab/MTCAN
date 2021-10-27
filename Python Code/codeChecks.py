import scipy.io as sp
import numpy as np
import matplotlib.pyplot as plt
import math

allTrials = sp.loadmat('Meetz_2021_1021_4.mat', squeeze_me = True)
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
        earlyWithRT.append(currTrial['reactTimeMS'].item()['data'].item())

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
    for r in range(0, count):
        xDegConvert = (eyeRX[r] * eyeRightCal['m11'].item()) + (eyeRY[r] * eyeRightCal['m21'].item()) + eyeRightCal['tX'].item()
        rightXDeg.append(xDegConvert)
        yDegConvert = (eyeRX[r] * eyeRightCal['m12'].item()) + (eyeRY[r] * eyeRightCal['m22'].item()) + eyeRightCal['tY'].item()
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
    sacFromTrialStart = currTrial['saccade'].item()['timeMS'].item() - trialStart
    diff = sacFromTrialStart - eyePosSacTimeMS
    print(trialNumber, diff)

    targetOn = currTrial['targetOn'].item()['timeMS'].item()
    trialStart = currTrial['trialStart'].item()['timeMS'].item()
    targetOnTimeEyePos = (targetOn - trialStart)/2

    plt.figure()
    plt.plot(leftXDeg, color = 'olive')
    plt.plot(leftYDeg, color = 'green')
    plt.plot(rightXDeg, color = 'blue')
    plt.plot(rightYDeg, color = 'red')
    plt.axvline(x = targetOnTimeEyePos)
    plt.title(reactTime)