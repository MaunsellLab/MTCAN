import scipy.io as sp
import numpy as np
import matplotlib.pyplot as plt
import math
import os
from collections import defaultdict
from usefulFns import *

allTrialsData, header = loadMatFile('Meetz_2021_1028_1.mat')



# code to check whether valid RT are being rewarded or not
# check EOT to make sure that trialEnd is not including 'wrong/distracted' 
correcttNoRwrd = []
for n,currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data']
    trial = currTrial['trial'].item()['data'].item()
    if extendedEOT == 6 and fieldInTrial(['reactTimeMS'], currTrial):
        reactTimeMS = currTrial['reactTimeMS'].item()['data'].item()
        if reactTimeMS > 0 and reactTimeMS < 500:
            correcttNoRwrd.append(n)

# List of trials with early's but Target appears on screen
RTs = []
for n,currTrial in enumerate(allTrialsData.item()[0]):
    trial = currTrial['trial'].item()['data'].item()
    extendedEOT = currTrial['extendedEOT'].item()['data'].item()
    if extendedEOT == 6 and fieldInTrial(['targetOn'], currTrial):
        RTs.append(n)

# trial RTs with early's but Target appears on screen
earlyWithRT = []
for n, currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data']
    if extendedEOT == 6 and fieldInTrial(['targetOn', 'reactTimeMS'], currTrial):
        earlyWithRT.append((n, currTrial['reactTimeMS'].item()['data'].item()))

# correct trials with short RT
correctShortRT = []
for n, currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data']
    trial = currTrial['trial'].item()['data'].item()
    if extendedEOT == 0 and trial['catchTrial'] != 1:
        reactTimeMS = currTrial['reactTimeMS'].item()['data'].item()
        if reactTimeMS <100:
            correctShortRT.append((n, reactTimeMS))

# EyePos for a trial
currTrial = allTrialsData.item()[0][i]
eyesXYDeg = defaultdict(list)
eyeLX = currTrial['eyeLXData'].item()['data'].item() 
eyeLY = currTrial['eyeLYData'].item()['data'].item()
eyeRX = currTrial['eyeRXData'].item()['data'].item()
eyeRY = currTrial['eyeRYData'].item()['data'].item()
eyeLeftCal = currTrial['eyeLeftCalibrationData'].item()['data'].item()['cal'].item()
eyeRightCal = currTrial['eyeRightCalibrationData'].item()['data'].item()['cal'].item()
count = min([len(eyeLX), len(eyeLY), len(eyeRX), len(eyeRY)])

for s in range(0, count):
    xDegConvert = (eyeLX[s] * eyeLeftCal['m11'].item()) + (eyeLY[s] * eyeLeftCal['m21'].item()) + eyeLeftCal['tX'].item()
    eyesXYDeg['leftX'].append(xDegConvert)
    yDegConvert = (eyeLX[s] * eyeLeftCal['m12'].item()) + (eyeLY[s] * eyeLeftCal['m22'].item()) + eyeLeftCal['tY'].item()
    eyesXYDeg['leftY'].append(yDegConvert)
    xDegConvert = (eyeRX[s] * eyeRightCal['m11'].item()) + (eyeRY[s] * eyeRightCal['m21'].item()) + eyeRightCal['tX'].item()
    eyesXYDeg['rightX'].append(xDegConvert)
    yDegConvert = (eyeRX[s] * eyeRightCal['m12'].item()) + (eyeRY[s] * eyeRightCal['m22'].item()) + eyeRightCal['tY'].item()
    eyesXYDeg['rightY'].append(yDegConvert)

# EyePos for a list of trials indexes
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

'''
OR
'''
for n,currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data'].item()
    if extendedEOT == 0:
        plt.scatter(n+1,extendedEOT, color = 'green')
    elif extendedEOT == 6:
        plt.scatter(n+1, extendedEOT, color = 'orange')
    else:
        plt.scatter(n+1, extendedEOT, color = 'blue')
plt.show()

# diff in the len of the stimOn and stimOff list
stimDiff = []
for currTrial in allTrialsData.item()[0]:
    trial = currTrial['trial'].item()['data'].item()
    if trial['instructTrial'] != 1:
        if fieldInTrial(['stimulusOn', 'stimulusOff'], currTrial):   
            stimulusOn = currTrial['stimulusOn'].item()['timeMS'].item()
            stimulusOff = currTrial['stimulusOff'].item()['timeMS'].item()
            if type(stimulusOn) == int:
                stimDiff.append(0)
            elif type(stimulusOff) == int:
                stimDiff.append(len(stimulusOn) - 1)
            else:
                stimDiff.append(len(stimulusOn) - len(stimulusOff))

dist = {}
for i in stimDiff:
    dist[i] = dist.get(i,0)+1

#hist
plt.hist(np.array(distRT).flatten(), bins = 20)

#weird Trials
weirdTrials = []
for n, currTrial in enumerate(allTrialsData.item()[0]):
    extendedEOT = currTrial['extendedEOT'].item()['data']
    if extendedEOT == 6 and fieldInTrial(['reactTimeMS'], currTrial):
        weirdTrials.append(n)

#stim diff for trials preceding the weird trials
stimDiff = []
for i in weirdTrials:
    currTrial = allTrialsData.item()[0][i-1]
    trial = currTrial['trial'].item()['data'].item()
    if trial['instructTrial'] != 1:
        if fieldInTrial(['stimulusOn', 'stimulusOff'], currTrial):  
            stimulusOn = currTrial['stimulusOn'].item()['timeMS'].item()
            stimulusOff = currTrial['stimulusOff'].item()['timeMS'].item()
            if type(stimulusOn) == int:
                stimDiff.append(0)
            elif type(stimulusOff) == int:
                stimDiff.append(len(stimulusOn) - 1)
            else:
                stimDiff.append(len(stimulusOn) - len(stimulusOff))


'''
MTNAN
'''

#code will generate a list of interstim frame diff for task gabor
stimDescInterstim = []
for count, currTrial in enumerate(allTrialsData.item()[0]):
    print(count)
    trial = currTrial['trial'].item()['data'].item()
    stimDesc = currTrial['stimDesc'].item()['data'].item()
    taskStimCount = 0
    for stim in stimDesc:
        if stim['stimLoc'] == 2:
            if taskStimCount == 0:
                frameOff = stim['stimOffFrame']
                taskStimCount += 1
            else:
                frameDiff = stim['stimOnFrame'] - frameOff
                frameOff = stim['stimOffFrame']
                stimDescInterstim.append(frameDiff)



targetOnTimes = []

for currTrial in allTrialsData.item()[0]:
    trial = currTrial['trial'].item()['data'].item()
    if trial['instructTrial'] != 1 and trial['catchTrial'] != 1:
        targetOnTimes.append(trial['targetOnTimeMS'].tolist())
    




# for count, i in enumerate(stimDesc['listType']):
#     if i == 2:
#         stimSeqLen = count
#         break
# targetOnFrame = stimDesc['stimOnFrame'][stimSeqLen]
# targetOnTimeMs = targetOnFrame * (1000/75)
# targetOnTimes.append(targetOnTimeMS)