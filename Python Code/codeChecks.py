
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

# convert Left eyeXY to eyePos
eyeLX = currTrial['eyeLXData'].item()['data'].item() 
eyeLY = currTrial['eyeLYData'].item()['data'].item()
eyeRX = currTrial['eyeRXData'].item()['data'].item()
eyeRY = currTrial['eyeRYData'].item()['data'].item()
eyeLeftCal = currTrial['eyeLeftCalibrationData'].item()['data'].itm()['cal'].item()
eyeRightCal = currTrial['eyeRightCalibrationData'].item()['data'].item()['cal'].item()

for n,eyePos in enumerate(eyeLX):
    xDeg = (eyePos * eyeLeftCal['m11'].item()) + (eyeLY[n] * eyeLeftCal['m21'].item()) + eyeLeftCal['tX'].item()
    yDeg = (eyePos * eyeLeftCal['m12'].item()) + (eyeLY[n]) * eyeLeftCal['m22'].item() + eyeLeftCal['tY'].item()

for n, eyePos in enumerate(eyeRX):
    xDeg = (eyePos * eyeRightCal['m11'].item()) + (eyeRY[n]) * eyeRightCal['m21'].item()) + eyeRightCal['tX'].item()
    yDeg = (eyePos * eyeRightCal['m12'].item()) + (eyeRY[n]) * eyeRightCal['m22'].item()) + eyeRightCal['tY'].item()


