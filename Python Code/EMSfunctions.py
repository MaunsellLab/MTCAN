def emsNormGenCondensed0(fixed, LN, LP, S, al0, al1, sig, m):
    """
    scaled loc is 1
    loc 0 has stronger response
    """

    c0, c1, l0, l1 = fixed
    L = np.array([LN, LP])

    loc0 = (c0 * (l0 * L)).sum(-1) / (
            c0[:, 0] + (c1[:, 0] * al1) + sig)

    loc1 = (c1 * S * (l1 * L)).sum(-1) / (
           (c0[:, 0] * al0) + c1[:, 0] + sig)

    return loc0 + loc1 + m


def emsNormGenCondensed1(fixed, LN, LP, S, al0, al1, sig, m):
    """
    scaled loc is 1
    loc 0 has stronger response
    """

    c0, c1, l0, l1 = fixed
    L = np.array([LN, LP])

    loc0 = (c0 * S * (l0 * L)).sum(-1) / (
            c0[:, 0] + (c1[:, 0] * al1) + sig)

    loc1 = (c1 * (l1 * L)).sum(-1) / (
           (c0[:, 0] * al0) + c1[:, 0] + sig)

    return loc0 + loc1 + m


def fixedValsForEMSGenCondensed(stimMatReIndex, stimIndexDict, prefDir, nullDir):
    '''
    function will create the fixed values for running the
    generic EMS normazliation curve_fit
    '''

    c0s, c1s, l0s, l1s = [], [], [], []
    direction_set = np.arange(0, 360, 60)
    directionSet = np.array([nullDir, prefDir])
    if type(stimMatReIndex) == int:
        c0 = stimIndexDict[stimMatReIndex][0]['contrast']
        l0 = stimIndexDict[stimMatReIndex][0]['direction']
        c1 = stimIndexDict[stimMatReIndex][1]['contrast']
        l1 = stimIndexDict[stimMatReIndex][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(2)
        l1_oh = np.zeros(2)
        l0_oh[np.argwhere(directionSet == l0).squeeze()] = 1
        l1_oh[np.argwhere(directionSet == l1).squeeze()] = 1
        c0s.append(np.repeat(c0, 2))
        c1s.append(np.repeat(c1, 2))
        l0s.append(l0_oh)
        l1s.append(l1_oh)
    else:
        for i in stimMatReIndex:
            c0 = stimIndexDict[i][0]['contrast']
            l0 = stimIndexDict[i][0]['direction']
            c1 = stimIndexDict[i][1]['contrast']
            l1 = stimIndexDict[i][1]['direction']

            # Make one-hot encoding of l0 and l1
            l0_oh = np.zeros(2)
            l1_oh = np.zeros(2)
            l0_oh[np.argwhere(directionSet == l0).squeeze()] = 1
            l1_oh[np.argwhere(directionSet == l1).squeeze()] = 1
            c0s.append(np.repeat(c0, 2))
            c1s.append(np.repeat(c1, 2))
            l0s.append(l0_oh)
            l1s.append(l1_oh)

    return (np.array(c0s), np.array(c1s), np.array(l0s), np.array(l1s))



def emsNormGen0(fixed, L0, L60, L120, L180, L240, L300, S, al0, al1, sig, m):
    '''
    scaled loc is 1
    loc 0 has stronger response
    '''
    c0, c1, l0, l1 = fixed
    L = np.array([L0, L60, L120, L180, L240, L300])
    loc0 = (c0 * (l0 * L)).sum(-1) / (
            c0[:, 0] + (c1[:, 0] * al1) + sig)

    loc1 = (c1 * S * (l1 * L)).sum(-1) / (
           (c0[:, 0] * al0) + c1[:, 0] + sig)

    return loc0 + loc1 + m


def emsNormGen1(fixed, L0, L60, L120, L180, L240, L300, S, al0, al1, sig, m):
    '''
    scaled loc is 0
    loc 1 has stronger response
    '''

    c0, c1, l0, l1 = fixed
    L = np.array([L0, L60, L120, L180, L240, L300])
    loc0 = (c0 * S * (l0 * L)).sum(-1) / (
            c0[:, 0] + (c1[:, 0] * al1) + sig)

    loc1 = (c1 * (l1 * L)).sum(-1) / (
           (c0[:, 0] * al0) + c1[:, 0] + sig)

    return loc0 + loc1 + m


def fixedValsForEMSGen(stimMatReIndex, stimIndexDict):
    '''
    function will create the fixed values for running the
    generic EMS normazliation curve_fit
    '''

    c0s, c1s, l0s, l1s = [], [], [], []
    direction_set = np.arange(0, 360, 60)
    for i in stimMatReIndex.reshape(49):
        c0 = stimIndexDict[i][0]['contrast']
        l0 = stimIndexDict[i][0]['direction']
        c1 = stimIndexDict[i][1]['contrast']
        l1 = stimIndexDict[i][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(6)
        l1_oh = np.zeros(6)
        l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
        l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
        c0s.append(np.repeat(c0,6))
        c1s.append(np.repeat(c1,6))
        l0s.append(l0_oh)
        l1s.append(l1_oh)

    return (np.array(c0s), np.array(c1s), np.array(l0s), np.array(l1s))


r2 = []
for unitCount, unit in enumerate(units):
    b = meanSpikeReshaped[unitCount].reshape(7,7) * 1000/trueStimDurMS
    bReIndex = np.zeros((7,7))

    # find direction tested that is closest to the pref dir
    # and reindex around this so that it is in the middle of the grid
    prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
    nullDirIndex = np.where(dirArray == nullDir)[0][0]
    reIndex = (np.array([0,1,2,3,4,5])+nullDirIndex) % 6

    tempMain = b[:6, :6][:, reIndex]
    tempMain = tempMain[:6, :6][reIndex, :]
    temp0Blank = b[:6, 6][reIndex]
    temp1Blank = b[6, :6][reIndex]
    bReIndex[:6, :6] = tempMain
    bReIndex[:6, 6] = temp0Blank
    bReIndex[6, :6] = temp1Blank
    bReIndex[6, 6] = b[6, 6]

    # fixed (independent) variables - matrix of corresponding stim Indexes
    stimMat = np.zeros((7, 7))
    stimMat[:6, :6] = np.arange(36).reshape(6, 6)
    stimMat[6, :6] = np.arange(36, 42)
    stimMat[:, 6] = np.arange(42, 49)

    # reshape fixed variables to match dependent variables
    stimMatReIndex = np.zeros((7, 7))
    tempMain = stimMat[:6, :6][:,reIndex]
    tempMain = np.squeeze(tempMain[:6, :6][reIndex, :])
    temp0Blank = stimMat[:6, 6][reIndex]
    temp1Blank = stimMat[6, :6][reIndex]
    stimMatReIndex[:6, :6] = tempMain
    stimMatReIndex[:6, 6] = temp0Blank
    stimMatReIndex[6, :6] = temp1Blank
    stimMatReIndex[6, 6] = stimMat[6, 6]

    ## for fitting across mean spike count
    resp = bReIndex.reshape(49)
    fixedVals = fixedValsForEMSGen(stimMatReIndex, stimIndexDict)
    if max(bReIndex[:6, 6])-min(bReIndex[:6, 6]) > max(bReIndex[6, :6])-min(bReIndex[6, :6]):
        pOpt, pCov = curve_fit(emsNormGen1, fixedVals, resp, bounds=(
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 1, 10, 10, 1, np.inf)))
        print(unit, pOpt)
        y_pred = emsNormGen1(fixedVals, *pOpt)
        print(r2_score(resp.squeeze(), y_pred))
    else:
        print('using norm func 0')
        pOpt, pCov = curve_fit(emsNormGen0, fixedVals, resp, bounds=(
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 1, 10, 10, 1, np.inf)))
        print(unit, pOpt)
        y_pred = emsNormGen0(fixedVals, *pOpt)
        print(r2_score(resp.squeeze(), y_pred))

    r2.append(r2_score(resp.squeeze(), y_pred))

r2 = np.array(r2)
np.mean(r2)



# step 1: fitting gaussian function to single presentations
# to extract gaussian fit and scaling factor (scalar on loc0)
# if scalar >1 then loc0 has stronger resp
resp = np.concatenate((bReIndex[6, :6], bReIndex[:6, 6]), axis=0)
singleStim = np.concatenate((stimMatReIndex[6, :6], stimMatReIndex[:6, 6]),
                            axis=0)
fixedVals = fixedValsForCurveFit(prefDir, singleStim, stimIndexDict)

pOpt, pCov = curve_fit(gaussNormFunc, fixedVals, resp.squeeze(),
                       bounds=((0, 0, 0, 0, 0),
                               (np.inf, np.inf, 360, 360, np.inf)))

y_pred = gaussNormFunc(fixedVals, *pOpt)
print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')

scalar = pOpt[4]

# step 2, fitting normalization eqn with fixed gauss tuning curve
# only two free parameters (a0, a1)

resp = bReIndex[:6, :6].reshape(36)
pairedStim = stimMatReIndex[:6, :6].reshape(36)
fixedVals = fixedValsCurveFitForPairedStim(prefDir, pairedStim, stimIndexDict, pOpt)
pOpt, pCov = curve_fit(EMSNormStepFunc, fixedVals, resp.squeeze(),
                       bounds=((-10, -10),
                               (10, 10)))

y_pred = EMSNormStepFunc(fixedVals, *pOpt)
print(r2_score(resp.squeeze(), y_pred))
r2 = r2_score(resp.squeeze(), y_pred)
print(pOpt)

# step 3, fit normalization eqn with fixed guass tuning curve to
# p+n, n+p conditions only, only two free parameters (a0,a1)
condensedMat = np.zeros((2, 2))
condensedMat[0, :] = bReIndex[0, [0, 3]]
condensedMat[1, :] = bReIndex[3, [0, 3]]
resp = condensedMat.reshape(4)

condensedStimMat = np.zeros((2, 2))
condensedStimMat[0, :] = stimMatReIndex[0, [0, 3]]
condensedStimMat[1, :] = stimMatReIndex[3, [0, 3]]
fixedVals = fixedValsCurveFitForPairedStim(prefDir, condensedStimMat.reshape(4),
                                           stimIndexDict, pOpt)

if max(bReIndex[:6, 6]) - min(bReIndex[:6, 6]) > max(bReIndex[6, :6]) - min(bReIndex[6, :6]):
    pOpt, pCov = curve_fit(EMSNormStepFunc, fixedVals, resp.squeeze(), bounds=(
        (-100, -100), (100, 100)))
    print(pOpt)
    y_pred = EMSNormStepFunc(fixedVals, *pOpt)
    print(r2_score(resp.squeeze(), y_pred))
    r2 = r2_score(resp.squeeze(), y_pred)
    sLoc = 0
else:
    print('using norm func 0')
    pOpt, pCov = curve_fit(EMSNormStepFunc, fixedVals, resp.squeeze(), bounds=(
        (-100, -100), (100, 100)))
    print(pOpt)
    y_pred = EMSNormStepFunc(fixedVals, *pOpt)
    print(r2_score(resp.squeeze(), y_pred))
    r2 = r2_score(resp.squeeze(), y_pred)
    sLoc = 1

unitAlphaIndex[unitCount] = (pOpt[0] + pOpt[1]) / 2

# step 1: fitting L0-L6 to single presentations (scalar on loc0)
# if scalar >1 then loc0 has stronger resp than loc1
resp = np.concatenate((bReIndex[6, :6], bReIndex[:6, 6]), axis=0)
singleStim = np.concatenate((stimMatReIndex[6, :6], stimMatReIndex[:6, 6]),
                            axis=0)
fixedVals = fixedValsForEMSGen(singleStim, stimIndexDict)
pOpt, pCov = curve_fit(emsSingleStim, fixedVals, resp,
                       bounds=((0, 0, 0, 0, 0, 0, 0),
                               (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
y_pred = emsSingleStim(fixedVals, *pOpt)
print(f'the L0-L6 version fit R2 {r2_score(resp.squeeze(), y_pred)}')

# step 2: fitting normalization eqn with fixed L0-L6
# only two free parameters (a0,a1)
resp = bReIndex[:6, :6].reshape(36)
pairedStim = stimMatReIndex[:6, :6].reshape(36)
fixedVals = fixedValsForPairedStimL0L6(pairedStim, stimIndexDict, pOpt)
pOpt, pCov = curve_fit(EMSNormStepFuncGen, fixedVals, resp, bounds=(
    (-np.inf, -np.inf), (np.inf, np.inf)))
y_pred = EMSNormStepFuncGen(fixedVals, *pOpt)
print(f'EMS L0-L300 version fit {r2_score(resp.squeeze(), y_pred)}')
print(pOpt)

if scalar <= 1:
    sLoc = 0
    unitAlphaIndex[unitCount] = pOpt[0] / (pOpt[0] + pOpt[1])
else:
    sLoc = 1
    unitAlphaIndex[unitCount] = pOpt[1] / (pOpt[0] + pOpt[1])
scaledLoc[unitCount] = sLoc
unitR2[unitCount] = r2

"""
# BO, A, MU, SIG, S, al0, al1, c50, m


# using EMS Norm with Gauss Tuning

# fixedVals = fixedValsForCurveFit(prefDir, stimMatReIndex.reshape(49), stimIndexDict)
# if max(bReIndex[:6, 6])-min(bReIndex[:6, 6]) > max(bReIndex[6, :6])-min(bReIndex[6, :6]):
#     pOpt, pCov = curve_fit(emsNormFunc1, fixedVals, resp.squeeze(),
#                         bounds=((0, 0, prefDir-15, 0, 0, 0, 0, 0, 0),
#                         (np.inf, max(bReIndex[:6, 6]), prefDir+15, 360, 1, 10, 10, 1, 100)))
#     print(unit,pOpt)
#     y_pred = emsNormFunc1(fixedVals, *pOpt)
#     print(r2_score(resp.squeeze(), y_pred))
#     r2 = r2_score(resp.squeeze(), y_pred)
#     sLoc = 0
# else:
#     pOpt, pCov = curve_fit(emsNormFunc0, fixedVals, resp.squeeze(),
#                         bounds=((0, 0, prefDir-15, 0, 0, 0, 0, 0, 0),
#                         (np.inf, max(bReIndex[6, :6]), prefDir+15, 360, 1, 10, 10, 1, 100)))
#     print('using normFunc0')
#     print(unit,pOpt)
#     y_pred = emsNormFunc0(fixedVals, *pOpt)
#     print(r2_score(resp.squeeze(), y_pred))
#     r2 = r2_score(resp.squeeze(), y_pred)
#     sLoc = 1


# ## fitting just the pref+null and diff combs of that

# condensedMat = np.zeros((3, 3))
# condensedMat[0, :] = bReIndex[0, [0, 3, 6]]
# condensedMat[1, :] = bReIndex[3, [0, 3, 6]]
# condensedMat[2, :] = bReIndex[6, [0, 3, 6]]
# resp = condensedMat.reshape(9)
#
# condensedStimMat = np.zeros((3, 3))
# condensedStimMat[0, :] = stimMatReIndex[0, [0, 3, 6]]
# condensedStimMat[1, :] = stimMatReIndex[3, [0, 3, 6]]
# condensedStimMat[2, :] = stimMatReIndex[6, [0, 3, 6]]
# fixedVals = fixedValsForEMSGenCondensed(condensedStimMat.reshape(9),
#                                         stimIndexDict, prefDir, nullDir)
#
# if max(bReIndex[:6, 6])-min(bReIndex[:6, 6]) > max(bReIndex[6, :6])-min(bReIndex[6, :6]):
#     pOpt, pCov = curve_fit(emsNormGenCondensed1, fixedVals, resp, bounds=(
#         (0, 0, 0, 0, 0, 0, 0),
#         (np.inf, np.inf, 1, 10, 10, 1, np.inf)))
#     print(unit, pOpt)
#     y_pred = emsNormGenCondensed1(fixedVals, *pOpt)
#     print(r2_score(resp.squeeze(), y_pred))
#     r2 = r2_score(resp.squeeze(), y_pred)
#     sLoc = 0
# else:
#     print('using norm func 0')
#     pOpt, pCov = curve_fit(emsNormGenCondensed0, fixedVals, resp, bounds=(
#         (0, 0, 0, 0, 0, 0, 0),
#         (np.inf, np.inf, 1, 10, 10, 1, np.inf)))
#     print(unit, pOpt)
#     y_pred = emsNormGenCondensed0(fixedVals, *pOpt)
#     print(r2_score(resp.squeeze(), y_pred))
#     r2 = r2_score(resp.squeeze(), y_pred)
#     sLoc = 1
#
# scaledLoc[unitCount] = sLoc
# unitNormFitEstimate[unitCount] = pOpt
# unitR2[unitCount] = r2
# unitAlphaMax[unitCount] = max([pOpt[3], pOpt[4]])
# unitAlphaMin[unitCount] = min([pOpt[3], pOpt[4]])


# # using EMS Norm with L0-L6
#
# ## for fitting across mean spike count
# resp = bReIndex.reshape(49)
# fixedVals = fixedValsForEMSGen(stimMatReIndex.reshape(49), stimIndexDict)
#
# if max(bReIndex[:6, 6])-min(bReIndex[:6, 6]) > max(bReIndex[6, :6])-min(bReIndex[6, :6]):
#     pOpt, pCov = curve_fit(emsNormGen1, fixedVals, resp, bounds=(
#         (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#         (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 10, 10, 1, np.inf)))
#     print(unit, pOpt)
#     y_pred = emsNormGen1(fixedVals, *pOpt)
#     print(r2_score(resp.squeeze(), y_pred))
#     r2 = r2_score(resp.squeeze(), y_pred)
#     sLoc = 0
# else:
#     print('using norm func 0')
#     pOpt, pCov = curve_fit(emsNormGen0, fixedVals, resp, bounds=(
#         (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#         (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 10, 10, 1, np.inf)))
#     print(unit, pOpt)
#     y_pred = emsNormGen0(fixedVals, *pOpt)
#     print(r2_score(resp.squeeze(), y_pred))
#     r2 = r2_score(resp.squeeze(), y_pred)
#     sLoc = 1


#
# fixedVals = []
# for i in stimMatReIndex.reshape(49):
#     c0 = stimIndexDict[i][0]['contrast']
#     l0 = stimIndexDict[i][0]['direction']
#     if abs(l0-prefDir) > 180:
#         l0 = prefDir - (360 - (l0 - prefDir))
#     c1 = stimIndexDict[i][1]['contrast']
#     l1 = stimIndexDict[i][1]['direction']
#     if abs(l1-prefDir) > 180:
#         l1 = prefDir - (360 - (l1 - prefDir))
#     fixedVals.append((c0,l0,c1,l1))
# fixedVals = np.array(fixedVals)

# # for fitting across every trial spike count
# fixedVals = []
# for i in np.tile(stimMatReIndex.reshape(49),blocksDone):
#     c0 = stimIndexDict[i][0]['contrast']
#     l0 = stimIndexDict[i][0]['direction']
#     c1 = stimIndexDict[i][1]['contrast']
#     l1 = stimIndexDict[i][1]['direction']
#     fixedVals.append((c0,l0,c1,l1))
# fixedVals = np.array(fixedVals)
# tempMat = spikeCountMat[unitCount,:blocksDone,:][:,stimMatReIndex.reshape(49).astype(int)]
# resp = tempMat.reshape(49*blocksDone) * 1000/trueStimDurMS
#
# if max(bReIndex[:6,6])-min(bReIndex[:6,6]) > max(bReIndex[6,:6])-min(bReIndex[6,:6]):
#     pOpt, pCov = curve_fit(normFunc1, fixedVals, resp.squeeze(), bounds=((0,0,0,0,0,0,0,0),
#     (np.inf,np.inf,360,360,1,1,0.10,b[6,6]+1)))
#     print(unit,pOpt)
# else:
#     pOpt, pCov = curve_fit(normFunc0, fixedVals, resp.squeeze(), bounds=((0,0,0,0,0,0,0,0),
#     (np.inf,np.inf,360,360,1,1,0.10,b[6,6]+1)))
#     print('using normFunc0')
#     print(unit,pOpt)
"""




