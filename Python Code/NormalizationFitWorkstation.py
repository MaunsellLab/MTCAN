    ## scipy curvefit Normalization parameters
    unitAlphaIndex = np.zeros(len(units))
    unitGaussFit = []
    unitPairedNormFit = []
    unitGaussR2 = np.zeros(len(units))
    unitPairedEV = np.zeros(len(units))
    unitPairedNormR2 = np.zeros(len(units))
    unitNormFitEstimate = [[] for i in range(len(units))]
    scaledLoc = np.zeros(len(units))
    scalar = np.zeros(len(units))
    for unitCount, unit in enumerate(units):
        b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000/trueStimDurMS
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
        tempMain = stimMat[:6, :6][:, reIndex]
        tempMain = np.squeeze(tempMain[:6, :6][reIndex, :])
        temp0Blank = stimMat[:6, 6][reIndex]
        temp1Blank = stimMat[6, :6][reIndex]
        stimMatReIndex[:6, :6] = tempMain
        stimMatReIndex[:6, 6] = temp0Blank
        stimMatReIndex[6, :6] = temp1Blank
        stimMatReIndex[6, 6] = stimMat[6, 6]

        # step 1: fitting gaussian function to single presentations
        # to extract gaussian fit and scaling factor (scalar on loc0)
        # if scalar >1 then loc0 has stronger resp
        resp = np.concatenate((bReIndex[6, :6], bReIndex[:6, 6]), axis=0)
        singleStim = np.concatenate((stimMatReIndex[6, :6], stimMatReIndex[:6, 6]),
                                    axis=0)
        fixedVals = fixedValsForCurveFit(prefDir, singleStim, stimIndexDict)
        if max(bReIndex[6, :6]) > max(bReIndex[:6, 6]):
            pOptGauss, pCov = curve_fit(gaussNormFunc0, fixedVals, resp.squeeze(),
                                        bounds=((0, 0, 0, 0, 0),
                                        (np.inf, max(resp)*1.5, 360, 360, 1)))
            y_pred = gaussNormFunc0(fixedVals, *pOptGauss)
            print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')
            sLoc = 1
        else:
            pOptGauss, pCov = curve_fit(gaussNormFunc1, fixedVals, resp.squeeze(),
                                        bounds=((0, 0, 0, 0, 0),
                                        (np.inf, max(resp)*1.5, 360, 360, 1)))
            y_pred = gaussNormFunc1(fixedVals, *pOptGauss)
            print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')
            print('scaled loc=0')
            sLoc = 0
        unitGaussR2[unitCount] = r2_score(resp.squeeze(), y_pred)
        unitGaussFit.append(pOptGauss)
        scalar[unitCount] = pOptGauss[4]

        # step 2: fitting  EMS normalization eqn with fixed gauss tuning curve
        # to paired stimuli
        # only two free parameters (a0, a1)
        resp = bReIndex[:6, :6].reshape(36)
        pairedStim = stimMatReIndex[:6, :6].reshape(36)
        pairedStim = pairedStim.astype(int)
        fixedVals = fixedValsCurveFitForPairedStim(prefDir, pairedStim, stimIndexDict, pOptGauss)
        if sLoc == 1:
            pOpt, pCov = curve_fit(EMSGaussNormStepFunc0, fixedVals, resp.squeeze(),
                                   bounds=((-1, -1), (5, 5)))
            y_pred = EMSGaussNormStepFunc0(fixedVals, *pOpt)
            print(r2_score(resp.squeeze(), y_pred))
            r2 = r2_score(resp.squeeze(), y_pred)
            EV = explained_variance_score(resp.squeeze(), y_pred)
            print(f'explained variance {EV}')
            print(pOpt)
            print('using EMS Func0, sLoc==1')
        else:
            pOpt, pCov = curve_fit(EMSGaussNormStepFunc1, fixedVals, resp.squeeze(),
                                   bounds=((-1, -1), (5, 5)))
            y_pred = EMSGaussNormStepFunc1(fixedVals, *pOpt)
            print(r2_score(resp.squeeze(), y_pred))
            r2 = r2_score(resp.squeeze(), y_pred)
            EV = explained_variance_score(resp.squeeze(), y_pred)
            print(f'explained variance {EV}')
            print(pOpt)

        # # step 3: fitting EMS normalization eqn with fixed gauss tuning curve
        # # but to only the non-preferred stimuli (null dir + rest in both locations)
        # resp = bReIndex[:6, 0]
        # pairedStim = stimMatReIndex[:6, 0]
        # pairedStim = pairedStim.astype(int)
        # fixedVals = fixedValsCurveFitForPairedStim(prefDir, pairedStim, stimIndexDict, pOptGauss)
        # if sLoc == 1:
        #     pOpt, pCov = curve_fit(EMSGaussNormStepFunc0, fixedVals, resp.squeeze(),
        #                            bounds=((-1, -1), (5, 5)))
        #     y_pred = EMSGaussNormStepFunc0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)
        #     print('using EMS Func0, sLoc==1')
        # else:
        #     pOpt, pCov = curve_fit(EMSGaussNormStepFunc1, fixedVals, resp.squeeze(),
        #                            bounds=((-1, -1), (5, 5)))
        #     y_pred = EMSGaussNormStepFunc1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)
        # loc0NonPrefAlpha = pOpt[0]
        #
        # #loc 1 non-pref
        # resp = bReIndex[0, :6]
        # pairedStim = stimMatReIndex[0, :6]
        # pairedStim = pairedStim.astype(int)
        # fixedVals = fixedValsCurveFitForPairedStim(prefDir, pairedStim, stimIndexDict, pOptGauss)
        # if sLoc == 1:
        #     pOpt, pCov = curve_fit(EMSGaussNormStepFunc0, fixedVals, resp.squeeze(),
        #                            bounds=((-1, -1), (5, 5)))
        #     y_pred = EMSGaussNormStepFunc0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)
        #     print('using EMS Func0, sLoc==1')
        # else:
        #     pOpt, pCov = curve_fit(EMSGaussNormStepFunc1, fixedVals, resp.squeeze(),
        #                            bounds=((-1, -1), (5, 5)))
        #     y_pred = EMSGaussNormStepFunc1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)
        # loc1NonPrefAlpha = pOpt[1]


        # # step 3, fit normalization eqn with fixed guass tuning curve to
        # # p+n, n+p conditions only, only two free parameters (a0,a1)
        # condensedMat = np.zeros((2, 2))
        # condensedMat[0, :] = bReIndex[0, [0, 3]]
        # condensedMat[1, :] = bReIndex[3, [0, 3]]
        # resp = condensedMat.reshape(4)
        #
        # condensedStimMat = np.zeros((2, 2))
        # condensedStimMat[0, :] = stimMatReIndex[0, [0, 3]]
        # condensedStimMat[1, :] = stimMatReIndex[3, [0, 3]]
        # fixedVals = fixedValsCurveFitForPairedStim(prefDir, condensedStimMat.reshape(4),
        #                                            stimIndexDict, pOpt)
        #
        # if sLoc == 1:
        #     pOpt, pCov = curve_fit(EMSGaussNormStepFunc0, fixedVals, resp.squeeze())
        #     print(pOpt)
        #     y_pred = EMSGaussNormStepFunc0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        # else:
        #     print('using norm func 0')
        #     pOpt, pCov = curve_fit(EMSGaussNormStepFunc1, fixedVals, resp.squeeze())
        #     print(pOpt)
        #     y_pred = EMSGaussNormStepFunc1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')

        # # Fitting two step EMS Normalization Equation with L0-L6
        # # instead of guassian
        # # step 1, fitting single stimulus with L0-L6 for direction tuning
        # resp = np.concatenate((bReIndex[6, :6], bReIndex[:6, 6]), axis=0)
        # singleStim = np.concatenate((stimMatReIndex[6, :6], stimMatReIndex[:6, 6]),
        #                             axis=0)
        # fixedVals = fixedValsForEMSGen(singleStim, stimIndexDict)
        # if max(bReIndex[6, :6]) > max(bReIndex[:6, 6]):
        #     pOpt, pCov = curve_fit(emsSingleStim0, fixedVals, resp.squeeze(),
        #                            bounds=((0, 0, 0, 0, 0, 0, 0),
        #                                    (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
        #     y_pred = emsSingleStim0(fixedVals, *pOpt)
        #     print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')
        #     sLoc = 1
        # else:
        #     pOpt, pCov = curve_fit(emsSingleStim1, fixedVals, resp.squeeze(),
        #                            bounds=((0, 0, 0, 0, 0, 0, 0),
        #                                    (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
        #     y_pred = emsSingleStim1(fixedVals, *pOpt)
        #     print(f'first step fit R2 {r2_score(resp.squeeze(), y_pred)}')
        #     print('scaled loc=0')
        #     sLoc = 0
        #
        # # step 2, fitting paired stimulus to EMS normalization eqn with
        # # fixed L0-L6 and S values. Only two free parameters (a0,a1)
        # resp = bReIndex[:6, :6].reshape(36)
        # pairedStim = stimMatReIndex[:6, :6].reshape(36)
        # pairedStim = pairedStim.astype(int)
        # fixedVals = fixedValsForPairedStimL0L6(pairedStim, stimIndexDict, pOpt)
        # if sLoc == 1:
        #     pOpt, pCov = curve_fit(EMSGenNormPaired0, fixedVals, resp.squeeze())
        #     y_pred = EMSGenNormPaired0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)
        #     print('using EMS Func0, sLoc==1')
        # else:
        #     pOpt, pCov = curve_fit(EMSGenNormPaired1, fixedVals, resp.squeeze())
        #     y_pred = EMSGenNormPaired1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     EV = explained_variance_score(resp.squeeze(), y_pred)
        #     print(f'explained variance {EV}')
        #     print(pOpt)

        unitPairedEV[unitCount] = EV
        unitPairedNormR2[unitCount] = r2
        scaledLoc[unitCount] = sLoc
        unitPairedNormFit.append(pOpt)
        unitAlphaIndex[unitCount] = (pOpt[0] + pOpt[1]) / 2


        # for fitting across mean spike count
        # resp = bReIndex.reshape(49)[:-1]
        # fixedVals = fixedValsForEMSGen(stimMatReIndex.reshape(49)[:-1], stimIndexDict)
        #
        # if max(bReIndex[:6, 6])-min(bReIndex[:6, 6]) > max(bReIndex[6, :6])-min(bReIndex[6, :6]):
        #     pOpt, pCov = curve_fit(emsNormGen1, fixedVals, resp, bounds=(
        #         (0, 0, 0, 0, 0, 0, 0, 0, 0),
        #         (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
        #     print(unit, pOpt)
        #     y_pred = emsNormGen1(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     sLoc = 0
        # else:
        #     print('using norm func 0')
        #     pOpt, pCov = curve_fit(emsNormGen0, fixedVals, resp, bounds=(
        #         (0, 0, 0, 0, 0, 0, 0, 0, 0),
        #         (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
        #     print(unit, pOpt)
        #     y_pred = emsNormGen0(fixedVals, *pOpt)
        #     print(r2_score(resp.squeeze(), y_pred))
        #     r2 = r2_score(resp.squeeze(), y_pred)
        #     sLoc = 1

def genericNormNoScalar(fixed, L0_0, L0_60, L0_120, L0_180, L0_240, L0_300,
                        L1_0, L1_60, L1_120, L1_180, L1_240, L1_300,
                        al0, al1, sig):
    """
    this function applies the generic normalization equation without a scalar
    at the weaker location. Instead, it gives that location its own L0-L6 value
    """
    c0, c1, l0, l1 = fixed
    L0 = np.array([L0_0, L0_60, L0_120, L0_180, L0_240, L0_300])
    L1 = np.array([L1_0, L1_60, L1_120, L1_180, L1_240, L1_300])
    # num = ((c0 * L0 * l0).sum(-1) + (c1 * L1 * l1).sum(-1))
    # denom = ((al0 * c0[:, 0]) + (al1 * c1[:, 0]) + sig)
    #
    # return num/denom
    loc0 = ((c0 * L0 * l0).sum(-1)) / (c0[:, 0] + (al1 * c1[:, 0]) + sig)
    loc1 = ((c1 * L1 * l1).sum(-1)) / ((al0 * c0[:, 0]) + c1[:, 0] + sig)

    return loc0 + loc1

pOpt, pCov = curve_fit(genericNormNoScalar, fixedVals, resp.squeeze(),
                       bounds=((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                               (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
                                np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
                                np.inf, np.inf, .10)))