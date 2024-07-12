"""
Functions for trying different normalization models
"""

import numpy as np

def genericNorm(fixed, L0_0, L0_1, L1_0, L1_1, sig):
    """
    generic normalization
    """

    c0, c1, l0_0, l0_1, l1_0, l1_1 = fixed
    L0 = np.tile(np.array([L0_0, L0_1]), (c0.shape[0], 1))
    L1 = np.tile(np.array([L1_0, L1_1]), (c1.shape[0], 1))

    # generic norm_num
    c0_num = c0[:,None] * L0 * np.array([l0_0, l0_1]).T
    c1_num = c1[:,None] * L1 * np.array([l1_0, l1_1]).T
    num = c0_num.sum(-1) + c1_num.sum(-1)
    denom = c0 + c1 + sig

    return num / denom


def weightedNorm(fixed, L0_0, L0_1, L1_0, L1_1, sig):
    """
    RF weighted normalization
    """

    c0, c1, l0_0, l0_1, l1_0, l1_1 = fixed
    L0 = np.tile(np.array([L0_0, L0_1]), (c0.shape[0], 1))
    L1 = np.tile(np.array([L1_0, L1_1]), (c1.shape[0], 1))

    # generic norm_num
    c0_num = c0[:,None] * L0 * np.array([l0_0, l0_1]).T
    c1_num = c1[:,None] * L1 * np.array([l1_0, l1_1]).T * L1_0/L0_0
    num = c0_num.sum(-1) + c1_num.sum(-1)
    denom = c0 + (c1*L1_0/L0_0) + sig

    return num / denom


def heegerNorm(fixed, Lp, Lnp, sig):
    """
    generic normalization Heeger model
    """

    c0, c1, l0_0, l0_1, l1_0, l1_1 = fixed
    L0 = np.tile(np.array([Lp, Lnp]), (c0.shape[0], 1))
    L1 = np.tile(np.array([Lp, Lnp]), (c1.shape[0], 1))

    # generic norm_num
    c0_num = c0[:, None] * L0 * np.array([l0_0, l0_1]).T
    c1_num = c1[:, None] * L1 * np.array([l1_0, l1_1]).T
    num = c0_num.sum(-1) + c1_num.sum(-1)
    denom = c0 + c1 + sig

    return num / denom


def fixedValsForMTND(pResp, nResp, pnResp, npResp, ppResp, nnResp, sponResp, offPos):
    """
    this function will return the fixed values and the responses so that they can be
    modelled using curve_fit. An array of responses are made and so are an array of the
    fixed values (such as contrast and which stimulus was where during each condition)
    fixedVals = [C0, C1, L0, L1]
    """

    responses = [pResp[0], pResp[offPos], nResp[0], nResp[offPos],
                 pnResp[offPos-1], npResp[offPos-1], ppResp[offPos-1],
                 nnResp[offPos-1], sponResp]

    fixedVals = [[1, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0],
                 [1, 0, 0, 1, 0, 0], [0, 1, 0, 0, 0, 1],
                 [1, 1, 1, 0, 0, 1], [1, 1, 0, 1, 1, 0],
                 [1, 1, 1, 0, 1, 0], [1, 1, 0, 1, 0, 1],
                 [0, 0, 0, 0, 0, 0]]

    return np.array(responses), np.array(fixedVals)






