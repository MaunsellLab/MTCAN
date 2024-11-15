"""
Functions for trying different normalization models
"""
# imports
import numpy as np

# for MTND plug-in
def genericNorm(fixed, L0_0, L0_1, L1_0, L1_1, sig):
    """
    generic normalization
    """

    c0, c1, l0_0, l0_1, l1_0, l1_1 = fixed
    L0 = np.tile(np.array([L0_0, L0_1]), (c0.shape[0], 1))
    L1 = np.tile(np.array([L1_0, L1_1]), (c1.shape[0], 1))

    # generic norm_num
    c0_num = c0[:, None] * L0 * np.array([l0_0, l0_1]).T
    c1_num = c1[:, None] * L1 * np.array([l1_0, l1_1]).T
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
    c0_num = c0[:, None] * L0 * np.array([l0_0, l0_1]).T
    c1_num = c1[:, None] * L1 * np.array([l1_0, l1_1]).T  # * L1_0/L0_0
    num = c0_num.sum(-1) + c1_num.sum(-1)
    denom = c0 + (c1*L1_0/L0_0) + sig

    return num / denom


def weightedNormFewerParams(fixed, Lp, Lnp, sig):
    """
    RF weighted normalization but fewer params such that the RF weight is fixed and
    sent in through the fixed variables
    """

    c0, c1, l0_0, l0_1, l1_0, l1_1, x1 = fixed
    L0 = np.tile(np.array([Lp, Lnp]), (c0.shape[0], 1))
    L1 = np.tile(np.array([Lp, Lnp]), (c1.shape[0], 1))

    c0_num = c0[:, None] * L0 * np.array([l0_0, l0_1]).T
    c1_num = c1[:, None] * L1 * np.array([l1_0, l1_1]).T * x1[:, None]
    num = c0_num.sum(-1) + c1_num.sum(-1)
    denom = c0 + (c1*x1) + sig

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


def emsNorm(fixed, Lp, Lnp, a0, a1, sig):
    """
    EMS normalization Heeger model
    """

    c0, c1, l0_0, l0_1, l1_0, l1_1 = fixed
    L0 = np.tile(np.array([Lp, Lnp]), (c0.shape[0], 1))
    L1 = np.tile(np.array([Lp, Lnp]), (c1.shape[0], 1))

    loc0Num = (c0[:, None] * L0 * np.array([l0_0, l0_1]).T).sum(-1)
    loc1Denom = c0 + (c1*a1) + sig
    loc0 = loc0Num/loc1Denom

    loc1Num = (c1[:, None] * L1 * np.array([l1_0, l1_1]).T).sum(-1)
    loc1Denom = (c0*a0) + c1 + sig
    loc1 = loc1Num + loc1Denom

    return loc0 + loc1


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


def fixedValsForRFWeightFewerParams(pResp, nResp, pnResp, npResp, ppResp, nnResp, sponResp, offPos):
    """
    this function will return the fixed values and the responses so that they can be
    modelled using curve_fit. An array of responses are made and so are an array of the
    fixed values (such as contrast and which stimulus was where during each condition)
    fixedVals = [C0, C1, L0, L1]
    """

    responses = [pResp[0], pResp[offPos], nResp[0], nResp[offPos],
                 pnResp[offPos-1], npResp[offPos-1], ppResp[offPos-1],
                 nnResp[offPos-1], sponResp]

    x1 = pResp[offPos] / pResp[0]

    fixedVals = [[1, 0, 1, 0, 0, 0, x1], [0, 1, 0, 0, 1, 0, x1],
                 [1, 0, 0, 1, 0, 0, x1], [0, 1, 0, 0, 0, 1, x1],
                 [1, 1, 1, 0, 0, 1, x1], [1, 1, 0, 1, 1, 0, x1],
                 [1, 1, 1, 0, 1, 0, x1], [1, 1, 0, 1, 0, 1, x1],
                 [0, 0, 0, 0, 0, 0, x1]]

    return np.array(responses), np.array(fixedVals)


# Updated response model for element-wise operation
def fullRFWeightedNorm(contrast_center, contrast_periphery, location, stim_type_center,
                       stim_type_periphery, Lp, Lnp, W1, W2, W3, W4, sigma):
    """
    RF weighted norm fit for full stimulus set, fitting all 26 conditions at once (excluding baseline and mapping)
    """

    # Determine weight based on peripheral location
    if location == -1:  # No peripheral stimulus
        W_periphery = 0
    else:
        W_periphery = [W1, W2, W3, W4][int(location)]  # Choose weight based on location

    # Center stimulus response
    if stim_type_center == 1:  # Preferred
        L_center = Lp
    elif stim_type_center == 0:  # Non-preferred
        L_center = Lnp
    else:  # No center stimulus
        L_center = 0

    # Peripheral stimulus response
    if stim_type_periphery == 1:  # Preferred
        L_periphery = Lp
    elif stim_type_periphery == 0:  # Non-preferred
        L_periphery = Lnp
    else:  # No periphery stimulus
        L_periphery = 0

    # Equation for the response R
    numerator = (contrast_center * L_center) + (contrast_periphery * L_periphery * W_periphery)
    denominator = contrast_center + (W_periphery * contrast_periphery) + sigma

    return numerator / denominator


# Wrapper function to pass all data element-wise to curve_fit
def model_wrapper(data, Lp, Lnp, W1, W2, W3, W4, sigma):
    contrast_center, contrast_periphery, location, stim_type_center, stim_type_periphery = data
    # Apply the response model element-wise
    return [fullRFWeightedNorm(c_center, c_periph, loc, stim_c, stim_p, Lp, Lnp, W1, W2, W3, W4, sigma)
            for c_center, c_periph, loc, stim_c, stim_p in
            zip(contrast_center, contrast_periphery, location, stim_type_center, stim_type_periphery)]


# Define the function to calculate predicted responses using the fitted parameters
def apply_fitted_model(contrast_center, contrast_periphery, location, stim_type_center,
                       stim_type_periphery, Lp, Lnp, W1, W2, W3, W4, sigma):
    # Same logic as in response_model, but handles arrays
    predicted_responses = []
    for c_center, c_periph, loc, stim_c, stim_p in zip(contrast_center, contrast_periphery, location,
                                                       stim_type_center, stim_type_periphery):
        # Call the response model with the fitted parameters for each data point
        pred = fullRFWeightedNorm(c_center, c_periph, loc, stim_c, stim_p, Lp, Lnp, W1, W2, W3, W4, sigma)
        predicted_responses.append(pred)
    return predicted_responses

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# for MTNC plug-in

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def fixedValsForRFWeightMTNC(stimMatReIndex, stimIndexDict, rfweight):
    """
    function will create the fixed values for running the
    RF weighted normalization curve_fit, where the weight is a fixed value
    """

    c0s, c1s, l0s, l1s = [], [], [], []
    direction_set = np.arange(0, 360, 60)
    if type(stimMatReIndex) == int:
        c0 = stimIndexDict[stimMatReIndex][0]['contrast']
        l0 = stimIndexDict[stimMatReIndex][0]['direction']
        c1 = stimIndexDict[stimMatReIndex][1]['contrast']
        l1 = stimIndexDict[stimMatReIndex][1]['direction']

        # Make one-hot encoding of l0 and l1
        l0_oh = np.zeros(6)
        l1_oh = np.zeros(6)
        l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
        l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
        c0s.append(np.repeat(c0, 6))
        c1s.append(np.repeat(c1, 6))
        l0s.append(l0_oh)
        l1s.append(l1_oh)
    else:
        for i in stimMatReIndex:
            c0 = stimIndexDict[i][0]['contrast']
            l0 = stimIndexDict[i][0]['direction']
            c1 = stimIndexDict[i][1]['contrast']
            l1 = stimIndexDict[i][1]['direction']

            # Make one-hot encoding of l0 and l1
            l0_oh = np.zeros(6)
            l1_oh = np.zeros(6)
            l0_oh[np.argwhere(direction_set == l0).squeeze()] = 1
            l1_oh[np.argwhere(direction_set == l1).squeeze()] = 1
            c0s.append(np.repeat(c0, 6))
            c1s.append(np.repeat(c1, 6))
            l0s.append(l0_oh)
            l1s.append(l1_oh)

    rfw = np.repeat(rfweight, 294).reshape(49, 6)

    return np.array(c0s), np.array(c1s), np.array(l0s), np.array(l1s), np.array(rfw)


def rfWeightCondensed(fixed, L0_0, L0_1, L1_0, L1_1, sig):
    """
    generic normalization condensed
    """

    c0, c1, l0, l1, weight = fixed
    L0 = np.array([L0_0, L0_1])
    L1 = np.array([L1_0, L1_1])

    # generic norm
    num = ((c0 * (L0 ** 2) * l0).sum(-1) + (c1 * (L1 ** 2) * l1).sum(-1))
    denom = (c0[:, 0] + (c1[:, 0] *(L1_0 / L0_0)) + sig)

    return num / denom


def rfWeightMTNC(fixed, L0_0, L0_60, L0_120, L0_180, L0_240, L0_300,
                 L1_0, L1_60, L1_120, L1_180, L1_240, L1_300, sig, base):
    """
    this function applies the RF weighted normalization equation.
    It gives each location its own L0-L6 value
    """

    c0, c1, l0, l1, weight = fixed
    L0 = np.array([L0_0, L0_60, L0_120, L0_180, L0_240, L0_300])
    L1 = np.array([L1_0, L1_60, L1_120, L1_180, L1_240, L1_300])

    # generic norm
    num = ((c0 * (L0 ** 2) * l0).sum(-1) + (c1 * (L1 ** 2) * l1).sum(-1))
    denom = c0[:, 0] + (c1[:, 0] * weight[:, 0]) + sig

    return (num/denom) + base


def rfWeightWithAlphaMTNC(fixed, L0_0, L0_60, L0_120, L0_180, L0_240, L0_300,
                          L1_0, L1_60, L1_120, L1_180, L1_240, L1_300, al1, sig, base):
    """
    this function applies the RF weighted normalization equation.
    It gives each location its own L0-L6 value
    """

    c0, c1, l0, l1, weight = fixed
    L0 = np.array([L0_0, L0_60, L0_120, L0_180, L0_240, L0_300])
    L1 = np.array([L1_0, L1_60, L1_120, L1_180, L1_240, L1_300])

    # generic norm
    num = ((c0 * L0 * l0).sum(-1) + (c1 * L1 * l1).sum(-1))
    denom = c0[:, 0] + (c1[:, 0] * weight[:, 0] * al1) + sig

    return (num/denom) + base


def EMSNormMTNC(fixed, L0, L60, L120, L180, L240, L300, al0, al1, sig, base):
    """
    this function applies the RF weighted normalization equation.
    It gives each location its own L0-L6 value
    """

    c0, c1, l0, l1 = fixed
    L = np.array([L0, L60, L120, L180, L240, L300])

    # ems
    loc0 = ((c0 * L * l0).sum(-1)) / (c0[:, 0] + (al1 * c1[:, 0]) + sig)

    loc1 = ((c1 * L * l1).sum(-1)) / ((al0 * c0[:, 0]) + c1[:, 0] + sig)

    return loc0 + loc1 + base


def bramNormMTNC(fixed, L0_0, L0_60, L0_120, L0_180, L0_240, L0_300,
                 L1_0, L1_60, L1_120, L1_180, L1_240, L1_300, al0,
                 al1, sig, base):
    """
    this function applies the generic normalization equation without a scalar
    at the weaker location. Instead, it gives that location its own L0-L6 value
    """

    c0, c1, l0, l1 = fixed
    L0 = np.array([L0_0, L0_60, L0_120, L0_180, L0_240, L0_300])
    L1 = np.array([L1_0, L1_60, L1_120, L1_180, L1_240, L1_300])

    # generic norm
    num = ((c0 * (L0 * l0)).sum(-1) + (c1 * (L1 * l1)).sum(-1))
    denom = ((al0 * c0[:, 0]) + (al1 * c1[:, 0]) + sig)

    return num/denom + base


def rfWeightedNormCondensed(fixed, L0_0, L0_1, L1_0, L1_1, w, sig):
    """
    RF weighted normalization for smaller stimulus set (condensed)
    """

    c0, c1, l0, l1 = fixed
    loc0 = np.array([L0_0, L0_1])
    loc1 = np.array([L1_0, L1_1])


    # generic norm
    num = ((c0 * (loc0 ** 2) * l0).sum(-1) + (c1 * (loc1 ** 2) * l1).sum(-1))
    denom = ((1 * c0[:, 0]) + ((w ** 2) * c1[:, 0])) + (sig ** 2)

    return num / denom


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# for MTSIGMA plug-in

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def normMTSIG(c, L, sigma, n):
    """
    equation to fit single gabor presentation of a stimulus at different contrasts to get
    the L term and the sigma term
    """

    return (c**n * L**n) / (c**n + sigma**n)


def weightedNormMTSIG(C, L, w, sigma, condition):
    """
    weighted normalization to fit single gabor presentation of a stimulus at different contrasts
    at two different locations (center vs periphery), will return the L term, weight (w) and sigma
    """

    if condition:
        return (C * L * w) / ((C * w) + sigma)
    else:
        return (C * L) / (C + sigma)




