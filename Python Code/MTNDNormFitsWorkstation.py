import matplotlib.pyplot as plt
import numpy as np


def fullClassicWeightedNorm(contrast_center, contrast_periphery, location, stim_type_center,
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
    denominator = contrast_center + contrast_periphery + sigma

    return numerator / denominator


# Wrapper function to pass all data element-wise to curve_fit
def model_wrapperClassicWeightedNorm(data, Lp, Lnp, W1, W2, W3, W4, sigma):
    contrast_center, contrast_periphery, location, stim_type_center, stim_type_periphery = data
    # Apply the response model element-wise
    return [fullClassicWeightedNorm(c_center, c_periph, loc, stim_c, stim_p, Lp, Lnp, W1, W2, W3, W4, sigma)
            for c_center, c_periph, loc, stim_c, stim_p in
            zip(contrast_center, contrast_periphery, location, stim_type_center, stim_type_periphery)]


# Define the function to calculate predicted responses using the fitted parameters
def apply_fitted_modelClassicWeighted(contrast_center, contrast_periphery, location, stim_type_center,
                                      stim_type_periphery, Lp, Lnp, W1, W2, W3, W4, sigma):
    # Same logic as in response_model, but handles arrays
    predicted_responses = []
    for c_center, c_periph, loc, stim_c, stim_p in zip(contrast_center, contrast_periphery, location,
                                                       stim_type_center, stim_type_periphery):
        # Call the response model with the fitted parameters for each data point
        pred = fullClassicWeightedNorm(c_center, c_periph, loc, stim_c, stim_p, Lp, Lnp, W1, W2, W3, W4, sigma)
        predicted_responses.append(pred)
    return predicted_responses

# known variables
contrast_center = [1, 0, 0, 0, 0,
                   1, 0, 0, 0, 0,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1]  # First response has contrast in the center, second in periphery
contrast_periphery = [0, 1, 1, 1, 1,
                      0, 1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1]  # First has no periphery, second has periphery with contrast
locations = np.array([-1, 0, 1, 2, 3,
                      -1, 0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3])  # First has no peripheral stimulus, second is at location 1
stim_type_center = np.array([1, -1, -1, -1, -1,
                             0, -1, -1, -1, -1,
                             1, 1, 1, 1,
                             0, 0, 0, 0,
                             1, 1, 1, 1,
                             0, 0, 0, 0])
stim_type_periphery = np.array([-1, 1, 1, 1, 1,
                                -1, 0, 0, 0, 0,
                                0, 0, 0, 0,
                                1, 1, 1, 1,
                                1, 1, 1, 1,
                                0, 0, 0, 0])

# individual neurons fit
fullClassicNormPred = []
measuredResp = []
r2ScoresClassicNormFull = []
for i in range(len(prefNormalized)):
    resp = np.concatenate((prefNormalized[i], nonprefNormalized[i], pnNormalized[i],
                           npNormalized[i], ppNormalized[i], nnNormalized[i]), axis=0)

    initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10]  # Lp, Lnp, W1, W2, W3, W4, sigma

    popt, pcov = curve_fit(model_wrapperClassicWeightedNorm, (contrast_center, contrast_periphery, locations, stim_type_center,
                                                              stim_type_periphery),
                           resp, p0=initial_guess, maxfev=10000000)
    y_pred1 = apply_fitted_modelClassicWeighted(contrast_center, contrast_periphery, locations,
                                                stim_type_center, stim_type_periphery, *popt)
    r2 = r2_score(resp, y_pred1)
    r2ScoresClassicNormFull.append(r2)

    # single presentation prediction
    pCenterSinglePred = np.array(fullClassicWeightedNorm(1, 0, -1, 1, -1, *popt))
    pCenterSinglePred = pCenterSinglePred[np.newaxis]
    pPeriSinglePred = np.array([fullClassicWeightedNorm(0, 1, j, -1, 1, *popt) for j in range(4)])
    pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

    npCenterSinglePred = np.array(fullClassicWeightedNorm(1, 0, -1, 0, -1, *popt))
    npCenterSinglePred = npCenterSinglePred[np.newaxis]
    npPeriSinglePred = np.array([fullClassicWeightedNorm(0, 1, j, -1, 0, *popt) for j in range(4)])
    npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

    # paired presentation prediction
    pnFitPred = [fullClassicWeightedNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
    npFitPred = [fullClassicWeightedNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
    ppFitPred = [fullClassicWeightedNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
    nnFitPred = [fullClassicWeightedNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

    predResp = np.concatenate((pSingleFitPred, npSingleFitPred, pnFitPred,
                               npFitPred, ppFitPred, nnFitPred), axis=0)

    for j in range(len(resp)):
        fullClassicNormPred.append(predResp[j])
        measuredResp.append(resp[j])

    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')
    popt_str = (f"Lp = {popt[0]:.2f}\n"
                f"Lnp = {popt[1]:.2f}\n"
                f"W1 = {popt[2]:.2f}\n"
                f"W2 = {popt[3]:.2f}\n"
                f"W3 = {popt[4]:.2f}\n"
                f"W4 = {popt[5]:.2f}\n"
                f"sigma = {popt[6]:.2f}\n"
                f"r2={r2:.3f}")
    mSize = 70
    x = np.arange(0, numSteps + 1)

    ax1.plot(x, prefNormalized[i], color='black', label='PO', linewidth=1.5)
    ax1.scatter(x, prefNormalized[i], s=mSize,
                facecolor='white', edgecolor='black')
    ax1.plot(x, nonprefNormalized[i], color='grey', label='NO', linewidth=1.5)
    ax1.scatter(x, nonprefNormalized[i], color='grey', s=mSize,
                facecolor='white', edgecolor='gray')
    ax1.plot(xNorm, pnNormalized[i], color='green', label='P0 N1', linewidth=1.5)
    ax1.scatter(xNorm, pnNormalized[i], color='green', s=mSize)
    ax1.plot(xNorm, npNormalized[i], color='magenta', label='N0 P1', linewidth=1.5)
    ax1.scatter(xNorm, npNormalized[i], color='magenta', s=mSize)
    ax1.plot(xNorm, ppNormalized[i], color='black', label='P0 P1',
             linewidth=1.5)
    ax1.scatter(xNorm, ppNormalized[i], color='black', s=mSize)
    ax1.plot(xNorm, nnNormalized[i], color='grey', label='N0 N1',
             linewidth=1.5)
    ax1.scatter(xNorm, nnNormalized[i], color='grey', s=mSize)
    ax1.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=mSize)
    ax1.plot(xNorm, npFitPred, color='magenta', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, npFitPred, facecolors='magenta', edgecolors='m', s=mSize)
    ax1.plot(xNorm, ppFitPred, color='black', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=mSize)
    ax1.plot(xNorm, nnFitPred, color='grey', linestyle='--', linewidth=1.5)
    ax1.scatter(xNorm, nnFitPred, facecolors='gray', edgecolors='gray', s=mSize)
    ax1.plot(x, pSingleFitPred, color='black', linestyle='--', linewidth=1.5)
    ax1.scatter(x, pSingleFitPred, facecolors='white', edgecolors='black', s=mSize)
    ax1.plot(x, npSingleFitPred, color='grey', linestyle='--', linewidth=1.5)
    ax1.scatter(x, npSingleFitPred, facecolors='white', edgecolors='gray', s=mSize)
    ax1.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    ax1.set_ylabel('Response Rate (spikes/s)', fontsize=25, **hfont)
    ax1.set_ylim([0, 1.6])
    ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
             transform=ax1.transAxes)
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        f'{round(normVal[i][0])}' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.set_xticks([0, 1, 2, 3, 4])
    ax1.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4.5])
    ax1.tick_params(axis='both', width=2, length=8)
    ax1.tick_params(axis='both', width=2, length=8)

    # Save the figure to the specified path
    filePath = f'/Users/chery/Desktop/Project 2 Figs/individualNeuronsFits/{masterGoodUnits[i]}_weightedHeeger.pdf'
    fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    plt.close('all')

# population fits
resp = np.concatenate((prefMean, nonprefMean, pnMean, npMean, ppMean, nnMean), axis=0)
initial_guess = [1.0, 0.5, 1.0, 0.75, 0.50, .10, 0.10]  # Lp, Lnp, W1, W2, W3, W4, sigma
popt, pcov = curve_fit(model_wrapperClassicWeightedNorm, (contrast_center, contrast_periphery,
                                                          locations, stim_type_center,
                                                          stim_type_periphery),
                       resp, p0=initial_guess)

y_pred1 = apply_fitted_modelClassicWeighted(contrast_center, contrast_periphery, locations,
                                            stim_type_center, stim_type_periphery, *popt)
r2 = r2_score(resp, y_pred1)

pnFitPred = [fullClassicWeightedNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [fullClassicWeightedNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullClassicWeightedNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullClassicWeightedNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullClassicWeightedNorm(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullClassicWeightedNorm(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullClassicWeightedNorm(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullClassicWeightedNorm(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# main heeger weighed plot
hfont = {'fontname': 'Arial'}
for filler in range(1):
    x = np.arange(0, numSteps + 1)
    fig, ax2 = plt.subplots(figsize=(7, 7), layout='constrained')

    ax2.plot(x, prefMean, color='black', label='PO', linewidth=1.5)
    ax2.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
                 markerfacecolor='none', markeredgecolor='black', markersize=7)
    ax2.plot(x, nonprefMean, color='grey', label='NO', linewidth=1.5)
    ax2.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
                 markerfacecolor='none', markeredgecolor='grey', markersize=7)
    ax2.plot(xNorm, pnMean, color='green', label='P0 N1', linewidth=1.5)
    ax2.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
                 color='green', markersize=7)
    ax2.plot(xNorm, npMean, color='red', label='N0 P1', linewidth=1.5)
    ax2.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='red',
                 color='red', markersize=7)
    ax2.plot(xNorm, ppMean, color='black', label='P0 P1', linewidth=1.5)
    ax2.errorbar(xNorm, ppMean, yerr=ppSEM, fmt='o', ecolor='black',
                 color='black', markersize=7)
    ax2.plot(xNorm, nnMean, color='grey', label='N0 N1', linewidth=1.5)
    ax2.errorbar(xNorm, nnMean, yerr=nnSEM, fmt='o', ecolor='grey',
                 color='grey', markersize=7)
    ax2.plot(xNorm, pnFitPred, color='green', linestyle='--', linewidth=1.5)
    ax2.scatter(xNorm, pnFitPred, facecolors='green', edgecolors='g', s=50)
    ax2.plot(xNorm, npFitPred, color='red', linestyle='--', linewidth=1.5)
    ax2.scatter(xNorm, npFitPred, facecolors='red', edgecolors='r', s=50)
    ax2.plot(xNorm, ppFitPred, color='black', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax2.scatter(xNorm, ppFitPred, facecolors='black', edgecolors='black', s=50)
    ax2.plot(xNorm, nnFitPred, color='grey', linestyle='--', alpha=0.8,
             linewidth=1.5)
    ax2.scatter(xNorm, nnFitPred, facecolors='gray', edgecolors='gray', s=50)
    ax2.plot(x, pSingleFitPred, color='black', linestyle='--', linewidth=1.5)
    ax2.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
    ax2.plot(x, npSingleFitPred, color='grey', linestyle='--', linewidth=1.5)
    ax2.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
    ax2.set_xlabel('Receptive Field Location', fontsize=25, **hfont)
    ax2.set_ylabel('Normalized Response Rate', fontsize=25, **hfont)
    ax2.set_ylim([0, 1.50])
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_ticks = ax2.get_yticks()
    y_tick_labels = [
        '0' if tick == 0 else
        '1' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax2.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax2.set_xticks([0, 1, 2, 3, 4])
    ax2.set_xticklabels(['0', '1', '2', '3', '4'], **hfont, fontsize=20)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['left'].set_linewidth(2)
    ax2.set_xlim([-0.5, 4.5])
    ax2.tick_params(axis='both', width=2, length=8)

    plt.show()
    # # Save the figure to the specified path
    # filePath = f'/Users/chery/Documents/Grad School/Maunsell Lab/Manuscripts/Python Figures/populationWeightedHeeger.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')


fig, ax2 = plt.subplots()
ax2.scatter(measuredResp, fullClassicNormPred, color='black', facecolors='none')
ax2.set_ylim([0, 2])
ax2.set_xlim([0, 2])
ax2.set_xlabel('Measured Response', fontsize=15, **hfont)
ax2.set_ylabel('Predicted Response', fontsize=15, **hfont)
line = lines.Line2D([0, 1], [0, 1], color='black')
transform = ax2.transAxes
line.set_transform(transform)
ax2.add_line(line)

plt.show()


########################################################################################################################

############################################################ Classic Heeger ############################################

########################################################################################################################


def fullClassicHeegerNorm(contrast_center, contrast_periphery, stim_type_center,
                          stim_type_periphery, Lp, Lnp, sigma):
    """
    classic heeger norm fit for full stimulus set, fitting all 26 conditions at once (excluding baseline and mapping)
    """

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
    numerator = (contrast_center * L_center) + (contrast_periphery * L_periphery)
    denominator = contrast_center + contrast_periphery + sigma

    return numerator / denominator


# Wrapper function to pass all data element-wise to curve_fit
def model_wrapperClassicHeegerNorm(data, Lp, Lnp, sigma):
    contrast_center, contrast_periphery, stim_type_center, stim_type_periphery = data
    # Apply the response model element-wise
    return [fullClassicHeegerNorm(c_center, c_periph, stim_c, stim_p, Lp, Lnp, sigma)
            for c_center, c_periph, stim_c, stim_p in
            zip(contrast_center, contrast_periphery, stim_type_center, stim_type_periphery)]


# Define the function to calculate predicted responses using the fitted parameters
def apply_fitted_modelClassicHeeger(contrast_center, contrast_periphery, stim_type_center,
                                    stim_type_periphery, Lp, Lnp, sigma):
    predicted_responses = []
    for c_center, c_periph, stim_c, stim_p in zip(contrast_center, contrast_periphery,
                                                  stim_type_center, stim_type_periphery):
        # Call the response model with the fitted parameters for each data point
        pred = fullClassicHeegerNorm(c_center, c_periph, stim_c, stim_p, Lp, Lnp, sigma)
        predicted_responses.append(pred)
    return predicted_responses

# known variables
contrast_center = [1, 0, 0, 0, 0,
                   1, 0, 0, 0, 0,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1]  # First response has contrast in the center, second in periphery
contrast_periphery = [0, 1, 1, 1, 1,
                      0, 1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1]  # First has no periphery, second has periphery with contrast
stim_type_center = np.array([1, -1, -1, -1, -1,
                             0, -1, -1, -1, -1,
                             1, 1, 1, 1,
                             0, 0, 0, 0,
                             1, 1, 1, 1,
                             0, 0, 0, 0])
stim_type_periphery = np.array([-1, 1, 1, 1, 1,
                                -1, 0, 0, 0, 0,
                                0, 0, 0, 0,
                                1, 1, 1, 1,
                                1, 1, 1, 1,
                                0, 0, 0, 0])

resp = np.concatenate((prefMean, nonprefMean, pnMean, npMean, ppMean, nnMean), axis=0)
initial_guess = [1.0, 0.5, 0.10]  # Lp, Lnp, sigma
popt, pcov = curve_fit(model_wrapperClassicHeegerNorm, (contrast_center, contrast_periphery,
                                                        stim_type_center, stim_type_periphery),
                       resp, p0=initial_guess)

y_pred1 = apply_fitted_modelClassicHeeger(contrast_center, contrast_periphery,
                                          stim_type_center, stim_type_periphery, *popt)
r2 = r2_score(resp, y_pred1)

pnFitPred = [fullClassicHeegerNorm(1, 1, 1, 0, *popt) for j in range(4)]
npFitPred = [fullClassicHeegerNorm(1, 1, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullClassicHeegerNorm(1, 1, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullClassicHeegerNorm(1, 1, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullClassicHeegerNorm(1, 0, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullClassicHeegerNorm(0, 1, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullClassicHeegerNorm(1, 0, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullClassicHeegerNorm(0, 1, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# plotting
fig, ax2 = plt.subplots()
x = np.arange(0, numSteps + 1)
ax2.plot(x, prefMean, color='black', label='PO')
ax2.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
             color='black', markersize=7)
ax2.plot(x, nonprefMean, color='grey', label='NO')
ax2.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
             color='grey', markersize=7)
ax2.plot(xNorm, pnMean, color='green', label='P0 N1')
ax2.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
             color='green', markersize=7)
ax2.plot(xNorm, npMean, color='red', label='N0 P1')
ax2.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='red',
             color='red', markersize=7)
ax2.plot(xNorm, ppMean, color='black', label='P0 P1')
ax2.errorbar(xNorm, ppMean, yerr=pnSEM, fmt='o', ecolor='black',
             color='black', markersize=7)
ax2.plot(xNorm, nnMean, color='grey', label='N0 N1')
ax2.errorbar(xNorm, nnMean, yerr=npSEM, fmt='o', ecolor='grey',
             color='grey', markersize=7)
ax2.plot(xNorm, pnFitPred, color='green', linestyle='--')
ax2.scatter(xNorm, pnFitPred, facecolors='none', edgecolors='g', s=50)
ax2.plot(xNorm, npFitPred, color='red', linestyle='--')
ax2.scatter(xNorm, npFitPred, facecolors='none', edgecolors='r', s=50)
ax2.plot(xNorm, ppFitPred, color='black', linestyle=':', alpha=0.8)
ax2.scatter(xNorm, ppFitPred, facecolors='none', edgecolors='black', s=50)
ax2.plot(xNorm, nnFitPred, color='grey', linestyle=':', alpha=0.8)
ax2.scatter(xNorm, nnFitPred, facecolors='none', edgecolors='gray', s=50)
ax2.plot(x, pSingleFitPred, color='black', linestyle='--')
ax2.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
ax2.plot(x, npSingleFitPred, color='grey', linestyle='--')
ax2.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
ax2.set_xlabel('Receptive Field Offset', fontsize=15)
ax2.set_ylabel('Normalized Response Rate', fontsize=15)
ax2.set_ylim([0, 1.50])
ax2.axhline(y=meanSpon, linestyle='--', color='blue', alpha=0.8)
ax2.fill_between(x=x, y1=meanSpon + semSpon, y2=meanSpon - semSpon,
                 color='blue', alpha=0.2)
ax2.set_title('Normalization vs Distance PN/NP')
plt.show()

# bar plots for PN
for filler in range(1):
    hfont = {'fontname': 'Arial'}
    prefCenterBar = prefMean[0]
    barResponsesPN = {'PN measured': [pnMean[0], pnMean[3]],
                      'PN predicted': [(prefMean[0] + nonprefMean[1]) / 2,
                                       (prefMean[0] + nonprefMean[4]) / 2],
                      'Non-pref offset': [nonprefMean[1], nonprefMean[3]]}

    fig, ax1 = plt.subplots(layout='constrained')
    fig.set_size_inches(10, 8)

    width = 0.20
    colorList = ['tab:green', 'tab:green', 'grey']
    cluster_positions = [0, 1, 2]
    x_single = [cluster_positions[0]]  # For Position 0
    x_multiple = [cluster_positions[1], cluster_positions[2]]  # For Positions 1 and 4

    ax1.bar(x_single, prefCenterBar, width, label="Pref Center", color="black")

    multiplier = 0
    for attribute, measurement in barResponsesPN.items():
        for i, x_pos in enumerate(x_multiple):
            offset = width * multiplier
            if multiplier == 1:
                rects = ax1.bar(x_pos + offset, measurement[i], width, label=attribute if i == 0 else "",
                                color=colorList[multiplier], hatch='//', alpha=0.5)
            else:
                rects = ax1.bar(x_pos + offset, measurement[i], width, label=attribute if i == 0 else "",
                                color=colorList[multiplier])
        multiplier += 1

    # Calculate the center of each cluster
    tick_positions = [cluster_positions[0],
                      cluster_positions[1] + width,
                      cluster_positions[2] + width]

    # Set custom x-axis tick labels at the center of each cluster
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(['Position 0', 'Position 1', 'Position 4'],
                        **hfont, fontsize=15)
    ax1.set_ylim([0, 1.6])
    ax1.set_ylabel('Normalized Responses', **hfont, fontsize=15)

    ax1.legend()
    plt.show()


# bar plot for individual neurons (PN)
for filler in range(1):
    hfont = {'fontname': 'Arial'}

    # find unit 230626_71 (reference unit)
    j = np.where(masterGoodUnits == '230627_71')[0][0]

    prefCenterBar = prefNormalized[j][0]
    barResponsesPN = {'PN measured': [pnNormalized[j][0], pnNormalized[j][3]],
                      'PN predicted': [(prefNormalized[j][0] + nonprefNormalized[j][1]) / 2,
                                       (prefNormalized[j][0] + nonprefNormalized[j][4]) / 2],
                      'Non-pref offset': [nonprefNormalized[j][1], nonprefNormalized[j][3]]}

    fig, ax1 = plt.subplots(figsize=(7, 7), layout='constrained')

    width = 0.25
    colorList = ['tab:green', 'tab:green', 'grey']
    cluster_positions = [0, 1, 3]
    x_single = [cluster_positions[0]]  # For Position 0
    x_multiple = [cluster_positions[1], cluster_positions[2]]  # For Positions 1 and 4

    ax1.bar(x_single, prefCenterBar, width, label="Pref Center", color="black")
    multiplier = 0
    for attribute, measurement in barResponsesPN.items():
        for i, x_pos in enumerate(x_multiple):
            offset = width * multiplier
            if multiplier == 1:
                rects = ax1.bar(x_pos + offset, measurement[i], width, label=attribute if i == 0 else "",
                                color=colorList[multiplier], hatch='//', alpha=0.5)
            else:
                rects = ax1.bar(x_pos + offset, measurement[i], width, label=attribute if i == 0 else "",
                                color=colorList[multiplier])
        multiplier += 1

    # Calculate the center of each cluster
    tick_positions = [cluster_positions[0],
                      cluster_positions[1] + width,
                      cluster_positions[2] + width]

    # Set custom x-axis tick labels at the center of each cluster
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(['0', '1', '4'], fontsize=20, **hfont)
    ax1.set_xlabel('Receptive Field Location', **hfont, fontsize=25)
    ax1.set_ylim([0, 1.6])
    ax1.set_ylabel('Response Rate (spikes/s)', **hfont, fontsize=25)
    ax1.legend()
    y_ticks = ax1.get_yticks()
    # Create a list of y-tick labels with only 0 and 1.0 labeled
    y_tick_labels = [
        '0' if tick == 0 else
        f'{round(normVal[j][0])}' if tick == 1.0 else
        ''  # blank for other ticks
        for tick in y_ticks
    ]
    ax1.set_yticklabels(y_tick_labels, fontsize=20, **hfont)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.set_xlim([-0.5, 4])
    ax1.tick_params(axis='both', width=2, length=8)

    plt.show()
    # # Save as a PDF with tight bounding box and high DPI
    # filePath = f'/Users/chery/Documents/Grad School/Maunsell Lab/Manuscripts/Python Figures/{unitList[j]}.pdf'
    # fig.savefig(filePath, format="pdf", bbox_inches='tight', dpi=300)
    # plt.close('all')


# scatter of real vs predicted responses at position 1 and 4



plt.scatter([1, 4], [prefMean[1], prefMean[3]], color='black')
plt.scatter([1, 4], pSingleFitPred[:2], facecolors='none', edgecolors='black')
plt.scatter([1, 4], [nonprefMean[1], nonprefMean[3]], color='gray')
plt.scatter([1, 4], npSingleFitPred[:2], facecolors='none', edgecolors='gray')
plt.scatter([1, 4], [pnMean[1], pnMean[3]], color='green')
plt.scatter([1, 4], pnFitPred[:2], facecolors='none', edgecolors='green')

plt.xlim([0, 5])
plt.show()



########################################################################################################################

############################################################ EMS #######################################################

########################################################################################################################


def fullEMSNorm(contrast_center, contrast_periphery, location, stim_type_center,
                stim_type_periphery, Lp, Lnp, W0, W1, W2, W3, W4, sigma):
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
    center = (contrast_center * L_center) / (contrast_center + (contrast_periphery * W_periphery) + sigma)
    periphery = (contrast_periphery * L_periphery) / ((contrast_center * W0) + contrast_periphery + sigma)

    return center + periphery


# Wrapper function to pass all data element-wise to curve_fit
def model_wrapperEMSNorm(data, Lp, Lnp, W0, W1, W2, W3, W4, sigma):
    contrast_center, contrast_periphery, location, stim_type_center, stim_type_periphery = data
    # Apply the response model element-wise
    return [fullEMSNorm(c_center, c_periph, loc, stim_c, stim_p, Lp, Lnp, W0, W1, W2, W3, W4, sigma)
            for c_center, c_periph, loc, stim_c, stim_p in
            zip(contrast_center, contrast_periphery, location, stim_type_center, stim_type_periphery)]


# Define the function to calculate predicted responses using the fitted parameters
def apply_fitted_modelEMS(contrast_center, contrast_periphery, location, stim_type_center,
                          stim_type_periphery, Lp, Lnp, W0, W1, W2, W3, W4, sigma):
    # Same logic as in response_model, but handles arrays
    predicted_responses = []
    for c_center, c_periph, loc, stim_c, stim_p in zip(contrast_center, contrast_periphery, location,
                                                       stim_type_center, stim_type_periphery):
        # Call the response model with the fitted parameters for each data point
        pred = fullEMSNorm(c_center, c_periph, loc, stim_c, stim_p, Lp, Lnp, W0, W1, W2, W3, W4, sigma)
        predicted_responses.append(pred)
    return predicted_responses


# known variables
contrast_center = [1, 0, 0, 0, 0,
                   1, 0, 0, 0, 0,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1,
                   1, 1, 1, 1]  # First response has contrast in the center, second in periphery
contrast_periphery = [0, 1, 1, 1, 1,
                      0, 1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1,
                      1, 1, 1, 1]  # First has no periphery, second has periphery with contrast
locations = np.array([-1, 0, 1, 2, 3,
                      -1, 0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3])  # First has no peripheral stimulus, second is at location 1
stim_type_center = np.array([1, -1, -1, -1, -1,
                             0, -1, -1, -1, -1,
                             1, 1, 1, 1,
                             0, 0, 0, 0,
                             1, 1, 1, 1,
                             0, 0, 0, 0])
stim_type_periphery = np.array([-1, 1, 1, 1, 1,
                                -1, 0, 0, 0, 0,
                                0, 0, 0, 0,
                                1, 1, 1, 1,
                                1, 1, 1, 1,
                                0, 0, 0, 0])


# individual neurons fit
fullEMSNormPred = []
measuredResp = []
r2ScoresEMSNormFull = []
hfont = {'fontname': 'Arial'}
for i in range(len(prefNormalized)):
    resp = np.concatenate((prefNormalized[i], nonprefNormalized[i], pnNormalized[i],
                           npNormalized[i], ppNormalized[i], nnNormalized[i]), axis=0)

    initial_guess = [1.0, 0.5, 1.0, 1.0, 0.75, 0.50, .10, 0.10]  # Lp, Lnp, W1, W2, W3, W4, sigma

    popt, pcov = curve_fit(model_wrapperEMSNorm, (contrast_center, contrast_periphery, locations, stim_type_center,
                                                  stim_type_periphery), resp,
                           p0=initial_guess)
                           # bounds=[[0, 0, 0, 0, 0, 0, 0, 0],
                           #         [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]])
    y_pred1 = apply_fitted_modelEMS(contrast_center, contrast_periphery, locations,
                                    stim_type_center, stim_type_periphery, *popt)
    r2 = r2_score(resp, y_pred1)
    r2ScoresEMSNormFull.append(r2)

    # single presentation prediction
    pCenterSinglePred = np.array(fullEMSNorm(1, 0, -1, 1, -1, *popt))
    pCenterSinglePred = pCenterSinglePred[np.newaxis]
    pPeriSinglePred = np.array([fullEMSNorm(0, 1, j, -1, 1, *popt) for j in range(4)])
    pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

    npCenterSinglePred = np.array(fullEMSNorm(1, 0, -1, 0, -1, *popt))
    npCenterSinglePred = npCenterSinglePred[np.newaxis]
    npPeriSinglePred = np.array([fullEMSNorm(0, 1, j, -1, 0, *popt) for j in range(4)])
    npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

    # paired presentation prediction
    pnFitPred = [fullEMSNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
    npFitPred = [fullEMSNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
    ppFitPred = [fullEMSNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
    nnFitPred = [fullEMSNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

    predResp = np.concatenate((pSingleFitPred, npSingleFitPred, pnFitPred,
                               npFitPred, ppFitPred, nnFitPred), axis=0)

    for j in range(len(resp)):
        fullEMSNormPred.append(predResp[j])
        measuredResp.append(resp[j])

# population average fit
resp = np.concatenate((prefMean, nonprefMean, pnMean, npMean, ppMean, nnMean), axis=0)
initial_guess = [1.0, 0.5, 1.0, 1.0, 0.75, 0.50, .10, 0.10]  # Lp, Lnp, W1, W2, W3, W4, sigma
popt, pcov = curve_fit(model_wrapperEMSNorm, (contrast_center, contrast_periphery,
                                                  locations, stim_type_center,
                                                  stim_type_periphery), resp,
                       p0=initial_guess)
                       # bounds=[[0, 0, 0, 0, 0, 0, 0, 0],
                       #         [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]])

y_pred1 = apply_fitted_modelEMS(contrast_center, contrast_periphery, locations,
                                stim_type_center, stim_type_periphery, *popt)
r2 = r2_score(resp, y_pred1)

pnFitPred = [fullEMSNorm(1, 1, j, 1, 0, *popt) for j in range(4)]
npFitPred = [fullEMSNorm(1, 1, j, 0, 1, *popt) for j in range(4)]
ppFitPred = [fullEMSNorm(1, 1, j, 1, 1, *popt) for j in range(4)]
nnFitPred = [fullEMSNorm(1, 1, j, 0, 0, *popt) for j in range(4)]

pCenterSinglePred = np.array(fullEMSNorm(1, 0, -1, 1, -1, *popt))
pCenterSinglePred = pCenterSinglePred[np.newaxis]
pPeriSinglePred = np.array([fullEMSNorm(0, 1, j, -1, 1, *popt) for j in range(4)])
pSingleFitPred = np.concatenate((pCenterSinglePred, pPeriSinglePred), axis=0)

npCenterSinglePred = np.array(fullEMSNorm(1, 0, -1, 0, -1, *popt))
npCenterSinglePred = npCenterSinglePred[np.newaxis]
npPeriSinglePred = np.array([fullEMSNorm(0, 1, j, -1, 0, *popt) for j in range(4)])
npSingleFitPred = np.concatenate((npCenterSinglePred, npPeriSinglePred), axis=0)

# plotting
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.set_size_inches(16, 6)
popt_str = (f"Lp = {popt[0]:.2f}\n"
            f"Lnp = {popt[1]:.2f}\n"
            f"W0 = {popt[2]:.2f}\n"
            f"W1 = {popt[3]:.2f}\n"
            f"W2 = {popt[4]:.2f}\n"
            f"W3 = {popt[5]:.2f}\n"
            f"W4 = {popt[6]:.2f}\n"
            f"sigma = {popt[7]:.2f}")

ax1.plot(x, prefMean, color='black', label='PO')
ax1.errorbar(x, prefMean, yerr=prefSEM, fmt='o', ecolor='black',
             color='black', markersize=7)
ax1.plot(x, nonprefMean, color='grey', label='NO')
ax1.errorbar(x, nonprefMean, yerr=nonprefSEM, fmt='o', ecolor='grey',
             color='grey', markersize=7)
ax1.plot(xNorm, pnMean, color='green', label='P0 N1')
ax1.errorbar(xNorm, pnMean, yerr=pnSEM, fmt='o', ecolor='green',
             color='green', markersize=7)
ax1.plot(xNorm, npMean, color='red', label='N0 P1')
ax1.errorbar(xNorm, npMean, yerr=npSEM, fmt='o', ecolor='red',
             color='red', markersize=7)
ax1.plot(xNorm, ppMean, color='black', label='P0 P1')
ax1.errorbar(xNorm, ppMean, yerr=pnSEM, fmt='o', ecolor='black',
             color='black', markersize=7)
ax1.plot(xNorm, nnMean, color='grey', label='N0 N1')
ax1.errorbar(xNorm, nnMean, yerr=npSEM, fmt='o', ecolor='grey',
             color='grey', markersize=7)
ax1.plot(xNorm, pnFitPred, color='green', linestyle='--')
ax1.scatter(xNorm, pnFitPred, facecolors='none', edgecolors='g', s=50)
ax1.plot(xNorm, npFitPred, color='red', linestyle='--')
ax1.scatter(xNorm, npFitPred, facecolors='none', edgecolors='r', s=50)
ax1.plot(xNorm, ppFitPred, color='black', linestyle=':', alpha=0.8)
ax1.scatter(xNorm, ppFitPred, facecolors='none', edgecolors='black', s=50)
ax1.plot(xNorm, nnFitPred, color='grey', linestyle=':', alpha=0.8)
ax1.scatter(xNorm, nnFitPred, facecolors='none', edgecolors='gray', s=50)
ax1.plot(x, pSingleFitPred, color='black', linestyle='--')
ax1.scatter(x, pSingleFitPred, facecolors='none', edgecolors='black', s=50)
ax1.plot(x, npSingleFitPred, color='grey', linestyle='--')
ax1.scatter(x, npSingleFitPred, facecolors='none', edgecolors='gray', s=50)
ax1.set_xlabel('Receptive Field Offset', fontsize=15)
ax1.set_ylabel('Normalized Response Rate', fontsize=15)
ax1.set_ylim([0, 1.50])
ax1.set_title('Normalization vs Distance PN/NP')
ax1.text(0.95, 0.95, popt_str, fontsize=10, verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightyellow'),
         transform=ax1.transAxes)

ax2.scatter(measuredResp, fullEMSNormPred, color='black', facecolors='none')
ax2.set_ylim([0, 2])
ax2.set_xlim([0, 2])
ax2.set_xlabel('Measured Response', fontsize=15, **hfont)
ax2.set_ylabel('Predicted Response', fontsize=15, **hfont)
line = lines.Line2D([0, 1], [0, 1], color='black')
transform = ax2.transAxes
line.set_transform(transform)
ax2.add_line(line)

ax3.hist(r2ScoresEMSNormFull, bins=5)
ax3.set_xlim([0, 1])
ax3.set_xlabel('Explained R2', fontsize=15, **hfont)
ax3.set_ylabel('Frequency', fontsize=15, **hfont)

plt.tight_layout()
plt.show()


#############

# Create a figure and plot the violin plot
plt.figure(figsize=(8, 6))
sns.violinplot(data=corrLists)

# Set x-axis labels
plt.xticks(ticks=range(4), labels=['1', '2', '3', '4'])
plt.xlabel('List Number')
plt.ylabel('Values')
plt.title('Violin Plot of corrLists')

# Show plot
plt.show()


# Calculate the medians for each list
medians = [np.median(lst) for lst in corrLists]
# Create a figure and plot the box plot
plt.figure(figsize=(8, 6))
sns.boxplot(data=corrLists)

# Set x-axis labels
plt.xticks(ticks=range(4), labels=['1', '2', '3', '4'])
plt.xlabel('List Number')
plt.ylabel('Values')
plt.title('Box Plot of corrLists with Medians')

# Annotate the median for each box plot
for i, median in enumerate(medians):
    plt.text(i, median, f'{median:.2f}', ha='center', va='top', color='black', fontweight='bold')

# Show plot
plt.show()
