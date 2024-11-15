from scipy.optimize import curve_fit
import numpy as np
from usefulFns import *


# Define a function to calculate the bootstrap c50 values
def bootstrap_c50(response_matrix, contrasts, n_bootstraps):
    # Generate all bootstrap samples in a single operation
    bootstrap_samples = response_matrix[
        np.random.choice(response_matrix.shape[0],
                         size=(n_bootstraps, response_matrix.shape[0]),
                         replace=True)
    ]

    # Calculate mean across the rows for each bootstrap sample
    mean_counts = np.mean(bootstrap_samples, axis=1)

    # Fit curve for each mean count sample
    c50_values = []
    for mean_count in mean_counts:
        try:
            # Attempt to fit the curve with a lower maxfev
            pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, mean_count,
                                bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
                                maxfev=1000)  # Reduced maxfev for faster fitting
            c50_values.append(pOpt[1])  # Append fitted c50
        except RuntimeError:
            # If curve_fit fails, append the large c50 value
            c50_values.append(1e6)

    return np.array(c50_values)



########## extra code for bootstrapping  ####

# potentialGoodUnits.append(f'{seshDate}_{unit}')
# # bootstrap to estimate CI on c50
# prefCentRawSpikeCounts = spikeCountMat[unitID][:blocksDone, :6] * 1000 / trueStimDurMS
# prefCentRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
#                                             prefCentRawSpikeCounts))
# prefCentResp = prefCentRawSpikeCountsWithBase - baselineResp
#
# nonprefCentRawSpikeCounts = spikeCountMat[unitID][:blocksDone, 6:12] * 1000 / trueStimDurMS
# nonprefCentRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
#                                                nonprefCentRawSpikeCounts))
# nonprefCentResp = nonprefCentRawSpikeCountsWithBase - baselineResp
#
# prefPeriRawSpikeCounts = spikeCountMat[unitID][:blocksDone, 12:18] * 1000 / trueStimDurMS
# prefPeriRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
#                                             prefPeriRawSpikeCounts))
# prefPeriResp = prefPeriRawSpikeCountsWithBase - baselineResp
#
# nonprefPeriRawSpikeCounts = spikeCountMat[unitID][:blocksDone, 18:] * 1000 / trueStimDurMS
# nonprefPeriRawSpikeCountsWithBase = np.hstack((np.full((blocksDone, 1), baselineResp),
#                                                nonprefPeriRawSpikeCounts))
# nonprefPeriResp = nonprefPeriRawSpikeCountsWithBase - baselineResp
#
# n_bootstraps = 500
#
# # c50bootstrapPrefCent = []
# # c50bootstrapNonprefCent = []
# # c50bootstrapPrefPeri = []
# # c50bootstrapNonprefPeri = []
# # # Bootstrap loop
# # for _ in range(n_bootstraps):
# #     # pref center
# #     # Sample rows with replacement to get a bootstrap sample
# #     bootstrap_sample = prefCentResp[
# #         np.random.choice(prefCentResp.shape[0],
# #                          size=prefCentResp.shape[0],
# #                          replace=True)]
# #     meanCount = np.mean(bootstrap_sample, axis=0)
# #     pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
# #                         bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
# #                         maxfev=100000)
# #     c50bootstrapPrefCent.append(pOpt[1])
# #
# #     # nonpref cent
# #     bootstrap_sample = nonprefCentResp[
# #         np.random.choice(nonprefCentResp.shape[0],
# #                          size=nonprefCentResp.shape[0],
# #                          replace=True)]
# #     meanCount = np.mean(bootstrap_sample, axis=0)
# #     pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
# #                         bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
# #                         maxfev=100000)
# #     c50bootstrapNonprefCent.append(pOpt[1])
# #
# #     # pref peri
# #     bootstrap_sample = prefPeriResp[
# #         np.random.choice(prefPeriResp.shape[0],
# #                          size=prefPeriResp.shape[0],
# #                          replace=True)]
# #     meanCount = np.mean(bootstrap_sample, axis=0)
# #     pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
# #                         bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
# #                         maxfev=100000)
# #     c50bootstrapPrefPeri.append(pOpt[1])
# #
# #     # nonpref peri
# #     bootstrap_sample = nonprefPeriResp[
# #         np.random.choice(nonprefPeriResp.shape[0],
# #                          size=nonprefPeriResp.shape[0],
# #                          replace=True)]
# #     meanCount = np.mean(bootstrap_sample, axis=0)
# #     pOpt, _ = curve_fit(contrastFnNoBaseline, contrasts, meanCount,
# #                         bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
# #                         maxfev=100000)
# #     c50bootstrapNonprefPeri.append(pOpt[1])
# #
# # c50bootstrapPrefCent = np.array(c50bootstrapPrefCent)
# # c50bootstrapNonprefCent = np.array(c50bootstrapNonprefCent)
# # c50bootstrapPrefPeri = np.array(c50bootstrapPrefPeri)
# # c50bootstrapNonprefPeri = np.array(c50bootstrapNonprefPeri)
#
# # Compute bootstrap c50s for each condition using vector approach
# c50bootstrapPrefCent = bootstrap_c50(prefCentResp, contrasts, n_bootstraps)
# c50bootstrapNonprefCent = bootstrap_c50(nonprefCentResp, contrasts, n_bootstraps)
# c50bootstrapPrefPeri = bootstrap_c50(prefPeriResp, contrasts, n_bootstraps)
# c50bootstrapNonprefPeri = bootstrap_c50(nonprefPeriResp, contrasts, n_bootstraps)
#
# print(f'done bootstrap {unit}')
# print(time.time() - t1)
#
# # Calculate confidence interval (e.g., 95%)
# # confidence_interval = np.percentile(c50bootstrapPrefCent, [2.5, 97.5])
#
# prefCentCI = np.percentile(c50bootstrapPrefCent, [2.5, 97.5])
# nonprefCentCI = np.percentile(c50bootstrapNonprefCent, [2.5, 97.5])
# prefPeriCI = np.percentile(c50bootstrapPrefPeri, [2.5, 97.5])
# nonprefPeriCI = np.percentile(c50bootstrapNonprefPeri, [2.5, 97.5])
#
# potentialGoodUnitsCI.append([prefCentCI, nonprefCentCI, prefPeriCI, nonprefPeriCI])
#
# randVar = 0
# threshValue = 30
# if prefCentCI[1] / prefCentCI[0] > threshValue:
#     randVar = 1
# if nonprefCentCI[1] / nonprefCentCI[0] > threshValue:
#     randVar = 1
# if prefPeriCI[1] / prefPeriCI[0] > threshValue:
#     randVar = 1
# if nonprefPeriCI[1] / nonprefPeriCI[0] > threshValue:
#     randVar = 1
#
# if randVar == 0:
#     print(f'unit meets threshold cutoff')






