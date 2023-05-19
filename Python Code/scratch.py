import numpy as np
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Plot, compute arrays
ppArr = np.array(ppArr)
pnArr = np.array(pnArr)
npArr = np.array(npArr)
nnArr = np.array(nnArr)
p0Arr = np.array(p0Arr)
p1Arr = np.array(p1Arr)
n0Arr = np.array(n0Arr)
n1Arr = np.array(n1Arr)
baseArr = np.array(baseArr)
popRespHeatmap = np.array(popRespHeatmap)

meanPPArr = np.mean(ppArr, 0)
semPPArr = stats.sem(ppArr, 0)
meanPNArr = np.mean(pnArr, 0)
semPNArr = stats.sem(pnArr, 0)
meanNPArr = np.mean(npArr, 0)
semNPArr = stats.sem(npArr, 0)
meanNNArr = np.mean(nnArr, 0)
semNNArr = stats.sem(nnArr, 0)
meanBaseArr = np.mean(baseArr, 0)
semBaseArr = stats.sem(baseArr, 0)

# loc 0 P, loc 1 N
meanP0Arr = np.mean(p0Arr, 0)
semP0Arr = stats.sem(p0Arr, 0)
meanN1Arr = np.mean(n1Arr, 0)
semN1Arr = stats.sem(n1Arr, 0)

# loc 0 N, loc 1 P
meanP1Arr = np.mean(p1Arr, 0)
semP1Arr = stats.sem(p1Arr, 0)
meanN0Arr = np.mean(n0Arr, 0)
semN0Arr = stats.sem(n0Arr, 0)

smoothPlotPP = gaussian_filter1d(meanPPArr, 5)
semSmoothPP = gaussian_filter1d(semPPArr, 5)
smoothPlotPN = gaussian_filter1d(meanPNArr, 5)
semSmoothPN = gaussian_filter1d(semPNArr, 5)
smoothPlotNP = gaussian_filter1d(meanNPArr, 5)
semSmoothNP = gaussian_filter1d(semNPArr, 5)
smoothPlotNN = gaussian_filter1d(meanNNArr, 5)
semSmoothNN = gaussian_filter1d(semNNArr, 5)
smoothPlotP0 = gaussian_filter1d(meanP0Arr, 5)
semSmoothP0 = gaussian_filter1d(semP0Arr, 5)
smoothPlotP1 = gaussian_filter1d(meanP1Arr, 5)
semSmoothP1 = gaussian_filter1d(semP1Arr, 5)
smoothPlotN0 = gaussian_filter1d(meanN0Arr, 5)
semSmoothN0 = gaussian_filter1d(semN0Arr, 5)
smoothPlotN1 = gaussian_filter1d(meanN1Arr, 5)
semSmoothN1 = gaussian_filter1d(semN1Arr, 5)
smoothBase = gaussian_filter1d(meanBaseArr, 5)
semSmoothBase = gaussian_filter1d(semBaseArr, 5)

# figure
x = np.arange(0, 448)

fig = plt.figure()
fig.set_size_inches(15, 6)
gs0 = gridspec.GridSpec(1, 3)
gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0])
gs01 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[2])

ax = fig.add_subplot(gs00[0, 0])
ax.plot(smoothPlotPP, linestyle='solid', color='black', label='PP')
ax.fill_between(x, smoothPlotPP - semSmoothPP, smoothPlotPP + semSmoothPP,
                alpha=0.1, color='black')
ax.plot(smoothPlotPN, linestyle='dashed', color='green', label='PN')
ax.fill_between(x, smoothPlotPN - semSmoothPN, smoothPlotPN + semSmoothPN,
                alpha=0.1, color='green')
ax.plot(smoothPlotNP, linestyle='dashed', color='blue', label='NP')
ax.fill_between(x, smoothPlotNP - semSmoothNP, smoothPlotNP + semSmoothNP,
                alpha=0.1, color='blue')
ax.plot(smoothPlotNN, linestyle='solid', color='grey', label='NN')
ax.fill_between(x, smoothPlotNN - semSmoothNN,
                smoothPlotNN + semSmoothNN, alpha=0.1, color='grey')
ax.plot(smoothBase, linestyle='solid', color='black', label='Baseline', linewidth=0.5)
ax.fill_between(x, smoothBase - semSmoothBase,
                smoothBase + semSmoothBase, alpha=0.1, color='black')
ax.set_xticks([0,
               histPrePostMS, 150,
               histPrePostMS + trueStimDurMS,
               2 * histPrePostMS + trueStimDurMS])
ax.set_xticklabels([])
ax.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
           color='grey', alpha=0.3)
ax.set_ylabel('Normalized Firing Rate (spikes/sec)', fontsize=15)
ax.set_xlabel('Stimulus Duration (ms)', fontsize=15)
ax.set_xticklabels([-histPrePostMS, 0, 50, 0 + trueStimDurMS,
                    trueStimDurMS + histPrePostMS],
                    fontsize=7)
ax.set_ylim([0, 1])
ax.legend()

ax2 = fig.add_subplot(gs01[0, 0])
ax2.plot(smoothPlotP0, linestyle='solid', color='black', label='P0')
ax2.fill_between(x, smoothPlotP0 - semSmoothP0, smoothPlotP0 + semSmoothP0,
                 alpha=0.1, color='black')
ax2.plot(smoothPlotPN, linestyle='dashed', color='green', label='PN')
ax2.fill_between(x, smoothPlotPN - semSmoothPN, smoothPlotPN + semSmoothPN,
                 alpha=0.1, color='green')
ax2.plot(smoothPlotN1, linestyle='solid', color='grey', label='N1')
ax2.fill_between(x, smoothPlotN1 - semSmoothN1,
                 smoothPlotN1 + semSmoothN1, alpha=0.1, color='grey')
ax2.plot(smoothBase, linestyle='solid', color='black', label='Baseline', linewidth=0.5)
ax2.fill_between(x, smoothBase - semSmoothBase,
                 smoothBase + semSmoothBase, alpha=0.1, color='black')
ax2.set_xticks([0,
               histPrePostMS, 150,
               histPrePostMS + trueStimDurMS,
               2 * histPrePostMS + trueStimDurMS])
ax2.set_xticklabels([])
ax2.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
ax2.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
            color='grey', alpha=0.3)
ax2.set_xlabel('Stimulus Duration (ms)', fontsize=15)
ax2.set_xticklabels([-histPrePostMS, 0, 50, 0 + trueStimDurMS,
                    trueStimDurMS + histPrePostMS],
                    fontsize=7)
ax2.set_ylim([0, 1])
ax2.legend()

ax3 = fig.add_subplot(gs02[0, 0])
ax3.plot(smoothPlotP1, linestyle='solid', color='black', label='P1')
ax3.fill_between(x, smoothPlotP1 - semSmoothP1, smoothPlotP1 + semSmoothP1,
                 alpha=0.1, color='black')
ax3.plot(smoothPlotNP, linestyle='dashed', color='blue', label='NP')
ax3.fill_between(x, smoothPlotNP - semSmoothNP, smoothPlotNP + semSmoothNP,
                 alpha=0.1, color='green')
ax3.plot(smoothPlotN0, linestyle='solid', color='grey', label='N0')
ax3.fill_between(x, smoothPlotN0 - semSmoothN0,
                 smoothPlotN0 + semSmoothN0, alpha=0.1, color='grey')
ax3.plot(smoothBase, linestyle='solid', color='black', label='Baseline', linewidth=0.5)
ax3.fill_between(x, smoothBase - semSmoothBase,
                 smoothBase + semSmoothBase, alpha=0.1, color='black')
ax3.set_xticks([0,
               histPrePostMS, 150,
               histPrePostMS + trueStimDurMS,
               2 * histPrePostMS + trueStimDurMS])
ax3.set_xticklabels([])
ax3.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
ax3.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
            color='grey', alpha=0.3)
ax3.set_xlabel('Stimulus Duration (ms)', fontsize=15)
ax3.set_xticklabels([-histPrePostMS, 0, 50, 0 + trueStimDurMS,
                    trueStimDurMS + histPrePostMS],
                    fontsize=7)
ax3.legend()
ax3.set_ylim([0, 1])

fig.suptitle('gabor sigma sep > 25th percentile and < 50th gabor sigma sep', fontsize=20)
plt.show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# plot PSTH session by session
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
seshMeanA0 = []
seshMeanA1 = []
start = 0
for count, i in enumerate(totUnits):
    end = start + i

    meanPPArr = np.mean(ppArr[start:end, :], 0)
    semPPArr = stats.sem(ppArr[start:end, :], 0)
    meanPNArr = np.mean(pnArr[start:end, :], 0)
    semPNArr = stats.sem(pnArr[start:end, :], 0)
    meanNPArr = np.mean(npArr[start:end, :], 0)
    semNPArr = stats.sem(npArr[start:end, :], 0)
    meanNNArr = np.mean(nnArr[start:end, :], 0)
    semNNArr = stats.sem(nnArr[start:end, :], 0)
    meanBaseArr = np.mean(baseArr[start:end, :], 0)
    semBaseArr = stats.sem(baseArr[start:end, :], 0)
    # loc 0 P, loc 1 N
    meanP0Arr = np.mean(p0Arr[start:end, :], 0)
    semP0Arr = stats.sem(p0Arr[start:end, :], 0)
    meanN1Arr = np.mean(n1Arr[start:end, :], 0)
    semN1Arr = stats.sem(n1Arr[start:end, :], 0)
    # loc 0 N, loc 1 P
    meanP1Arr = np.mean(p1Arr[start:end, :], 0)
    semP1Arr = stats.sem(p1Arr[start:end, :], 0)
    meanN0Arr = np.mean(n0Arr[start:end, :], 0)
    semN0Arr = stats.sem(n0Arr[start:end, :], 0)

    # session average alpha at loc 0 and loc 1
    meanA0 = np.mean(alphaLoc0PrefOnly[start:end])
    seshMeanA0.append(meanA0)
    meanA1 = np.mean(alphaLoc1PrefOnly[start:end])
    seshMeanA1.append(meanA1)

    # update start position for next session
    start = end

    smoothPlotPP = gaussian_filter1d(meanPPArr, 5)
    semSmoothPP = gaussian_filter1d(semPPArr, 5)
    smoothPlotPN = gaussian_filter1d(meanPNArr, 5)
    semSmoothPN = gaussian_filter1d(semPNArr, 5)
    smoothPlotNP = gaussian_filter1d(meanNPArr, 5)
    semSmoothNP = gaussian_filter1d(semNPArr, 5)
    smoothPlotNN = gaussian_filter1d(meanNNArr, 5)
    semSmoothNN = gaussian_filter1d(semNNArr, 5)
    smoothPlotP0 = gaussian_filter1d(meanP0Arr, 5)
    semSmoothP0 = gaussian_filter1d(semP0Arr, 5)
    smoothPlotP1 = gaussian_filter1d(meanP1Arr, 5)
    semSmoothP1 = gaussian_filter1d(semP1Arr, 5)
    smoothPlotN0 = gaussian_filter1d(meanN0Arr, 5)
    semSmoothN0 = gaussian_filter1d(semN0Arr, 5)
    smoothPlotN1 = gaussian_filter1d(meanN1Arr, 5)
    semSmoothN1 = gaussian_filter1d(semN1Arr, 5)
    smoothBase = gaussian_filter1d(meanBaseArr, 5)
    semSmoothBase = gaussian_filter1d(semBaseArr, 5)

    # figure
    x = np.arange(0, 448)

    fig = plt.figure()
    fig.set_size_inches(15, 6)
    gs0 = gridspec.GridSpec(1, 3)
    gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0])
    gs01 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1])
    gs02 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[2])

    ax = fig.add_subplot(gs00[0, 0])
    ax.plot(smoothPlotPP, linestyle='solid', color='black', label='PP')
    ax.fill_between(x, smoothPlotPP - semSmoothPP, smoothPlotPP + semSmoothPP,
                    alpha=0.1, color='black')
    ax.plot(smoothPlotPN, linestyle='dashed', color='green', label='PN')
    ax.fill_between(x, smoothPlotPN - semSmoothPN, smoothPlotPN + semSmoothPN,
                    alpha=0.1, color='green')
    ax.plot(smoothPlotNP, linestyle='dashed', color='blue', label='NP')
    ax.fill_between(x, smoothPlotNP - semSmoothNP, smoothPlotNP + semSmoothNP,
                    alpha=0.1, color='blue')
    ax.plot(smoothPlotNN, linestyle='solid', color='grey', label='NN')
    ax.fill_between(x, smoothPlotNN - semSmoothNN,
                    smoothPlotNN + semSmoothNN, alpha=0.1, color='grey')
    ax.plot(smoothBase, linestyle='solid', color='black', label='Baseline', linewidth=0.5)
    ax.fill_between(x, smoothBase - semSmoothBase,
                    smoothBase + semSmoothBase, alpha=0.1, color='black')
    ax.set_xticks([0,
                   histPrePostMS, 150,
                   histPrePostMS + trueStimDurMS,
                   2 * histPrePostMS + trueStimDurMS])
    ax.set_xticklabels([])
    ax.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
    ax.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
               color='grey', alpha=0.3)
    ax.set_ylabel('Normalized Firing Rate (spikes/sec)', fontsize=15)
    ax.set_xlabel('Stimulus Duration (ms)', fontsize=15)
    ax.set_xticklabels([-histPrePostMS, 0, 50, 0 + trueStimDurMS,
                        trueStimDurMS + histPrePostMS],
                       fontsize=7)
    ax.set_ylim([0, 1])
    ax.legend()

    ax2 = fig.add_subplot(gs01[0, 0])
    ax2.plot(smoothPlotP0, linestyle='solid', color='black', label='P0')
    ax2.fill_between(x, smoothPlotP0 - semSmoothP0, smoothPlotP0 + semSmoothP0,
                     alpha=0.1, color='black')
    ax2.plot(smoothPlotPN, linestyle='dashed', color='green', label='PN')
    ax2.fill_between(x, smoothPlotPN - semSmoothPN, smoothPlotPN + semSmoothPN,
                     alpha=0.1, color='green')
    ax2.plot(smoothPlotN1, linestyle='solid', color='grey', label='N1')
    ax2.fill_between(x, smoothPlotN1 - semSmoothN1,
                     smoothPlotN1 + semSmoothN1, alpha=0.1, color='grey')
    ax2.plot(smoothBase, linestyle='solid', color='black', label='Baseline', linewidth=0.5)
    ax2.fill_between(x, smoothBase - semSmoothBase,
                     smoothBase + semSmoothBase, alpha=0.1, color='black')
    ax2.set_xticks([0,
                    histPrePostMS, 150,
                    histPrePostMS + trueStimDurMS,
                    2 * histPrePostMS + trueStimDurMS])
    ax2.set_xticklabels([])
    ax2.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
    ax2.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                color='grey', alpha=0.3)
    ax2.set_xlabel('Stimulus Duration (ms)', fontsize=15)
    ax2.set_xticklabels([-histPrePostMS, 0, 50, 0 + trueStimDurMS,
                         trueStimDurMS + histPrePostMS],
                        fontsize=7)
    ax2.set_ylim([0, 1])
    ax2.legend()

    ax3 = fig.add_subplot(gs02[0, 0])
    ax3.plot(smoothPlotP1, linestyle='solid', color='black', label='P1')
    ax3.fill_between(x, smoothPlotP1 - semSmoothP1, smoothPlotP1 + semSmoothP1,
                     alpha=0.1, color='black')
    ax3.plot(smoothPlotNP, linestyle='dashed', color='blue', label='NP')
    ax3.fill_between(x, smoothPlotNP - semSmoothNP, smoothPlotNP + semSmoothNP,
                     alpha=0.1, color='green')
    ax3.plot(smoothPlotN0, linestyle='solid', color='grey', label='N0')
    ax3.fill_between(x, smoothPlotN0 - semSmoothN0,
                     smoothPlotN0 + semSmoothN0, alpha=0.1, color='grey')
    ax3.plot(smoothBase, linestyle='solid', color='black', label='Baseline', linewidth=0.5)
    ax3.fill_between(x, smoothBase - semSmoothBase,
                     smoothBase + semSmoothBase, alpha=0.1, color='black')
    ax3.set_xticks([0,
                    histPrePostMS, 150,
                    histPrePostMS + trueStimDurMS,
                    2 * histPrePostMS + trueStimDurMS])
    ax3.set_xticklabels([])
    ax3.set_xlim([0, trueStimDurMS + (2 * histPrePostMS + 1)])
    ax3.axvspan(histPrePostMS, histPrePostMS + trueStimDurMS,
                color='grey', alpha=0.3)
    ax3.set_xlabel('Stimulus Duration (ms)', fontsize=15)
    ax3.set_xticklabels([-histPrePostMS, 0, 50, 0 + trueStimDurMS,
                         trueStimDurMS + histPrePostMS],
                        fontsize=7)
    ax3.legend()
    ax3.set_ylim([0, 1])

    fig.suptitle(f'Session {count+1} average PSTH', fontsize=20)
    plt.savefig(f'../Analysis Screenshots/session{count+1}average.pdf')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# plt scatter of session mean alphas
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
fig, ax = plt.subplots()
ax.scatter(seshMeanA0, seshMeanA1)
ax.set_ylabel('session mean Alpha 1')
ax.set_xlabel('session mean Alpha 0')
ax.set_title('scatter of mean A0 vs A1 for each session')
ax.set_xlim(left=0, right=3.5)
ax.set_ylim(bottom=0, top=3.5)
ax.set_aspect('equal', adjustable='box')
plt.show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# population normalized heatmap
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
tickLabels = np.array(['-180', '-120', '-60', '0', '60', '120', 'blank'])

fig, ax = plt.subplots()
meanPopResp = np.mean(popRespHeatmap, 0)
meanPopResp = meanPopResp.reshape(7, 7)
ax = sns.heatmap(meanPopResp, square=True, linewidths=0.2, vmin=0,
                 vmax=np.max(meanPopResp), annot=True, annot_kws={'fontsize': 13})

ax.set_xticks(np.arange(7) + 0.5)
ax.set_title(f'Population Normalized Response: Aligned to Preferred',
             y=-0.1, fontsize=15)
ax.set_xlabel('Location 0 Stimulus Direction', fontsize=15)
ax.xaxis.set_label_position('top')
ax.set_xticklabels(tickLabels, rotation=45, fontsize=15)
ax.set_ylabel('Location 1 Stimulus Direction', fontsize=15)
ax.xaxis.set_ticks_position("top")
ax.set_yticks(np.arange(7) + 0.5)
ax.set_yticklabels(tickLabels, rotation=0, fontsize=15)
plt.show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# binned plot of gabor distance from center of RF and mean alpha at that loc
# also can plot pref resp (normalized to max resp of that neuron)
# as a fn of distance from RF center
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
loc0DistFromRFCent = np.array(loc0DistFromRFCent)
loc1DistFromRFCent = np.array(loc1DistFromRFCent)
loc0PrefNormalized = np.array(loc0PrefNormalized)
loc1PrefNormalized = np.array(loc1PrefNormalized)
alphaLoc0 = np.array(alphaLoc0)
alphaLoc1 = np.array(alphaLoc1)
loc0NullNMI = np.array(loc0NullNMI)
loc1NullNMI = np.array(loc1NullNMI)
combDistFromRFCent = np.array([loc0DistFromRFCent, loc1DistFromRFCent]).flatten()
combAlpha = np.array([alphaLoc0, alphaLoc1]).flatten()
combPrefNormalized = np.array([loc0PrefNormalized, loc1PrefNormalized]).flatten()
combNMI = np.array([loc0NullNMI, loc1NullNMI]).flatten()

# exclude elements that are > 10 deg from the mean
filt = np.where(combDistFromRFCent < 5)[0]
filtDistFromRFCent = combDistFromRFCent[filt]
filtAlpha = combAlpha[filt]
filtPref = combPrefNormalized[filt]
filtNMI = combNMI[filt]

# exclude alphas > 4
filt = np.where(filtAlpha < 4)[0]
filtDistFromRFCent = filtDistFromRFCent[filt]
filtAlpha = filtAlpha[filt]
filtPref = filtPref[filt]
filtNMI = filtNMI[filt]

sortIndex = np.argsort(filtDistFromRFCent)
sortedDist = filtDistFromRFCent[sortIndex]
sortedAlpha = filtAlpha[sortIndex]
sortedPrefResp = filtPref[sortIndex]
sortedNMI = filtNMI[sortIndex]

# manual bins - equally populated bins
n = 10
equalBinsDist = [sortedDist[i:i + n] for i in range(0, len(sortedDist), n)]
equalBinsAlpha = [sortedAlpha[i:i + n] for i in range(0, len(sortedAlpha), n)]
binMeanDist = np.array([np.mean(i) for i in equalBinsDist])
binMeanAlpha = np.array([np.mean(i) for i in equalBinsAlpha])

# polynomial fit
a, b = np.polyfit(binMeanDist, binMeanAlpha, 1)

# figure
fig, ax = plt.subplots()
ax.scatter(binMeanDist, binMeanAlpha)
ax.plot(binMeanDist, (a*binMeanDist+b))
ax.set_xlabel('Gabor distance from RF center')
ax.set_ylabel('Normalization Alpha')
ax.set_xlim([0, 5])
ax.set_ylim(bottom=0)
ax.set_title(f'Equal sized bins ({n}) of gabor distance from RF center vs alpha')
plt.show()

# manual binning option 2 (defined by number of bins)
binSize = 10
binMedian, binEdges, binNum = stats.binned_statistic(sortedDist, sortedAlpha,
                                                     statistic='median', bins=binSize)
binCount, binEdges, binNum = stats.binned_statistic(sortedDist, sortedAlpha,
                                                    statistic='count', bins=binSize)
binSD, binEdges, binNum = stats.binned_statistic(sortedDist, sortedAlpha,
                                                 statistic='std', bins=binSize)
binSEM = binSD/binCount
binWidth = (binEdges[1] - binEdges[0])
binCenters = binEdges[1:] - binWidth/2

# polynomial fit
a, b = np.polyfit(binCenters, binMedian, 1)

# figure
fig, ax = plt.subplots()
ax.scatter(binCenters, binMedian)
ax.plot(binCenters, (a*binCenters+b))
ax.fill_between(binCenters, binMedian - binSEM, binMedian + binSEM,
                 alpha=0.1, color='black')
ax.set_xlabel('Gabor distance from RF center')
ax.set_ylabel('Normalization Alpha')
ax.set_ylim([0, np.max(binMedian)*1.2])
ax.set_xlim([0, 5])
ax.set_title('Binned plot of gabor distance from RF center vs alpha')
plt.show()

# plot use bin reg
est = binsreg(filtAlpha, filtDistFromRFCent, line=(5, 5), cb=(5, 5), nbins=5)

# figure for pref resp normalized as fn of distance from RF center
# manual bins - equally populated bins
n = 10
equalBinsDist = [sortedDist[i:i + n] for i in range(0, len(sortedDist), n)]
equalBinsPref = [sortedPrefResp[i:i + n] for i in range(0, len(sortedPrefResp), n)]
binMeanDist = np.array([np.mean(i) for i in equalBinsDist])
binMeanPref = np.array([np.mean(i) for i in equalBinsPref])

# polynomial fit
a, b = np.polyfit(binMeanDist, binMeanPref, 1)

# figure
fig, ax = plt.subplots()
ax.scatter(binMeanDist, binMeanPref)
ax.plot(binMeanDist, (a*binMeanDist+b))
ax.set_xlabel('Gabor distance from RF center')
ax.set_ylabel('Pref resp normalized to stim with strongest resp')
ax.set_xlim([0, 5])
ax.set_ylim([0, 1])
ax.set_title(f'Equal sized bins ({n}) of gabor distance from RF center vs pref resp normalized')
plt.show()

# figure for NMI as a function of distance from RF center
# manual bins - equally populated bins
n = 20
equalBinsDist = [sortedDist[i:i + n] for i in range(0, len(sortedDist), n)]
equalBinsNMI = [sortedNMI[i:i + n] for i in range(0, len(sortedNMI), n)]
binMeanDist = np.array([np.mean(i) for i in equalBinsDist])
binMeanNMI = np.array([np.mean(i) for i in equalBinsNMI])

# polynomial fit
a, b = np.polyfit(binMeanDist, binMeanNMI, 1)

# figure
fig, ax = plt.subplots()
ax.scatter(binMeanDist, binMeanNMI)
ax.plot(binMeanDist, (a*binMeanDist+b))
ax.set_xlabel('Gabor distance from RF center')
ax.set_ylabel('NMI')
ax.set_xlim([0, 5])
ax.set_ylim(bottom=0)
ax.set_title(f'Equal sized bins ({n}) of gabor distance from RF center vs NMI')
plt.show()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# plot resp pref/null ratio vs alpha1/alpha0 ratio
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
alphaLoc0 = np.array(alphaLoc0)
alphaLoc1 = np.array(alphaLoc1)
loc0to1RespRatio = np.array(loc0to1RespRatio)
alpha1to0Ratio = alphaLoc0/alphaLoc1

filt = np.where(alpha1to0Ratio < 8.50)[0]
filteredAlphaRatio = alpha1to0Ratio[filt]
filteredRespRatio = loc0to1RespRatio[filt]

fig, ax = plt.subplots()
ax.scatter(filteredRespRatio, filteredAlphaRatio)
ax.set_ylim([0, 5])
plt.show()