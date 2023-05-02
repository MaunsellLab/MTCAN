
# get p, n, p+n psth for every unit - normalized and aligned to each unit's
# pref direction
sli = [0, 3, 6]
for unitCount, unit in enumerate(units):
    # find direction tested that is closest to the pref dir
    # and reindex around this so that it is in the middle of the grid
    prefDir, nullDir = dirClosestToPref(unitGaussMean[unitCount])
    nullDirIndex = np.where(dirArray == nullDir)[0][0]
    reIndex = (np.array([0, 1, 2, 3, 4, 5]) + nullDirIndex) % 6

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

    # condensed matrices for 2 directions only
    stimMatCond = stimMatReIndex[sli, :]
    stimMatCond = stimMatCond[:, sli]

    # matrix of each unit's normalized resp to p,n, p+n (9 conditions)
    unitNormHist = []
    yMax = 0
    for count, i in enumerate(stimMatCond.reshape(9)):
        dirPlot = spikeHists[unitCount, int(i), :] * 1000 / stimIndexCount[int(i)]
        smoothPlot = gaussian_filter1d(dirPlot, 5)
        if max(smoothPlot) > yMax:
            yMax = max(smoothPlot)
        unitNormHist.append(dirPlot)

    unitNormHist = np.array(unitNormHist) / yMax

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
ax.set_ylim([0, 70])
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
ax2.set_ylim([0, 70])
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
ax3.set_ylim([0, 70])

plt.show()
