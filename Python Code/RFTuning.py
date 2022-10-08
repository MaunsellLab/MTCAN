'''
RFTuning generates heatmaps and histogram of each unit's RF location.
The script will save the plots for each unit as a PDF in a folder specific
to the day's dataset. 

to do:
2D gaussian fit to isolate RF size with graphic display, find the center of the rfeceptive fiel
width, and center, elipsis on the heatmap

Chery March 2022
Modified to save plots as pngs and incorporated changes to track stimCounts
- 04/26/22
Modified to incorporate frame counter, fixed bugs relating to unit indexing
within the stimDesc, added 1-d gaussian filter - 06/10/22
Flipped heatmap and corresponding PSTHs, - 06/12/22
'''

from usefulFns import *
import scipy.io as sp
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import seaborn as sns


### for testing purposes, to make unit field similar to real data
for currTrial in allTrials:
    if 'spikeData' in currTrial:
        currTrial['spikeData']['unit'] = currTrial['spikeData']['unit'].tolist()
        for i in range(0,len(currTrial['spikeData']['channel'])):
            a = str(int(currTrial['spikeData']['channel'][i]))
            b = str(int(currTrial['spikeData']['unit'][i]))
            c = a + '_' + b
            currTrial['spikeData']['unit'][i] = c
        currTrial['spikeData']['unit'] = np.array(currTrial['spikeData']['unit'])
###

## Start here
# load data
allTrials, header = loadMatFile73('Meetz', '221005', 'Meetz_221005_GRF1_Spikes.mat')

# create folder and change dir to save PDF's and np.array
if not os.path.exists('RFLoc Tuning'):
    os.makedirs('RFLoc Tuning')
os.chdir('RFLoc Tuning/')


## Tuning code ##
units = activeUnits('spikeData', allTrials)
correctTrials = correctTrialsGRF(allTrials)

frameRateHz = header['displayCalibration']['data']['frameRateHz'].tolist()
numAzi = np.int32(header['map0Settings']['data']['azimuthDeg']['n'])
minAzi = np.float32(header['map0Settings']['data']['azimuthDeg']['minValue'])
maxAzi = np.float32(header['map0Settings']['data']['azimuthDeg']['maxValue'])
numEle = np.int32(header['map0Settings']['data']['elevationDeg']['n'])
minEle = np.float32(header['map0Settings']['data']['elevationDeg']['minValue'])
maxEle = np.float32(header['map0Settings']['data']['elevationDeg']['maxValue'])
stimDurMS = np.int32(header['mapStimDurationMS']['data'])
allTuningMat = np.zeros((len(units),numEle,numAzi))
histPrePostMS = 100
sponWindowMS = 100

# assert frame consistency for frame duration
stimDurFrame = []
for corrTrial in correctTrials:
    currTrial = allTrials[corrTrial]
    stimDesc = currTrial['stimDesc']['data']
    if 'numMap0Stim' in currTrial:
        map0StimLim = int(currTrial['numMap0Stim']['data'].tolist())
        map0Count = 0
        for stim in stimDesc:
            if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                frameDiff = stim['stimOffFrame'].tolist() - stim['stimOnFrame'].tolist()
                stimDurFrame.append(frameDiff)
if len(set(stimDurFrame)) != 1:
    print('stimulus frame duration not consistent across mapping stimuli')
else: 
    trueStimDurMS = np.around(1000/frameRateHz * stimDurFrame[0])


for uCount, unit in enumerate(units):
    stimCount = np.zeros((numEle,numAzi))
    spikeCountMat = np.zeros((25,numEle, numAzi))
    spikeCountMat[:,:,:] = np.nan
    spikeHists = np.zeros((stimDurMS + 2*histPrePostMS+12, numEle, numAzi))
    sponSpikesArr = []

    for corrTrial in correctTrials:
        currTrial = allTrials[corrTrial]
        # stimCount verify 
        if 'numMap0Stim' in currTrial:
            map0StimLim = int(currTrial['numMap0Stim']['data'].tolist())
            map0Count = 0
            stimDesc = currTrial['stimDesc']['data']
            spikeData = currTrial['spikeData']
            stim1TimeS = currTrial['taskEvents']['stimulusOn'][0]['time'].tolist()
            for stim in stimDesc:
                #extract spikes from each stimulus aligned to frame counter (on/off)
                if stim['gaborIndex'] == 1 and map0Count < map0StimLim:
                    aziIndex = int(stim['azimuthIndex'])
                    eleIndex = int(stim['elevationIndex'])
                    stCount = int(stimCount[eleIndex][aziIndex])
                    stimOnTimeS = ((1000/frameRateHz * stim['stimOnFrame'].tolist())
                                   /1000) + stim1TimeS
                    stimOffTimeS = ((1000/frameRateHz * stim['stimOffFrame'].tolist())
                                   /1000) + stim1TimeS
                    stimCount[eleIndex][aziIndex] += 1
                    map0Count += 1
                    if unit in spikeData['unit']:
                        spikeData['timeStamp'] = np.around(spikeData['timeStamp'], 3)
                        unitIndex = np.where(spikeData['unit'] == unit)[0]
                        unitTimeStamps = spikeData['timeStamp'][unitIndex]
                        stimSpikes = np.where((unitTimeStamps >= stimOnTimeS) & 
                                        (unitTimeStamps <= stimOffTimeS))
                        spikeCountMat[stCount][eleIndex][aziIndex] = len(stimSpikes[0])
                        #spontaneous spike count (x ms before stim on)
                        sponSpikes = np.where((unitTimeStamps>=(stimOnTimeS-(sponWindowMS/1000))) &
                                     (unitTimeStamps <= stimOnTimeS))
                        sponSpikesArr.extend([len(sponSpikes[0])])

                        #histograms
                        stimOnPreSNEV = stimOnTimeS - histPrePostMS/1000
                        stimOffPostSNEV = stimOffTimeS + histPrePostMS/1000
                        histStimSpikes = unitTimeStamps[((unitTimeStamps >= stimOnPreSNEV)\
                                            & (unitTimeStamps <= stimOffPostSNEV))] - stimOnPreSNEV
                        histStimSpikes = np.int32(histStimSpikes*1000)
                        spikeHists[histStimSpikes, eleIndex, aziIndex] += 1


    spikeCountMean = ma.mean(ma.masked_invalid(spikeCountMat), axis = 0)
    sponSpikesArr = np.array(sponSpikesArr)
    sponSpikesMean = np.mean(sponSpikesArr)
    allTuningMat[uCount] = spikeCountMean * 1000/trueStimDurMS

    ## Plot figure ##
    fig = plt.figure()
    fig.set_size_inches(6,8)
    aziLabel = np.around(np.linspace(minAzi, maxAzi, numAzi),2)
    eleLabel = np.around(np.linspace(minEle, maxEle, numEle),2)

    date = header['date']

    text = fig.text(0.05, 0.85, f'RF tuning for unit {unit}\n{date}\n- - - - -\n\
    Stimulus Duration = {stimDurMS} ms\nNumber of Blocks = {int(stimCount[0][0])}',\
                    size=10, fontweight='bold')
    text.set_path_effects([path_effects.Normal()])

    # heatmap
    ax_row1 = []
    ax_row1.append(plt.subplot2grid((numEle+3,6), (0,3), colspan = 3, rowspan = 3)) # ax2
    # bSmooth = gaussian_filter1d(spikeCountMean, sigma=0.5)
    bSmooth = gaussian_filter(spikeCountMean,sigma=1)
    # ax_row1[0] = sns.heatmap(spikeCountMean)
    ax_row1[0] = sns.heatmap(bSmooth)
    ax_row1[0].contour(np.arange(.5, bSmooth.shape[1]), np.arange(.5, 
                       bSmooth.shape[0]), bSmooth, colors='yellow')
    ax_row1[0].set_xlabel('azimith (˚)', fontsize=8)
    ax_row1[0].set_ylabel('elevation (˚)', fontsize = 8)
    ax_row1[0].set_title('Heatmap of unit RF location', fontsize=9)
    ax_row1[0].set_yticklabels(eleLabel, fontsize=5)
    ax_row1[0].set_xticklabels(aziLabel, fontsize=5)
    ax_row1[0].invert_yaxis()
    # plt.show()

    ax_row2 = []
    for i in np.arange(numEle+2, 2, -1):
        ax = []
        for j in np.arange(0, numAzi):
            ax.append(plt.subplot2grid((numEle+3,6), (i,j)))
        ax_row2.append(np.array(ax))

    ax_row2 = np.array(ax_row2) # 6 x 6

    # PSTHs
    yMax = 0
    for i in range(numEle):
        for j in range(numAzi):
            spikeHist = spikeHists[:,i,j] * 1000/stimCount[i,j]
            gaussSmooth = gaussian_filter1d(spikeHist, 5)
            if max(gaussSmooth) > yMax:
                yMax = max(gaussSmooth)
            ax_row2[i,j].plot(gaussSmooth)
            # ax_row2[i,j].set_title(f"Ele:{eleLabel[i]},Azi:{aziLabel[j]}", fontsize=4)
            # ax_row2[i,j].set_ylim(bottom=0)
            # ax_row2[i,j].set_ylim([0, 100])
            # ax_row2[i,j].set_yticks([0,50,100])
            # ax_row2[i,j].set_yticklabels([0,50,100], fontsize=5)
            ax_row2[i,j].tick_params(axis='y', which='major', labelsize=4)
            ax_row2[i,j].set_xticks([0,histPrePostMS,histPrePostMS+stimDurMS,2*histPrePostMS+stimDurMS])
            ax_row2[i,j].set_xticklabels([-(histPrePostMS), 0, 0+stimDurMS, stimDurMS+histPrePostMS], fontsize=3)
            ax_row2[i,j].axvspan(histPrePostMS, histPrePostMS+stimDurMS, color='grey', alpha=0.2)
            ax_row2[i,j].axhline(y=sponSpikesMean*1000/sponWindowMS, linestyle='--', color='grey')
            if i == 0 and j == 0:
                ax_row2[i,j].set_xlabel('Time (ms)', fontsize=7)
                ax_row2[i,j].set_ylabel('Firing Rate (spikes/sec)', fontsize=7)
                ax_row2[i,j].yaxis.set_label_coords(-0.5,1.75)
    plt.tight_layout(pad=0.1, w_pad=-0.25, h_pad=-1.45)

    for i in range(numEle):
        for j in range(numAzi):
            ax_row2[i,j].set_ylim([0, yMax*1.1])
    
    # saves plot as pdf
    plt.savefig(f'{unit}.pdf')
    continue

plt.close('all')
np.save('unitsRFLocMat', allTuningMat)



######## 2D gauss fit
'''
'''
allTuningMat
p_init = models.Gaussian2D(amplitude=1, x_mean=0, y_mean=0, x_stddev=None, 
                           y_stddev=None, theta=None, cov_matrix=None)
yi, xi = np.indices(spikeCountMean.shape)
fit_p = fitting.LevMarLSQFitter()
p = fit_p(p_init,xi,yi,spikeCountMean)
plt.imshow(p(xi,yi))
plt.show()

RFLocMat = spikeCountMean

#flips the spikeCountMean on its horizontal axis
for i in range(len(RFLocMat)):
    a = np.flip(RFLocMat[i], axis=0)
    com = ndimage.center_of_mass(a)
    p_init = models.Gaussian2D(amplitude=1, x_mean=com[1], y_mean=com[0], x_stddev=1, 
                                y_stddev=1, theta=None, cov_matrix=None)
    yi, xi = np.indices(a.shape)
    fit_p = fitting.LevMarLSQFitter()
    p = fit_p(p_init,xi,yi,a)

    theta = p.theta[0] * 180/np.pi
    xStdDev = p.x_stddev[0]
    yStdDev = p.y_stddev[0]
    xMean = p.x_mean[0]
    yMean = p.y_mean[0]
    amp = p.amplitude[0]

    rho = np.cos(theta)
    covMat = np.array([[xStdDev**2,rho*xStdDev*yStdDev],
                    [rho*xStdDev*yStdDev,yStdDev**2]])
    meanVec = np.array([[xMean],[yMean]])

    newMat = np.zeros((numEle*40,numAzi*40))
    yi,xi = np.indices(newMat.shape)
    p.x_mean = p.x_mean[0] * 40
    p.y_mean = p.y_mean[0] * 40
    p.x_stddev = p.x_stddev[0] * 40
    p.y_stddev = p.y_stddev[0] * 40
    modelData = p(xi,yi)
    sns.heatmap(modelData, cmap = sns.color_palette('mako', as_cmap=True),)
    plt.show()




a = np.flip(spikeCountMean, axis=0)
com = ndimage.center_of_mass(a)
p_init = models.Gaussian2D(amplitude=1, x_mean=com[1], y_mean=com[1], x_stddev=1, 
                            y_stddev=1, theta=None, cov_matrix=None)
yi, xi = np.indices(a.shape)
fit_p = fitting.LevMarLSQFitter()
p = fit_p(p_init,xi,yi,a)

theta = p.theta[0] * 180/np.pi
xStdDev = p.x_stddev[0]
yStdDev = p.y_stddev[0]
xMean = p.x_mean[0]
yMean = p.y_mean[0]
amp = p.amplitude[0]

rho = np.cos(theta)
covMat = np.array([[xStdDev**2,rho*xStdDev*yStdDev],
                [rho*xStdDev*yStdDev,yStdDev**2]])
meanVec = np.array([[xMean],[yMean]])

newMat = np.zeros((numEle*40,numAzi*40))
yi,xi = np.indices(newMat.shape)
p.x_mean = p.x_mean[0] * 40
p.y_mean = p.y_mean[0] * 40
p.x_stddev = p.x_stddev[0] * 40
p.y_stddev = p.y_stddev[0] * 40
modelData = p(xi,yi)


fig, ax = plt.subplots(figsize=(6, 6))
sns.heatmap(modelData, cmap = sns.color_palette('mako', as_cmap=True), ax=ax)
ax.contour(np.arange(.5, modelData.shape[1]), np.arange(.5, 
                       modelData.shape[0]), modelData, colors='yellow')

plt.show()




# cost2 = np.cos(theta) ** 2
# sint2 = np.sin(theta) ** 2
# sin2t = np.sin(2. * theta)
# xstd2 = xStdDev ** 2
# ystd2 = yStdDev ** 2
# a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2))
# b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2))
# c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2))


'''
#moving point average
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

spikeHist = spikeHists[:,1,1] * 1000/ stimCount[1,1]
histSmooth = smooth(spikeHist,75)
plt.plot(histSmooth)
plt.show()
'''

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)