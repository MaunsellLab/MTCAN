from pymatreader import read_mat
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# load .mat file
allTrials = read_mat('dotsAndST.mat')
trials = allTrials['dotsAndST']

# change trials to be list of dictionaries
nTrials = len(trials['spikeTimes'])
trials = [{k: v[i] for k, v in trials.items()} for i in range(nTrials)]

# filter for trials with EOT == 0
goodTrials = []
for i in range(len(trials)):
    if trials[i]['eotCode'] == 0:
        goodTrials.append(i)

bin1 = []
bin1dirs = []
bin1mag = []
bin2 = []
framePreSpike2 = []
framePreSpike3 = []
framePreSpike4 = []
framePreSpike5 = []

for i in goodTrials:

    xDots = trials[i]['dotXIn']
    yDots = trials[i]['dotYIn']
    st = np.where(trials[i]['spikeTimes'] == 1)[0]
    cohStepTime = trials[i]['stimOn'] + (trials[i]['frameCoh'] * (1000 / 75))
    preStepST = st[np.where((st > (trials[i]['stimOn']+100)) & (st < cohStepTime))[0]]
    for j in preStepST:
        framesBeforeSpike = np.where(trials[i]['frameTimesIn'].astype(int) - j <= 0)[0]
        frame1 = framesBeforeSpike[-2]
        # compute motion vector for 1 frame before spike
        xF0 = xDots[frame1-1, :].astype(int) * 34 / 1600
        xF1 = xDots[frame1, :].astype(int) * 34 / 1600
        yF0 = yDots[frame1-1, :].astype(int) * 26 / 1200
        yF1 = yDots[frame1, :].astype(int) * 26 / 1200
        xDif = xF1[:, None] - xF0[None, :]
        yDif = yF1[:, None] - yF0[None, :]
        magnitude = np.sqrt(xDif ** 2 + yDif ** 2).flatten() / 0.0266
        dirs = np.degrees(np.arctan2(yDif, xDif).flatten())
        dirBins = np.linspace(-180, 180, 13)
        magBins = np.linspace(0, 23.5, 9)
        hist2d = np.histogram2d(magnitude, dirs, bins=[magBins, dirBins])
        bin1.append(hist2d)

# calculate Average histogram
bin1 = np.array(bin1)
meanHist = np.mean(bin1[:, 0])
A, R = np.meshgrid(dirBins, magBins)

# plot
fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
pc = ax.pcolormesh(A, R, meanHist, cmap="magma_r")
fig.colorbar(pc)
plt.show(block=False)
# plt.pause(10)
# plt.close()


'''
def computeMVGrid(xD, yD, startFrame, xBounds, yBounds):
    """
    Given x and y coordinates for all dots in a trial,
    the function will compute all the motion vectors for one frame update that occur in a subgrid
    Inputs: xDots, yDots (2D array): rows=x/y position for each frame, columns = unique dot
    startFrame: the index of the frame after the mv update.
    xBounds, the x limits for the grid
    yBounds, the y limits for the grid
    """

    xF0 = xD[startFrame - 1, :].astype(float)
    xF1 = xD[startFrame, :].astype(float)
    yF0 = yD[startFrame - 1, :].astype(float)
    yF1 = yD[startFrame, :].astype(float)

    # mv mid points
    xMidPoint = ((xF1[:, None] + xF0[None, :]) / 2)
    yMidPoint = ((yF1[:, None] + yF0[None, :]) / 2)
    midPointInGrid = np.where((xMidPoint >= xBounds[0]) & (xMidPoint < xBounds[1]) &
                              (yMidPoint >= yBounds[0]) & (yMidPoint < yBounds[1]))

    # get dots in frame0 that start in the subgrid
    startInGrid = np.where((xF0 >= xBounds[0]) & (xF0 < xBounds[1]) &
                           (yF0 >= yBounds[0]) & (yF0 < yBounds[1]))[0]

    # get dots in frame1 that end in the subgrid
    endsInGrid = np.where((xF1 >= xBounds[0]) & (xF1 < xBounds[1]) &
                           (yF1 >= yBounds[0]) & (yF1 < yBounds[1]))[0]


    # create indices for mvs
    gridStartIndex = np.zeros((len(xF0), len(xF0)))
    gridStartIndex[:, startInGrid] = 1
    gridEndIndex = np.zeros((len(xF0), len(xF0)))
    gridEndIndex[endsInGrid, :] = 1
    midPointCentIndex = np.zeros((len(xF0), len(xF0)))
    midPointCentIndex[midPointInGrid] = 1

    validMV = ((gridStartIndex.astype('?') | gridEndIndex.astype('?'))
               & midPointCentIndex.astype('?')).flatten()
    validMVIndex = np.where(validMV == True)[0]

    # calculate MVs
    xDiff = (xF1[:, None] - xF0[None, :])
    yDiff = (yF1[:, None] - yF0[None, :])

    magnitude = np.sqrt((xDiff ** 2) + (yDiff ** 2)).flatten() / (2 / 75)
    dirs = np.arctan2(yDiff, xDiff).flatten()
    filteredMag = magnitude[validMVIndex]
    filteredDirs = dirs[validMVIndex]

    return filteredMag, filteredDirs


def computeMVCond(xD, yD, startFrame, orig, rad, subPatch=True):
    """
    Given x and y coordinates for all dots in a trial,
    the function will compute all the motion vectors for one frame update
    Inputs: xDots, yDots (2D array): rows=x/y position for each frame, columns = unique dot
    startFrame: the index of the frame after the mv update.
    orig: coordinates of the patch origin
    rad: radius of subpatch
    subPatch: true if center, false if annulus
    """

    xF0 = xD[startFrame - 1, :].astype(float)
    xF1 = xD[startFrame, :].astype(float)
    yF0 = yD[startFrame - 1, :].astype(float)
    yF1 = yD[startFrame, :].astype(float)

    # mv mid points
    xDifMid = (xF1[:, None] - xF0[None, :]) / 2
    yDifMid = (yF1[:, None] - yF0[None, :]) / 2
    xMidPoint = xF0[None, :] + xDifMid
    yMidPoint = yF0[None, :] + yDifMid
    midPointRad = np.sqrt((xMidPoint - orig[0]) ** 2 + (yMidPoint - orig[1]) ** 2)

    if subPatch == True:
        # get midpoints that are in the center
        midPointCent = np.where(midPointRad <= rad)

        # get central dots in frame0: dots that start in the center
        rad0 = np.sqrt((xF0 - orig[0]) ** 2 + (yF0 - orig[1]) ** 2)
        validDots0 = np.where(rad0 <= rad)[0]

        # get central dots in frame1: dots that end in the center
        rad1 = np.sqrt((xF1 - orig[0]) ** 2 + (yF1 - orig[1]) ** 2)
        validDots1 = np.where(rad1 <= rad)[0]

        # create indices for mvs
        centStartIndex = np.zeros((len(xF0), len(xF0)))
        centStartIndex[:, validDots0] = 1
        centEndIndex = np.zeros((len(xF0), len(xF0)))
        centEndIndex[validDots1, :] = 1
        midPointCentIndex = np.zeros((len(xF0), len(xF0)))
        midPointCentIndex[midPointCent] = 1

        validMV = ((centStartIndex.astype('?') | centEndIndex.astype('?'))
                   & midPointCentIndex.astype('?')).flatten()
        validMVIndex = np.where(validMV == True)[0]

        # calculate MVs
        xDiff = (xF1[:, None] - xF0[None, :])
        yDiff = (yF1[:, None] - yF0[None, :])

        magnitude = np.sqrt((xDiff ** 2) + (yDiff ** 2)).flatten() / (2 / 75)
        dirs = np.arctan2(yDiff, xDiff).flatten()
        filteredMag = magnitude[validMVIndex]
        filteredDirs = dirs[validMVIndex]

        return filteredMag, filteredDirs

    else:
        # get midpoints that are in the annulus
        midPointCent = np.where(midPointRad > rad)

        # get dots in frame0 that start in the annulus
        rad0 = np.sqrt((xF0 - orig[0]) ** 2 + (yF0 - orig[1]) ** 2)
        validDots0 = np.where(rad0 > rad)[0]

        # get dots in frame1 that end in the annulus
        rad1 = np.sqrt((xF1 - orig[0]) ** 2 + (yF1 - orig[1]) ** 2)
        validDots1 = np.where(rad1 > rad)[0]

        # create indices for mvs
        centStartIndex = np.zeros((len(xF0), len(xF0)))
        centStartIndex[:, validDots0] = 1
        centEndIndex = np.zeros((len(xF0), len(xF0)))
        centEndIndex[validDots1, :] = 1
        midPointCentIndex = np.zeros((len(xF0), len(xF0)))
        midPointCentIndex[midPointCent] = 1

        validMV = ((centStartIndex.astype('?') | centEndIndex.astype('?'))
                   & midPointCentIndex.astype('?')).flatten()
        validMVIndex = np.where(validMV == True)[0]

        # calculate MVs
        xDiff = (xF1[:, None] - xF0[None, :])
        yDiff = (yF1[:, None] - yF0[None, :])

        magnitude = np.sqrt((xDiff ** 2) + (yDiff ** 2)).flatten() / (2 / 75)
        dirs = np.arctan2(yDiff, xDiff).flatten()
        filteredMag = magnitude[validMVIndex]
        filteredDirs = dirs[validMVIndex]

        return filteredMag, filteredDirs

'''
########################################################################################################################

###
# 1. Compute whether each dot falls inside center for f0, f1

# 2. Generate indexing matrix 1: for every pair of dots (n x n),
# where each row/column is 1 if the start or end was in center and 0
# if neither

# 3. Generate indexing matrix 2: for every motion vector b/w every pair of dots
# (n x n), does the midpoint fall inside the center or not

# 4. Element-wise multiply the indexing matrices to generate a combination
# (1 where both (a) start or end falls in the center and (b) midpt of motion vector
# in the center, 0 otherwise)

#
framesBeforeSpike = np.where(trials[i]['frameTimesIn'].astype(int) - j <= 0)[0]
frame1 = framesBeforeSpike[-frameBeforeSpike]


xF0 = xDots[frame1 - 1, :].astype(int)
yF0 = yDots[frame1 - 1, :].astype(int)
xF1 = xDots[frame1, :].astype(int)
yF1 = yDots[frame1, :].astype(int)

# mv mid points
xDifMid = (xF1[:, None] - xF0[None, :]) / 2
yDifMid = (yF1[:, None] - yF0[None, :]) / 2
xMidPoint = xF0[None, :] + xDifMid
yMidPoint = yF0[None, :] + yDifMid
midPointRad = np.sqrt((xMidPoint - patchOrig[0]) ** 2 + (yMidPoint - patchOrig[1]) ** 2)

# get midpoints that are in the center
midPointCent = np.where(midPointRad <= spRad)

# get central dots in frame0: dots that start in the center
rad0 = np.sqrt((xDots[frame1 - 1] - patchOrig[0]) ** 2 + (yDots[frame1 - 1] - patchOrig[1]) ** 2)
validDots0 = np.where(rad0 <= spRad)[0]

# get central dots in frame1: dots that end in the center
rad1 = np.sqrt((xDots[frame1] - patchOrig[0]) ** 2 + (yDots[frame1] - patchOrig[1]) ** 2)
validDots1 = np.where(rad0 <= spRad)[0]

# create indices for mvs
centStartIndex = np.zeros((len(xF0), len(xF0)))
centStartIndex[:, validDots0] = 1
centEndIndex = np.zeros((len(xF0), len(xF0)))
centEndIndex[:, validDots1] = 1
midPointCentIndex = np.zeros((len(xF0), len(xF0)))
midPointCentIndex[midPointCent] = 1

validMV = (centStartIndex.astype('?') | centEndIndex.astype('?')
           & midPointCentIndex.astype('?')).flatten()
validMVIndex = np.where(validMV == True)[0]

# calculate MVs
xDif = (xF1[:, None] - xF0[None, :]) * 34 / 1600
yDif = (yF1[:, None] - yF0[None, :]) * 26 / 1200

magnitude = np.sqrt((xDif ** 2) + (yDif ** 2)).flatten() / (2 / 75)
dirs = np.arctan2(yDif, xDif).flatten()
filteredMag = magnitude[validMVIndex]
filteredDirs = dirs[validMVIndex]

dirBins = np.linspace(-np.pi, np.pi, numDirBins + 1)
magBins = np.linspace(0, 23.5, numSpeedBins + 1)
hist2d, _, _ = np.histogram2d(filteredMag, filteredDirs, bins=[magBins, dirBins])
bin1.append(hist2d)

########################################################################################################################
'''


# motion vectors
xDots = trials[0]['dotXIn']
yDots = trials[0]['dotYIn']
for i in range(len(xDots)-1):
    # convert to degrees
    xF0 = xDots[i, :].astype(int) * 34/1600
    xF1 = xDots[i+1, :].astype(int) * 34/1600
    yF0 = yDots[i, :].astype(int) * 26/1200
    yF1 = yDots[i+1, :].astype(int) * 26/1200

    xDif = xF1[:, None] - xF0[None, :]
    yDif = yF1[:, None] - yF0[None, :]
    magnitude = np.sqrt(xDif**2 + yDif**2).flatten() / 0.0266
    dirs = np.arctan2(yDif, xDif).flatten()

    dirBins = np.linspace(-np.pi, np.pi, 13)
    magBins = np.linspace(0, 23.5, 9)

    # calculate histogram
    hist, _, _ = np.histogram2d(dirs, magnitude, bins=[dirBins, magBins])
    A, R = np.meshgrid(dirBins, magBins)

    fig = plt.figure()
    ax = Axes3D(fig)
    plt.subplot(projection="polar")
    plt.pcolormesh(A, R, hist.T, cmap='magma_r')
    plt.plot(A, R, ls='none', color='k')
    plt.grid()
    plt.colorbar()
    plt.show(block=False)
    plt.pause(0.7)
    plt.close()

    # # plot
    # fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
    # pc = ax.pcolormesh(A, R, hist.T, cmap="magma_r")
    # fig.colorbar(pc)
    # plt.show(block=False)
    # plt.pause(.7)
    # plt.close()


#fake data
xDots = np.array([[0, 0], [1, 1]])
yDots = np.array([[0, 1], [0, 1]])

for i in range(len(xDots)-1):
    # convert to degrees
    xF0 = xDots[i, :].astype(int) * 34/1600
    xF1 = xDots[i+1, :].astype(int) * 34/1600
    yF0 = yDots[i, :].astype(int) * 26/1200
    yF1 = yDots[i+1, :].astype(int) * 26/1200

    xDif = xF1[:, None] - xF0[None, :]
    yDif = yF1[:, None] - yF0[None, :]
    magnitude = np.sqrt(xDif**2 + yDif**2).flatten()
    dirs = np.degrees(np.arctan2(yDif, xDif).flatten())

'''
% matplotlib
notebook
# compute STA CORRECTED CENTRAL DOTS
bin1 = []
spikeCount = 0
dotNumber = []
annXPos = []
annYPos = []
cenXPos = []
cenYPos = []

for i in goodTrials:

    xDots = trials[i]['dotXIn']
    yDots = trials[i]['dotYIn']
    st = np.where(trials[i]['spikeTimes'] == 1)[0] + 1
    cohStepTime = trials[i]['stimOn'] + (trials[i]['frameCoh'] * (1000 / 75))
    preStepST = st[np.where((st > (trials[i]['stimOn'] + spBuffer)) & (st < (cohStepTime - spBuffer)))[0]]

    # convert frame times to frame rate + stim on time
    frameUpdateRate = 2000 / 75
    numFrames = len(xDots)
    frameTimesInCorrected = np.array([ab * frameUpdateRate for ab in range(numFrames)])
    frameTimesInCorrected = frameTimesInCorrected + trials[i]['stimOn']

    for j in preStepST:
        framesBeforeSpike = np.where(trials[i]['frameTimesIn'].astype(int) - j <= 0)[0]
        #         framesBeforeSpike = np.where(frameTimesInCorrected - j <= 0)[0]
        frame1 = framesBeforeSpike[-frameBeforeSpike]

        # get central dots in frame0: dots that start in the center and end anywhere as long 
        # as their midway point is in the user defined smaller patch
        rad0 = np.sqrt((xDots[frame1 - 1] - patchOrig[0]) ** 2 + (yDots[frame1 - 1] - patchOrig[1]) ** 2)
        validDots0 = np.where(rad0 <= spRad)[0]

        xF0 = xDots[frame1 - 1, validDots0].astype(int)
        yF0 = yDots[frame1 - 1, validDots0].astype(int)
        xF1 = xDots[frame1, :].astype(int)
        yF1 = yDots[frame1, :].astype(int)
        xDifMid = (xF1[:, None] - xF0[None, :]) / 2
        yDifMid = (yF1[:, None] - yF0[None, :]) / 2

        # midpoints
        xMidPoint = xF0[None, :] + xDifMid
        yMidPoint = yF0[None, :] + yDifMid

        ###
        # 1. Compute whether each dot falls inside center for f0, f1

        # 2. Generate indexing matrix 1: for every pair of dots (n x n), 
        # where each row/column is 1 if the start or end was in center and 0
        # if neither

        # 3. Generate indexing matrix 2: for every motion vector b/w every pair of dots
        # (n x n), does the midpoint fall inside the center or not

        # 4. Element-wise multiply the indexing matrices to generate a combination 
        # (1 where both (a) start or end falls in the center and (b) midpt of motion vector
        # in the center, 0 otherwise)

        # 

        midPointRad = np.sqrt((xMidPoint - patchOrig[0]) ** 2 + (yMidPoint - patchOrig[1]) ** 2)

        validDotsCentStart = np.where(midPointRad <= spRad)
        vdcs = [[validDotsCentStart[0][ab], validDotsCentStart[1][ab]] for ab in range(len(validDotsCentStart[0]))]

        xDif = (xF1[:, None] - xF0[None, :])[validDotsCentStart] * 34 / 1600
        yDif = (yF1[:, None] - yF0[None, :])[validDotsCentStart] * 26 / 1200

        #         magnitude = np.sqrt((xDif ** 2) + (yDif ** 2)).flatten() / (2/75)
        #         dirs = np.arctan2(yDif, xDif).flatten()
        #         dirBins = np.linspace(-np.pi, np.pi, numDirBins+1)
        #         magBins = np.linspace(0, 23.5, numSpeedBins+1)
        #         hist2d, _, _ = np.histogram2d(magnitude, dirs, bins=[magBins, dirBins])
        #         bin1.append(hist2d)

        # get dots in frame1 that end somewhere in the center, doesn't matter where it starts as long 
        # as their midway point is in the user defined smaller patch
        rad1 = np.sqrt((xDots[frame1] - patchOrig[0]) ** 2 + (yDots[frame1] - patchOrig[1]) ** 2)
        validDots1 = np.where(rad1 <= spRad)[0]

        xF0 = xDots[frame1 - 1, :].astype(int)
        yF0 = yDots[frame1 - 1, :].astype(int)
        xF1 = xDots[frame1, validDots1].astype(int)
        yF1 = yDots[frame1, validDots1].astype(int)
        xDifMid = (xF1[:, None] - xF0[None, :]) / 2
        yDifMid = (yF1[:, None] - yF0[None, :]) / 2

        # midpoints
        xMidPoint = xF0[None, :] + xDifMid
        yMidPoint = yF0[None, :] + yDifMid

        midPointRad = np.sqrt((xMidPoint - patchOrig[0]) ** 2 + (yMidPoint - patchOrig[1]) ** 2)

        validDotsCentEnd = np.where(midPointRad <= spRad)
        vdce = [[validDotsCentEnd[0][ab], validDotsCentEnd[1][ab]] for ab in range(len(validDotsCentEnd[0]))]

        xDif = (xF1[:, None] - xF0[None, :])[validDotsCentEnd] * 34 / 1600
        yDif = (yF1[:, None] - yF0[None, :])[validDotsCentEnd] * 26 / 1200

        magnitude = np.sqrt((xDif ** 2) + (yDif ** 2)).flatten() / (2 / 75)
        dirs = np.arctan2(yDif, xDif).flatten()
        dirBins = np.linspace(-np.pi, np.pi, numDirBins + 1)
        magBins = np.linspace(0, 23.5, numSpeedBins + 1)
        hist2d, _, _ = np.histogram2d(magnitude, dirs, bins=[magBins, dirBins])
        bin1.append(hist2d)

        spikeCount += 1

#         # get x,y pos of middle dots and annulus dots for sanity check
#         for k in a:
#             cenXPos.append(xDots[frame1, k])
#             cenYPos.append(yDots[frame1, k])
#             cenXPos.append(xDots[frame1-1, k])
#             cenYPos.append(yDots[frame1-1, k])

#         # to get the dots for this frame update that lay in the annulus
#         rad1 = np.sqrt((xDots[frame1] - patchOrig[0])**2 + (yDots[frame1] - patchOrig[1])**2)
#         rad0 = np.sqrt((xDots[frame1-1] - patchOrig[0])**2 + (yDots[frame1-1] - patchOrig[1])**2)
#         validDots1 = np.where(rad1 > spRad)[0]
#         validDots0 = np.where(rad0 > spRad)[0]
#         b = np.intersect1d(validDots1, validDots0)

#         for c in b:
#             annXPos.append(xDots[frame1, c])
#             annYPos.append(yDots[frame1, c])
#             annXPos.append(xDots[frame1-1, c])
#             annYPos.append(yDots[frame1-1, c])


# plot STA
# calculate Average histogram
bin1 = np.array(bin1)
# meanHist = np.sum(bin1, axis=0) / spikeCount
meanHist = np.sum(bin1, axis=0) / np.sum(bin1)
A, R = np.meshgrid(dirBins, magBins)

# plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
pc = ax.pcolormesh(A, R, meanHist, cmap='magma')
ax.grid()
ax.set_theta_direction(-1)
ax.set_theta_zero_location('N')
plt.colorbar(pc, ax=ax)
plt.show()

'''

# Condensed Reverse
% matplotlib
notebook
# compute STA CORRECTED CENTRAL DOTS
bin1 = []
spikeCount = 0
dotNumber = []
annXPos = []
annYPos = []
cenXPos = []
cenYPos = []

for i in goodTrials:

    xDots = (trials[i]['dotXIn']).astype(int) * 34 / 1600
    yDots = (trials[i]['dotYIn']).astype(int) * 26 / 1200
    st = np.where(trials[i]['spikeTimes'] == 1)[0] + 1
    cohStepTime = trials[i]['stimOn'] + (trials[i]['frameCoh'] * (1000 / 75))
    preStepST = st[np.where((st > (trials[i]['stimOn'] + spBuffer)) & (st < (cohStepTime - spBuffer)))[0]]

    # convert frame times to frame rate + stim on time
    frameUpdateRate = 2000 / 75
    numFrames = len(xDots)
    frameTimesInCorrected = np.array([ab * frameUpdateRate for ab in range(numFrames)])
    frameTimesInCorrected = frameTimesInCorrected + trials[i]['stimOn']

    for j in preStepST:
        framesBeforeSpike = np.where(trials[i]['frameTimesIn'].astype(int) - j <= 0)[0]
        #         framesBeforeSpike = np.where(frameTimesInCorrected - j <= 0)[0]
        frame1 = framesBeforeSpike[-frameBeforeSpike]

        #         xF0 = xDots[frame1 - 1, :].astype(float)
        #         yF0 = yDots[frame1 - 1, :].astype(float)
        #         xF1 = xDots[frame1, :].astype(float)
        #         yF1 = yDots[frame1, :].astype(float)

        #         # mv mid points
        #         xDifMid = (xF1[:, None] - xF0[None, :]) / 2
        #         yDifMid = (yF1[:, None] - yF0[None, :]) / 2
        #         xMidPoint = xF0[None, :] + xDifMid
        #         yMidPoint = yF0[None, :] + yDifMid
        #         midPointRad = np.sqrt((xMidPoint - patchOrig[0]) ** 2 + (yMidPoint - patchOrig[1]) ** 2)

        #         # get midpoints that are in the center
        #         midPointCent = np.where(midPointRad <= spRad)

        #         # get central dots in frame0: dots that start in the center
        #         rad0 = np.sqrt((xDots[frame1 - 1] - patchOrig[0]) ** 2 + (yDots[frame1 - 1] - patchOrig[1]) ** 2)
        #         validDots0 = np.where(rad0 <= spRad)[0]

        #         # get central dots in frame1: dots that end in the center
        #         rad1 = np.sqrt((xDots[frame1] - patchOrig[0]) ** 2 + (yDots[frame1] - patchOrig[1]) ** 2)
        #         validDots1 = np.where(rad1 <= spRad)[0]

        #         # create indices for mvs
        #         centStartIndex = np.zeros((len(xF0), len(xF0)))
        #         centStartIndex[:, validDots0] = 1
        #         centEndIndex = np.zeros((len(xF0), len(xF0)))
        # #         centEndIndex[:, validDots1] = 1
        #         centEndIndex[validDots1, :] = 1
        #         midPointCentIndex = np.zeros((len(xF0), len(xF0)))
        #         midPointCentIndex[midPointCent] = 1

        #         validMV = ((centStartIndex.astype('?') | centEndIndex.astype('?'))
        #                    & midPointCentIndex.astype('?')).flatten()
        #         validMVIndex = np.where(validMV == True)[0]
        # #         validMV1 = centStartIndex.flatten().astype('?') & midPointCentIndex.flatten().astype('?')
        # #         indx1 = np.where(validMV1 == True)[0]
        # #         validMV2 = centEndIndex.flatten().astype('?') & midPointCentIndex.flatten().astype('?')
        # #         indx2 = np.where(validMV2 == True)[0]
        # #         fullIndx = np.unique(np.concatenate((indx1, indx2), axis=0))

        #         # calculate MVs
        #         xDif = (xF1[:, None] - xF0[None, :]) # * 34 / 1600
        #         yDif = (yF1[:, None] - yF0[None, :]) # * 26 / 1200

        #         magnitude = np.sqrt((xDif ** 2) + (yDif ** 2)).flatten() / (2 / 75)
        #         dirs = np.arctan2(yDif, xDif).flatten()
        #         filteredMag = magnitude[validMVIndex]
        #         filteredDirs = dirs[validMVIndex]

        filteredMag, filteredDirs = computeMVCond(xDots, yDots, frame1, patchOrig, spRad, subPatch=True)

        dirBins = np.linspace(-np.pi, np.pi, numDirBins + 1)
        magBins = np.linspace(0, 23.5, numSpeedBins + 1)
        hist2d, _, _ = np.histogram2d(filteredMag, filteredDirs, bins=[magBins, dirBins])
        bin1.append(hist2d)

        spikeCount += 1

#         # get x,y pos of middle dots and annulus dots for sanity check
#         for k in a:
#             cenXPos.append(xDots[frame1, k])
#             cenYPos.append(yDots[frame1, k])
#             cenXPos.append(xDots[frame1-1, k])
#             cenYPos.append(yDots[frame1-1, k])

#         # to get the dots for this frame update that lay in the annulus
#         rad1 = np.sqrt((xDots[frame1] - patchOrig[0])**2 + (yDots[frame1] - patchOrig[1])**2)
#         rad0 = np.sqrt((xDots[frame1-1] - patchOrig[0])**2 + (yDots[frame1-1] - patchOrig[1])**2)
#         validDots1 = np.where(rad1 > spRad)[0]
#         validDots0 = np.where(rad0 > spRad)[0]
#         b = np.intersect1d(validDots1, validDots0)

#         for c in b:
#             annXPos.append(xDots[frame1, c])
#             annYPos.append(yDots[frame1, c])
#             annXPos.append(xDots[frame1-1, c])
#             annYPos.append(yDots[frame1-1, c])


# plot STA
# calculate Average histogram
bin1 = np.array(bin1)
# meanHist = np.sum(bin1, axis=0) / spikeCount
meanHist = np.sum(bin1, axis=0) / np.sum(bin1)
A, R = np.meshgrid(dirBins, magBins)

# plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
pc = ax.pcolormesh(A, R, meanHist, cmap='magma')
ax.grid()
ax.set_theta_direction(-1)
ax.set_theta_zero_location('N')
plt.colorbar(pc, ax=ax)
plt.show()

# Condensed Forward

% matplotlib
notebook
forwardSTA = []
spikeCount = 0

for i in goodTrials:

    xDots = (trials[i]['dotXIn']).astype(int) * 34 / 1600
    yDots = (trials[i]['dotYIn']).astype(int) * 26 / 1200
    st = np.where(trials[i]['spikeTimes'] == 1)[0] + 1
    cohStepTime = trials[i]['stimOn'] + (trials[i]['frameCoh'] * (1000 / 75))
    preStepST = st[np.where((st > (trials[i]['stimOn'] + spBuffer)) & (st < (cohStepTime - spBuffer)))[0]]

    # convert frame times to frame rate + stim on time
    frameUpdateRate = 2000 / 75
    numFrames = len(xDots)
    frameTimesInCorrected = np.array([ab * frameUpdateRate for ab in range(numFrames)])
    frameTimesInCorrected = frameTimesInCorrected + trials[i]['stimOn']

    for j in preStepST:
        framesAfterSpike = np.where(trials[i]['frameTimesIn'].astype(int) - j > 0)[0]
        #         framesAfterSpike = np.where(frameTimesInCorrected - j > 0)[0]
        frame1 = framesAfterSpike[frameBeforeSpike - 1]

        #         # get only dots in frame update that are within the circle
        #         rad1 = np.sqrt((xDots[frame1] - patchOrig[0])**2 + (yDots[frame1] - patchOrig[1])**2)
        #         rad0 = np.sqrt((xDots[frame1-1] - patchOrig[0])**2 + (yDots[frame1-1] - patchOrig[1])**2)
        #         validDots1 = np.where(rad1 <= spRad)[0]
        #         validDots0 = np.where(rad0 <= spRad)[0]
        #         a = np.intersect1d(validDots1, validDots0)

        #         # compute motion vector for 2nd update after spike
        #         magnitude, dirs = computeMV(xDots[:, a], yDots[:, a], frame1)

        #         dirBins = np.linspace(-np.pi, np.pi, numDirBins+1)
        #         magBins = np.linspace(0, 23.5, numSpeedBins+1)
        #         hist2d, _, _ = np.histogram2d(magnitude, dirs, bins=[magBins, dirBins])
        #         forwardSTA.append(hist2d)
        #         spikeCount += 1

        #         xF0 = xDots[frame1 - 1, :].astype(float)
        #         yF0 = yDots[frame1 - 1, :].astype(float)
        #         xF1 = xDots[frame1, :].astype(float)
        #         yF1 = yDots[frame1, :].astype(float)

        #         # mv mid points
        #         xDifMid = (xF1[:, None] - xF0[None, :]) / 2
        #         yDifMid = (yF1[:, None] - yF0[None, :]) / 2
        #         xMidPoint = xF0[None, :] + xDifMid
        #         yMidPoint = yF0[None, :] + yDifMid
        #         midPointRad = np.sqrt((xMidPoint - patchOrig[0]) ** 2 + (yMidPoint - patchOrig[1]) ** 2)

        #         # get midpoints that are in the center
        #         midPointCent = np.where(midPointRad <= spRad)

        #         # get central dots in frame0: dots that start in the center
        #         rad0 = np.sqrt((xDots[frame1 - 1] - patchOrig[0]) ** 2 + (yDots[frame1 - 1] - patchOrig[1]) ** 2)
        #         validDots0 = np.where(rad0 <= spRad)[0]

        #         # get central dots in frame1: dots that end in the center
        #         rad1 = np.sqrt((xDots[frame1] - patchOrig[0]) ** 2 + (yDots[frame1] - patchOrig[1]) ** 2)
        #         validDots1 = np.where(rad1 <= spRad)[0]

        #         # create indices for mvs
        #         centStartIndex = np.zeros((len(xF0), len(xF0)))
        #         centStartIndex[:, validDots0] = 1
        #         centEndIndex = np.zeros((len(xF0), len(xF0)))
        # #         centEndIndex[:, validDots1] = 1
        #         centEndIndex[validDots1, :] = 1
        #         midPointCentIndex = np.zeros((len(xF0), len(xF0)))
        #         midPointCentIndex[midPointCent] = 1

        #         validMV = ((centStartIndex.astype('?') | centEndIndex.astype('?'))
        #                    & midPointCentIndex.astype('?')).flatten()
        #         validMVIndex = np.where(validMV == True)[0]
        # #         validMV1 = centStartIndex.flatten().astype('?') & midPointCentIndex.flatten().astype('?')
        # #         indx1 = np.where(validMV1 == True)[0]
        # #         validMV2 = centEndIndex.flatten().astype('?') & midPointCentIndex.flatten().astype('?')
        # #         indx2 = np.where(validMV2 == True)[0]
        # #         fullIndx = np.unique(np.concatenate((indx1, indx2), axis=0))

        #         # calculate MVs
        #         xDif = (xF1[:, None] - xF0[None, :]) # * 34 / 1600
        #         yDif = (yF1[:, None] - yF0[None, :]) # * 26 / 1200

        #         magnitude = np.sqrt((xDif ** 2) + (yDif ** 2)).flatten() / (2 / 75)
        #         dirs = np.arctan2(yDif, xDif).flatten()
        #         filteredMag = magnitude[validMVIndex]
        #         filteredDirs = dirs[validMVIndex]

        filteredMag, filteredDirs = computeMVCond(xDots, yDots, frame1, patchOrig, spRad, subPatch=True)

        dirBins = np.linspace(-np.pi, np.pi, numDirBins + 1)
        magBins = np.linspace(0, 23.5, numSpeedBins + 1)
        hist2d, _, _ = np.histogram2d(filteredMag, filteredDirs, bins=[magBins, dirBins])
        forwardSTA.append(hist2d)

        spikeCount += 1

# plot STA
# calculate Average histogram
forwardSTA = np.array(forwardSTA)
# forwardSTAHist = np.sum(forwardSTA, axis=0) / spikeCount
forwardSTAHist = np.sum(forwardSTA, axis=0) / np.sum(forwardSTA)
A, R = np.meshgrid(dirBins, magBins)

# plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
pc = ax.pcolormesh(A, R, forwardSTAHist, cmap='magma')
ax.grid()
ax.set_theta_direction(-1)
ax.set_theta_zero_location('N')
plt.colorbar(pc, ax=ax)
plt.show()

# annulus reverse
% matplotlib
notebook
# compute STA
bin1 = []
spikeCount = 0

for i in goodTrials:

    xDots = trials[i]['dotXIn'] * 34 / 1600
    yDots = trials[i]['dotYIn'] * 26 / 1200
    st = np.where(trials[i]['spikeTimes'] == 1)[0] + 1
    cohStepTime = trials[i]['stimOn'] + (trials[i]['frameCoh'] * (1000 / 75))
    preStepST = st[np.where((st > (trials[i]['stimOn'] + spBuffer)) & (st < (cohStepTime - spBuffer)))[0]]

    # convert frame times to frame rate + stim on time
    frameUpdateRate = 2000 / 75
    numFrames = len(xDots)
    frameTimesInCorrected = np.array([ab * frameUpdateRate for ab in range(numFrames)])
    frameTimesInCorrected = frameTimesInCorrected + trials[i]['stimOn']

    for j in preStepST:
        framesBeforeSpike = np.where(trials[i]['frameTimesIn'].astype(int) - j <= 0)[0]
        #         framesBeforeSpike = np.where(frameTimesInCorrected - j <= 0)[0]
        frame1 = framesBeforeSpike[-frameBeforeSpike]

        #         # get only dots in frame update that are within the circle
        #         rad1 = np.sqrt((xDots[frame1] - patchOrig[0])**2 + (yDots[frame1] - patchOrig[1])**2)
        #         rad0 = np.sqrt((xDots[frame1-1] - patchOrig[0])**2 + (yDots[frame1-1] - patchOrig[1])**2)
        #         validDots1 = np.where(rad1 > spRad)[0]
        #         validDots0 = np.where(rad0 > spRad)[0]
        #         a = np.intersect1d(validDots1, validDots0)

        #         # to get the dots for this frame update that lay in the annulus
        #         mask = np.ones(xDots[frame1].shape, dtype=bool)
        #         mask[a] = False
        #         b = np.where(mask==True)[0]

        #         # compute motion vector for 1 frame before spike
        #         magnitude, dirs = computeMV(xDots[:, b], yDots[:, b], frame1)

        #         dirBins = np.linspace(-np.pi, np.pi, numDirBins+1)
        #         magBins = np.linspace(0, 23.5, numSpeedBins+1)
        #         hist2d, _, _ = np.histogram2d(magnitude, dirs, bins=[magBins, dirBins])
        #         bin1.append(hist2d)

        #         xF0 = xDots[frame1 - 1, :].astype(float)
        #         yF0 = yDots[frame1 - 1, :].astype(float)
        #         xF1 = xDots[frame1, :].astype(float)
        #         yF1 = yDots[frame1, :].astype(float)

        #         # mv mid points
        #         xDifMid = (xF1[:, None] - xF0[None, :]) / 2
        #         yDifMid = (yF1[:, None] - yF0[None, :]) / 2
        #         xMidPoint = xF0[None, :] + xDifMid
        #         yMidPoint = yF0[None, :] + yDifMid
        #         midPointRad = np.sqrt((xMidPoint - patchOrig[0]) ** 2 + (yMidPoint - patchOrig[1]) ** 2)

        #         # get midpoints that are in the annulus
        #         midPointCent = np.where(midPointRad > spRad)

        #         # get annuluar dots in frame0: dots that start in the annulus
        #         rad0 = np.sqrt((xDots[frame1 - 1] - patchOrig[0]) ** 2 + (yDots[frame1 - 1] - patchOrig[1]) ** 2)
        #         validDots0 = np.where(rad0 > spRad)[0]

        #         # get annuluar dots in frame1: dots that end in the annulus
        #         rad1 = np.sqrt((xDots[frame1] - patchOrig[0]) ** 2 + (yDots[frame1] - patchOrig[1]) ** 2)
        #         validDots1 = np.where(rad1 > spRad)[0]

        #         # create indices for mvs
        #         centStartIndex = np.zeros((len(xF0), len(xF0)))
        #         centStartIndex[:, validDots0] = 1
        #         centEndIndex = np.zeros((len(xF0), len(xF0)))
        # #         centEndIndex[:, validDots1] = 1
        #         centEndIndex[validDots1, :] = 1
        #         midPointCentIndex = np.zeros((len(xF0), len(xF0)))
        #         midPointCentIndex[midPointCent] = 1

        #         validMV = ((centStartIndex.astype('?') | centEndIndex.astype('?'))
        #                    & midPointCentIndex.astype('?')).flatten()
        #         validMVIndex = np.where(validMV == True)[0]

        #         # calculate MVs
        #         xDif = (xF1[:, None] - xF0[None, :]) # * 34 / 1600
        #         yDif = (yF1[:, None] - yF0[None, :]) # * 26 / 1200

        #         magnitude = np.sqrt((xDif ** 2) + (yDif ** 2)).flatten() / (2 / 75)
        #         dirs = np.arctan2(yDif, xDif).flatten()
        #         filteredMag = magnitude[validMVIndex]
        #         filteredDirs = dirs[validMVIndex]

        filteredMag, filteredDirs = computeMVCond(xDots, yDots, frame1, patchOrig, spRad, subPatch=False)

        dirBins = np.linspace(-np.pi, np.pi, numDirBins + 1)
        magBins = np.linspace(0, 23.5, numSpeedBins + 1)
        hist2d, _, _ = np.histogram2d(filteredMag, filteredDirs, bins=[magBins, dirBins])
        bin1.append(hist2d)

        spikeCount += 1

# plot STA
# calculate Average histogram
bin1 = np.array(bin1)
# meanHist = np.sum(bin1, axis=0) / spikeCount
meanHist = np.sum(bin1, axis=0) / np.sum(bin1)
A, R = np.meshgrid(dirBins, magBins)

# plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
pc = ax.pcolormesh(A, R, meanHist, cmap='magma')
ax.grid()
ax.set_theta_direction(-1)
ax.set_theta_zero_location('N')
plt.colorbar(pc, ax=ax)
plt.show()

# annulus forward
% matplotlib
notebook
forwardSTA = []
spikeCount = 0

for i in goodTrials:

    xDots = trials[i]['dotXIn'] * 34 / 1600
    yDots = trials[i]['dotYIn'] * 26 / 1200
    st = np.where(trials[i]['spikeTimes'] == 1)[0] + 1
    cohStepTime = trials[i]['stimOn'] + (trials[i]['frameCoh'] * (1000 / 75))
    preStepST = st[np.where((st > (trials[i]['stimOn'] + spBuffer)) & (st < (cohStepTime - spBuffer)))[0]]

    # convert frame times to frame rate + stim on time
    frameUpdateRate = 2000 / 75
    numFrames = len(xDots)
    frameTimesInCorrected = np.array([ab * frameUpdateRate for ab in range(numFrames)])
    frameTimesInCorrected = frameTimesInCorrected + trials[i]['stimOn']

    for j in preStepST:
        framesAfterSpike = np.where(trials[i]['frameTimesIn'].astype(int) - j > 0)[0]
        #        framesAfterSpike = np.where(frameTimesInCorrected - j > 0)[0]
        frame1 = framesAfterSpike[frameBeforeSpike - 1]

        #         # get only dots in frame update that are within the circle
        #         rad1 = np.sqrt((xDots[frame1] - patchOrig[0])**2 + (yDots[frame1] - patchOrig[1])**2)
        #         rad0 = np.sqrt((xDots[frame1-1] - patchOrig[0])**2 + (yDots[frame1-1] - patchOrig[1])**2)
        #         validDots1 = np.where(rad1 > spRad)[0]
        #         validDots0 = np.where(rad0 > spRad)[0]
        #         a = np.intersect1d(validDots1, validDots0)

        #         # to get the dots for this frame update that lay in the annulus
        #         mask = np.ones(xDots[frame1].shape, dtype=bool)
        #         mask[a] = False
        #         b = np.where(mask==True)[0]

        #         # compute motion vector for 2nd update after spike
        #         magnitude, dirs = computeMV(xDots[:, b], yDots[:, b], frame1)

        #         dirBins = np.linspace(-np.pi, np.pi, numDirBins+1)
        #         magBins = np.linspace(0, 23.5, numSpeedBins+1)
        #         hist2d, _, _ = np.histogram2d(magnitude, dirs, bins=[magBins, dirBins])
        #         forwardSTA.append(hist2d)

        #         xF0 = xDots[frame1 - 1, :].astype(float)
        #         yF0 = yDots[frame1 - 1, :].astype(float)
        #         xF1 = xDots[frame1, :].astype(float)
        #         yF1 = yDots[frame1, :].astype(float)

        #         # mv mid points
        #         xDifMid = (xF1[:, None] - xF0[None, :]) / 2
        #         yDifMid = (yF1[:, None] - yF0[None, :]) / 2
        #         xMidPoint = xF0[None, :] + xDifMid
        #         yMidPoint = yF0[None, :] + yDifMid
        #         midPointRad = np.sqrt((xMidPoint - patchOrig[0]) ** 2 + (yMidPoint - patchOrig[1]) ** 2)

        #         # get midpoints that are in the annulus
        #         midPointCent = np.where(midPointRad > spRad)

        #         # get annuluar dots in frame0: dots that start in the annulus
        #         rad0 = np.sqrt((xDots[frame1 - 1] - patchOrig[0]) ** 2 + (yDots[frame1 - 1] - patchOrig[1]) ** 2)
        #         validDots0 = np.where(rad0 > spRad)[0]

        #         # get annuluar dots in frame1: dots that end in the annulus
        #         rad1 = np.sqrt((xDots[frame1] - patchOrig[0]) ** 2 + (yDots[frame1] - patchOrig[1]) ** 2)
        #         validDots1 = np.where(rad1 > spRad)[0]

        #         # create indices for mvs
        #         centStartIndex = np.zeros((len(xF0), len(xF0)))
        #         centStartIndex[:, validDots0] = 1
        #         centEndIndex = np.zeros((len(xF0), len(xF0)))
        # #         centEndIndex[:, validDots1] = 1
        #         centEndIndex[validDots1, :] = 1
        #         midPointCentIndex = np.zeros((len(xF0), len(xF0)))
        #         midPointCentIndex[midPointCent] = 1

        #         validMV = ((centStartIndex.astype('?') | centEndIndex.astype('?'))
        #                    & midPointCentIndex.astype('?')).flatten()
        #         validMVIndex = np.where(validMV == True)[0]

        #         # calculate MVs
        #         xDif = (xF1[:, None] - xF0[None, :]) # * 34 / 1600
        #         yDif = (yF1[:, None] - yF0[None, :]) # * 26 / 1200

        #         magnitude = np.sqrt((xDif ** 2) + (yDif ** 2)).flatten() / (2 / 75)
        #         dirs = np.arctan2(yDif, xDif).flatten()
        #         filteredMag = magnitude[validMVIndex]
        #         filteredDirs = dirs[validMVIndex]

        filteredMag, filteredDirs = computeMVCond(xDots, yDots, frame1, patchOrig, spRad, subPatch=False)

        dirBins = np.linspace(-np.pi, np.pi, numDirBins + 1)
        magBins = np.linspace(0, 23.5, numSpeedBins + 1)
        hist2d, _, _ = np.histogram2d(filteredMag, filteredDirs, bins=[magBins, dirBins])
        forwardSTA.append(hist2d)

        spikeCount += 1

# plot STA
# calculate Average histogram
forwardSTA = np.array(forwardSTA)
# forwardSTAHist = np.sum(forwardSTA, axis=0) / spikeCount
forwardSTAHist = np.sum(forwardSTA, axis=0) / np.sum(forwardSTA)
A, R = np.meshgrid(dirBins, magBins)

# plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
pc = ax.pcolormesh(A, R, forwardSTAHist, cmap='magma')
ax.grid()
ax.set_theta_direction(-1)
ax.set_theta_zero_location('N')
plt.colorbar(pc, ax=ax)
plt.show()
