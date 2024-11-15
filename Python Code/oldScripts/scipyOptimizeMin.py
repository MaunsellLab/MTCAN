import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt



# chery driver function



b = meanSpikeReshaped[unitCount].reshape(7, 7) * 1000 / trueStimDurMS

# fixed (independent) variables - matrix of corresponding stim Indexes
stimMat = np.zeros((7, 7))
stimMat[:6, :6] = np.arange(36).reshape(6, 6)
stimMat[6, :6] = np.arange(36, 42)
stimMat[:, 6] = np.arange(42, 49)

# Initial guess
p0 = np.concatenate((b[6, :6], b[:6, 6], [1, 1, 0.10]), axis=0)
# p0 = [10, 15, 20, 15, 10, 5, 10, 15, 20, 15, 10, 5, 0.5, 0.5, 0.01]
bnds = ((0, None), (0, None), (0, None), (0, None), (0, None), (0, None),
        (0, None), (0, None), (0, None), (0, None), (0, None), (0, None),
        (0, None), (0, None), (0, None))


# Generic Normalization (L1+L2)/(al1+al2+sig) w.o scalar
# fits L0-L6 for both locations separately loc0 and loc1 will
# have different L0-L6 values
resp = b.reshape(49)[:-1]
fixedVals = fixedValsForGenericNorm(stimMat.reshape(49)[:-1], stimIndexDict)
# Fit curve using SCIPY MINIMIZE
res = scipy.optimize.minimize(driverFunc, p0, args=(fixedVals, resp),
                              method='Powell', bounds=bnds)

y_pred = genericNormNoScalar(fixedVals, *res.x)
r2 = r2_score(resp.squeeze(), y_pred)
print(unit, r2)
