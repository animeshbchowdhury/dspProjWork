
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def calculateFrequency(d11, z):

    # Local Variables: diff, z, d11, w
    # Function calls: abs, calculateFrequency, round, mode
    #%     diff(1,:) = abs(d(z-6)-d(z-5));
    #%     diff(2,:) = abs(d(z-5)-d(z-4));
    #%     diff(3,:) = abs(d(z-4)-d(z-3));
    #%     diff(4,:) = abs(d(z-3)-d(z-2));
    diff[0,:] = np.abs(np.round((d11[int((z-3.))-1]-d11[int((z-2.))-1])))
    diff[1,:] = np.abs(np.round((d11[int((z-2.))-1]-d11[int((z-1.))-1])))
    diff[2,:] = np.abs(np.round((d11[int((z-1.))-1]-d11[int(z)-1])))
    #%     w = int16(diff(2,:));
    w = mode(diff)
    return [w]