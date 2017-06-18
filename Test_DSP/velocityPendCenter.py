
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def velocityPendCenter(pend_centers):

    # Local Variables: c, i, pend_centers, r, VelocityMarker, time
    # Function calls: velocityPendCenter, size
    time = 1./30.
    [r, c] = matcompat.size(pend_centers)
    for i in np.arange(2., (r)+1):
        VelocityMarker[int(i)-1,int((c-1.))-1] = matdiv(pend_centers[int(i)-1,int((c-1.))-1]-pend_centers[int((i-1.))-1,int((c-1.))-1], time)
        VelocityMarker[int(i)-1,int(c)-1] = matdiv(pend_centers[int(i)-1,int(c)-1]-pend_centers[int((i-1.))-1,int(c)-1], time)
        
    return [VelocityMarker]