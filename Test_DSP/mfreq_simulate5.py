
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def mfreq_simulate2(frequency):

    # Local Variables: nfreq, remainder2, i, frequencies, j, stri, sampling, nsampl, cam_fps, frequency, prime1, freq, remainder, sampling_option
    # Function calls: disp, rand, floor, fprintf, min, fclose, sort, s, num2str, mfreq_simulate2, input, mod
    nfreq = 4.
    nsampl = 33.
    sampling_option = 0.
    cam_fps = 15.
    prime1
    freq = np.random.rand(1., nfreq)
    #%  frequencies=100+ceil(freq*900);
    frequencies = np.array(np.hstack((70., 100., 170., 230.)))
    if sampling_option == 1.:
        for i in np.arange(1., (nsampl)+1):
            stri = np.array(np.hstack(('give ', num2str(i), 'th frequency of strobe')))
            np.disp(stri)
            sampling[int(i)-1] = input(\')
            
    else:
        sampling[0:nsampl] = prime1[0:nsampl]
        
    
    sampling
    for j in np.arange(1., (nsampl)+1):
        for i in np.arange(1., (nfreq)+1):
            if np.mod(np.floor(matdiv(matcompat.max(np.mod(frequencies[int(i)-1], sampling[int(j)-1]), (sampling[int(j)-1]-np.mod(frequencies[int(i)-1], sampling[int(j)-1]))), cam_fps)), 2.) == 0.:
                remainder[int(i)-1,int(j)-1] = np.mod(matcompat.max(np.mod(frequencies[int(i)-1], sampling[int(j)-1]), (sampling[int(j)-1]-np.mod(frequencies[int(i)-1], sampling[int(j)-1]))), cam_fps)
            else:
                remainder[int(i)-1,int(j)-1] = 15.-np.mod(matcompat.max(np.mod(frequencies[int(i)-1], sampling[int(j)-1]), (sampling[int(j)-1]-np.mod(frequencies[int(i)-1], sampling[int(j)-1]))), cam_fps)
                
            
            
        remainder2[:,int(j)-1] = np.sort(remainder[:,int(j)-1])
        
    #%remainder
    #%remainder2
    fprintf(s, '%f\n', nfreq)
    fprintf(s, '%f\n', nsampl)
    for j in np.arange(1., (nsampl)+1):
        fprintf(s, '%f ', sampling[int(j)-1])
        
    fprintf(s, '\n')
    for i in np.arange(1., (nfreq)+1):
        for j in np.arange(1., (nsampl)+1):
            fprintf(s, '%8.2f ', remainder2[int(i)-1,int(j)-1])
            
        fprintf(s, '\n')
        
    fclose(s)
    #%status=0;
    return [frequencies, sampling]