
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def FFT_MultiFrequency_update(s1, s2):

    # Local Variables: f1, NFFT1, Fs, g1, s2, s1, L, t, T, Y1, X1, d2, d1
    # Function calls: abs, fft, nextpow2, linspace, FFT_MultiFrequency_update, checkValidity_MultiFrquency, mean
    Fs = 31.
    #% Sampling frequency
    T = 1./Fs
    #% Sample time
    L = 512.
    #% Length of signal
    t = np.dot(np.arange(0., (L-1.)+1), T)
    #% Time vector
    #% % Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
    #% x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
    #% y = x + 2*randn(size(t));     % Sinusoids plus noise
    #% plot(Fs*t(1:50),y(1:50))
    #% title('Signal Corrupted with Zero-Mean Random Noise')
    #% xlabel('time (milliseconds)')
    NFFT1 = 2.**nextpow2(L)
    #% Next power of 2 from length of y
    s1[0,:] = s1[0,:]-np.mean(s1[0,:])
    s2[0,:] = s2[0,:]-np.mean(s2[0,:])
    X1 = np.abs(matdiv(np.fft(s1, NFFT1), L))
    Y1 = np.abs(matdiv(np.fft(s2, NFFT1), L))
    f1 = np.dot(Fs/2., np.linspace(0., 1., (NFFT1/2.+1.)))
    #% peak1 = sort(findpeaks(2*X1(1:NFFT1/2+1)),'descend');
    #% for i = 1:length(peak1)
    #% index1(i,:) = find(X1(1:NFFT1/2+1) == peak1(:,i)/2);
    #% d1(1,i) = f1(index1(i,:));
    #% end
    #% peak2 = sort(findpeaks(2*Y1(1:NFFT1/2+1)),'descend');
    #% for i = 1:length(peak2)
    #% index2(i,:) = find(Y1(1:NFFT1/2+1) == peak2(:,i)/2);
    #% d2(1,i) = f1(index2(i,:));
    #% end
    [g1, d1, d2] = checkValidity_MultiFrquency(X1, Y1, f1, NFFT1)
    #% for i = 1:length(g1)
    #%     if(g1(i,:) > 0.5)
    #%     [peaks1] = sort(findpeaks(2*X1(1:NFFT1/2+1)),'descend');
    #%     [peaks2] = sort(findpeaks(2*Y1(1:NFFT1/2+1)),'descend');
    #% 
    #%     index1(i,:) = find(X1(1:NFFT1/2+1) == peaks1(:,i)/2);
    #%     d1(1,i) = f1(index1(i,1));
    #% 
    #%     index2(i,:) = find(Y1(1:NFFT1/2+1) == peaks2(:,i)/2);
    #%     d2(1,i) = f1(index2(i,1));
    #%     else
    #%          d1(i,:) = 0;
    #%          d2(i,:) = 0;
    #%     end
    #% end
    return [d1, X1, Y1, f1, NFFT1, d2, g1]