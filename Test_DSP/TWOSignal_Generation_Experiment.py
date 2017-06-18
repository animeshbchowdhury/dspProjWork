
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def Signal_Generation_Experiment():

    # Local Variables: freq_variation2, phi, fs, i, signal, j, freq_variation1, T, amplitude, omega_variation, dur_variation, t
    # Function calls: sound, Signal_Generation_Experiment, length, pi, cos
    #% clc
    #% clear all
    #% close all
    #%% configure signal settings
    dur_variation = 5.
    #% duration
    freq_variation1 = np.array(np.hstack((10.)))
    #%frequency
    freq_variation2 = np.array(np.hstack((500.)))
    amplitude = 20.
    #% amplitude
    phi = np.dot(2.*np.pi, 0.5)
    #% phase offset, e.g.: 2*pi*0.25 = 1/4 cycle
    #%% configure output settings
    fs = 16000.
    #% sampling rate
    T = 1./fs
    #% sampling period
    for i in np.arange(1., (length(dur_variation))+1):
        for j in np.arange(1., (length(freq_variation1))+1):
            t = np.arange(0., (dur_variation[int(i)-1])+(T), T)
            #% time vector
            #%         mask=fs+1:3*fs:dur_variation(i)*fs;
            #%% create the signal
            omega_variation = np.dot(2.*np.pi, freq_variation1[int(j)-1])
            signal[0,:] = np.dot(np.cos((np.dot(omega_variation, t)+phi)), amplitude)
            #%         for k=1:length(mask)
            #%             sample_mask=mask(k):1:mask(k)+2*fs;
            #%             signal(sample_mask)=0;
            #%         end
            #%pause=zeros(1,length(t));
            #%signal=[signal1 pause];
            #%signal = cos(freq_variation(j)*t)*amplitude; 
            #% %          % plot the signal
            #%plot(t, signal);
            #% %  xlabel('Time (seconds)');
            #% %  ylabel('Amplitude');
            #% %  title('Complex Signal');
            #%% play the signal
            #%  sound(signal(1,:), fs);
            #%% save signal as stereo wave file
            #%  stereo_signal = [signal; signal]';
            #%  
            #%  name=sprintf('%dHz%dsec.wav',freq_variation(j),dur_variation(i));
            #%  wavwrite(stereo_signal, fs, 16, name);
            #%  clear signal
            
        
    for i in np.arange(1., (length(dur_variation))+1):
        for j in np.arange(1., (length(freq_variation2))+1):
            t = np.arange(0., (dur_variation[int(i)-1])+(T), T)
            #% time vector
            #%         mask=fs+1:3*fs:dur_variation(i)*fs;
            #%% create the signal
            omega_variation = np.dot(2.*np.pi, freq_variation2[int(j)-1])
            signal[1,:] = np.dot(np.cos((np.dot(omega_variation, t)+phi)), amplitude)
            #%         for k=1:length(mask)
            #%             sample_mask=mask(k):1:mask(k)+2*fs;
            #%             signal(sample_mask)=0;
            #%         end
            #%pause=zeros(1,length(t));
            #%signal=[signal1 pause];
            #%signal = cos(freq_variation(j)*t)*amplitude; 
            #% %          % plot the signal
            #%plot(t, signal);
            #% %  xlabel('Time (seconds)');
            #% %  ylabel('Amplitude');
            #% %  title('Complex Signal');
            #%% play the signal
            #%  sound(signal(2,:), fs);
            #%% save signal as stereo wave file
            #%  stereo_signal = [signal; signal]';
            #%  
            #%  name=sprintf('%dHz%dsec.wav',freq_variation(j),dur_variation(i));
            #%  wavwrite(stereo_signal, fs, 16, name);
            #%  clear signal
            
        
    sound(signal[0,:], fs)
    sound(signal[1,:], fs)
    return 