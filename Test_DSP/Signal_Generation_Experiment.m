function Signal_Generation_Experiment()

% clc
% clear all
% close all
%% configure signal settings
 dur_variation=5; % duration
 freq_variation=[500]; %frequency
 amplitude = 20;               % amplitude
 phi = 2*pi*0.5;              % phase offset, e.g.: 2*pi*0.25 = 1/4 cycle

 %% configure output settings
 fs = 16000;                    % sampling rate
 T = 1/fs;                      % sampling period
 
for i=1:length(dur_variation)
    for j=1:length(freq_variation)
        t = 0:T:dur_variation(i);              % time vector
        
%         mask=fs+1:3*fs:dur_variation(i)*fs;
        
        %% create the signal
        omega_variation = 2*pi*freq_variation(j);  
        signal = cos(omega_variation*t + phi)*amplitude; 
%         for k=1:length(mask)
%             sample_mask=mask(k):1:mask(k)+2*fs;
%             signal(sample_mask)=0;
%         end
        %pause=zeros(1,length(t));
        %signal=[signal1 pause];
       %signal = cos(freq_variation(j)*t)*amplitude; 
% %          % plot the signal
  %plot(t, signal);
% %  xlabel('Time (seconds)');
% %  ylabel('Amplitude');
% %  title('Complex Signal');

 %% play the signal
 sound(signal, fs);

 %% save signal as stereo wave file
%  stereo_signal = [signal; signal]';
%  
%  name=sprintf('%dHz%dsec.wav',freq_variation(j),dur_variation(i));
%  wavwrite(stereo_signal, fs, 16, name);
 clear signal
        
    end
end


