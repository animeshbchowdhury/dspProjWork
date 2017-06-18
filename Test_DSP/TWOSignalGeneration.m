clc;
clear all;

% Fs = 44100; 
% N = 40000; % sample size
% t = (1:N)*(1/Fs); 
% freq = 600; 
% 
% sound_array = zeros(N, 2); % initializes a N x 2 matrix
% sound_array(:,1) = sin(2*pi*freq*t); 
% sound_array(:,2) = sin(2*pi*2*freq*t); 
% 
% for i = 1:50
% % Sound 1 
% sound(sound_array(:,1),Fs); 
% 
% % Sound 2 
% sound(sound_array(:,2),Fs); 
% 
% % Play the first column at left channel and the second column at the right
% % channel
% sound(sound_array,Fs);
% end


Fs = 44100;
t = [0:4/Fs:4-4/Fs];
freq = 200;
f1 = sin(2*pi*freq*t);
f2 = sin(2*pi*4*freq*t);
%Sound 1
sound(f1,Fs);
%Sound 2
sound(f2,Fs);
%Play Consecutive
f12 = [f1 f2];
sound(f12,Fs);
%Play together
f_12 = [f1+f2];
sound(f_12,Fs);

