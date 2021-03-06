
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def calculateTripleFrequency(remainder, strobe):

    # Local Variables: a, est_fre, remainder1, remainder2, i, k, j, freq_minus, iter, t, freq_plus, fre_11, frequency1, frequency2, freq_plus_shift_negative2, fre_1, remainder_1, remainder, freq_plus_shift_negative1, strobe
    # Function calls: calculateTripleFrequency, floor, cat, length, round, size
    #% clear all
    #% clc
    #% 
    #% load('9950.mat');
    #% strobe = [61 67 71 73 79 83 89 97 101 103 107 109 113 127 131];
    #% %  strobe = [61 67 71];
    #% remainder = [2.8457 8.9 13.14 11.86;2.90 14.5312 13.017 10.35;1.6348 11.3848 9.99 7.45;3.45 13.8 10.11 7.08;10.35 9.38 11.62 13.44;10.47 12.29 13.5 14.34;10.9 0.61 12.3 1.63;2.48 10.41 7.75 1.33;1.57 10.71 0.42 9.38;3.75 14.9 11.08 6.84;12.83 7.81 10.35 11.62;10.41 11.8 1.08 12.71;6.72 14.04 10.35 3.69;7.93 14.59 10.05 10.65;12.04 10.83 9.86 13.19];
    #% % remainder = [2.8457 8.9 13.14 11.86;2.90 14.5312 13.017 10.35;1.6348 11.3848 9.99 7.45];
    a = np.array(np.hstack((1., 2., 3.)))
    remainder1 = remainder
    remainder2 = remainder1-0.6
    for i in np.arange(1., (length(strobe))+1):
        for j in np.arange(1., (length(a))+1):
            fre_1[int(i)-1,int(j)-1] = np.dot(strobe[int(i)-1], a[int(j)-1])
            
        
    t = matcompat.size(fre_1)
    iter = 1.
    for i in np.arange(1., (t[0,0])+1):
        for j in np.arange(1., (t[0,1])+1):
            fre_11 = fre_1[int(i)-1,int(j)-1]
            remainder_1 = remainder2[int(i)-1,:]
            for k in np.arange(1., (length(remainder_1))+1):
                if remainder_1[0,int(k)-1] > 0.:
                    freq_plus[int(iter)-1,int(k)-1] = fre_11+remainder_1[0,int(k)-1]
                    freq_minus[int(iter)-1,int(k)-1] = fre_11-remainder_1[0,int(k)-1]
                    #%         freq_plus_shift_positive1(iter,k) = fre_11 + (30 + remainder_1(1,k));
                    #%         freq_plus_shift_positive2(iter,k) = fre_11 - (30 + remainder_1(1,k));
                    freq_plus_shift_negative1[int(iter)-1,int(k)-1] = fre_11+30.-remainder_1[0,int(k)-1]
                    freq_plus_shift_negative2[int(iter)-1,int(k)-1] = fre_11-30.-remainder_1[0,int(k)-1]
                
                
                
            iter = iter+1.
            
        
    est_fre = cat(1., freq_plus, freq_minus, freq_plus_shift_negative1, freq_plus_shift_negative2)
    frequency1 = np.round(est_fre)
    frequency2 = np.floor(est_fre)
    #% freq_plus_shift_positive1
    #% freq_plus_shift_positive2,
    #% disp('done');
    #% t = size(fre_1);
    #% k = 1;
    #% dr = [];
    #% 
    #% est_fre = 0;
    #% 
    #% for count = 1:length(fre_1)
    #% fre_11 = fre_1(count,:);
    #% for j = 1:t(:,2)-1
    #% remainder_1 = remainder(1,j);
    #% for i = 1:length(fre_11)
    #%     
    #%     freq_plus(k,i) = fre_11(1,i) + remainder_1;
    #%     
    #%     freq_minus(k,i) = fre_11(1,i) - remainder_1;
    #%     
    #%     freq_plus_shift_positive1(k,i) = fre_11(1,i) + (30 + remainder_1);
    #%     
    #%     freq_plus_shift_positive2(k,i) = fre_11(1,i) - (30 + remainder_1);
    #%     
    #%     freq_plus_shift_negative1(k,i) = fre_11(1,i) + (30 - remainder_1);
    #%     
    #%     freq_plus_shift_negative2(k,i) = fre_11(1,i) - (30 - remainder_1);
    #%     
    #% end
    #% if(k < 4)
    #%     k = k + 1;
    #% else
    #%     k = 1;
    #%     est_fre = cat(1,freq_plus,freq_minus,freq_plus_shift_positive1,freq_plus_shift_positive2,freq_plus_shift_negative1,freq_plus_shift_negative2);
    #%     if(count == 1)
    #%         dr = est_fre;
    #%     else
    #%         dr = vertcat(dr,est_fre);
    #%     end
    #%     
    #%     
    #% end
    #% end
    #% end
    #% 
    #% frequency1 = round(dr);
    #% frequency2 = floor(dr);
    return [frequency1, frequency2]