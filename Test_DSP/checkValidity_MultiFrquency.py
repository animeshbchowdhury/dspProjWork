
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def checkValidity_MultiFrquency(X1, Y1, f1, NFFT1):

    # Local Variables: f1, NFFT1, index1, peak2, g2, g1, peak1, i, k, index2, chk_val, mean1, chk_val1, ratio2, Y1, mean2, X1, d2, d1
    # Function calls: sort, findpeaks, length, checkValidity_MultiFrquency, find, mean
    d1 = 0.
    d2 = 0.
    k = 1.
    g2 = 0.
    #% j = 1;
    #% flag1 = 0;
    peak1 = np.sort(findpeaks((2.*X1[0:NFFT1/2.+1.])), 'descend')
    for i in np.arange(1., (length(peak1))+1):
        chk_val1[int(i)-1,:] = matdiv(peak1[0,int(i)-1], np.mean(peak1[int(i+1.)-1:length(peak1)]))
        if chk_val1[int(i)-1,:] > 2.:
            index1[int(i)-1,:] = nonzero((X1[0:NFFT1/2.+1.] == peak1[:,int(i)-1]/2.))
            d1[0,int(k)-1] = f1[index1[int(i)-1,:]]
            k = k+1.
        
        
        
    mean1 = np.mean(peak1[0:length(peak1)])
    #% mean = mean(2*X1(1:NFFT1/2+1));
    #% for i = 1:length(peak1)
    #%     ratio1(i,:) = peak1(i)/mean1;
    #%     if(ratio1(i,:) > 2)
    #%         rat1(j,:) = ratio1(i,:);
    #%         j = j+1;
    #%         flag1 = 1;
    #%     end
    #% end
    k = 1.
    #% flag2 = 0;
    peak2 = np.sort(findpeaks((2.*Y1[0:NFFT1/2.+1.])), 'descend')
    for i in np.arange(1., (length(peak2)-1.)+1):
        chk_val[int(i)-1,:] = matdiv(peak2[0,int(i)-1], np.mean(peak2[int(i+1.)-1:length(peak2)]))
        if chk_val[int(i)-1,:] > 3.:
            g2[:,int(k)-1] = chk_val[int(i)-1,:]
            index2[int(i)-1,:] = nonzero((Y1[0:NFFT1/2.+1.] == peak2[:,int(i)-1]/2.))
            d2[0,int(k)-1] = f1[index2[int(i)-1,:]]
            k = k+1.
        
        
        
    mean2 = np.mean(peak2[0:length(peak2)])
    #% mean2(1,1:NFFT1/2+1) = 2*Y1(1,1:NFFT1/2+1);
    #% mean(mean2)
    ratio2 = matdiv(peak2[0,0], mean2)
    #% if(length(d1) == 0)
    #%     d1 = 0;
    #% end
    #% if(length(d2) == 0)
    #%     d2 = 0;
    #% end
    if g2 >= 0.:
        g1 = np.mean(g2)
    else:
        g1 = np.mean(g2)
        
    
    #% for i = 1:length(peak2)
    #%     ratio2(i,:) = peak2(i)/mean2;
    #%     if(ratio2 > 2)
    #%         rat2(k,:) = ratio2(i,:);
    #%         k = k+1;
    #%         flag2 = 1;
    #%     end
    #% end
    #% if((flag1 == 0) && (flag2 == 0))
    #%     g1 = 0;
    #% 
    #% elseif((flag1 == 1) && (flag2 == 1))
    #%     for i = 1:max(length(rat1),length(rat2))
    #%         if (max(length(rat1),length(rat2)) == length(rat1))
    #%             g1(i,:) = rat1(i,:);
    #%         else
    #%             g1(i,:) = rat2(i,:);
    #%         end
    #% % g1(i,:) = rat1(i,:)/min(length(rat1),length(rat2)) + rat2(i,:)/min(length(rat1),length(rat2));
    #%     end
    #%     
    #% elseif((flag1 == 1) && (flag2 == 0))
    #%     rat2(k,:) = 0;
    #%     for i = 1:max(length(rat1),length(rat2))
    #%         if (max(length(rat1),length(rat2)) == length(rat1))
    #%             g1(i,:) = rat1(i,:);
    #%         else
    #%             g1(i,:) = rat2(i,:);
    #%         end
    #% % g1(i,:) = rat1(i,:)/min(length(rat1),length(rat2)) + rat2(i,:)/min(length(rat1),length(rat2));
    #%     end
    #% 
    #% 
    #% else
    #%     rat1(k,:) = 0;
    #%     for i = 1:max(length(rat1),length(rat2))
    #%         if (max(length(rat1),length(rat2)) == length(rat1))
    #%             g1(i,:) = rat1(i,:);
    #%         else
    #%             g1(i,:) = rat2(i,:);
    #%         end
    #% % g1(i,:) = rat1(i,:)/min(length(rat1),length(rat2)) + rat2(i,:)/min(length(rat1),length(rat2));
    #%     end
    #% end
    #% 
    #% for i = 1:max(length(rat1),length(rat2))
    #%     if (max(length(rat1),length(rat2)) == length(rat1))
    #%         g1(i,:) = rat1(i,:);
    #%     else
    #%         g1(i,:) = rat2(i,:);
    #%     end
    #% % g1(i,:) = rat1(i,:)/min(length(rat1),length(rat2)) + rat2(i,:)/min(length(rat1),length(rat2));
    #% end
    return [g1, d1, d2]