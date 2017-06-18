
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

clear(all)
clc
matcompat.warning(off)

z = 1.
z2 = 1.
strobeFrequency = 50.
prime1 = np.array(np.hstack((251., 257., 263., 269., 271., 277., 281., 283., 293., 307., 311., 313., 317., 331., 337., 347., 349., 353., 359., 367., 373., 379., 383., 389., 397., 401., 409., 419., 421., 431., 433., 439., 443., 449., 457., 461., 463., 467., 479., 487., 491., 499., 503., 509., 521., 523., 541., 547., 557., 563., 569., 571., 577., 587., 593., 599., 601.)))
noOfIteration = 44.
pAll = prime1[0,0:noOfIteration]
iteration_primeNo = 0.
nsampl = 4.
nfreq = 1.
strobeText = fopen('strobe_file1.txt', 'w')
fprintf(strobeText, '%f\n', nfreq)
fprintf(strobeText, '%f\n', nsampl)
s = serial('/dev/ttyACM0')
fopen(s)
if strobeFrequency != prime1[0,int(z)-1]:
    b2 = prime1[0,int(z)-1]-strobeFrequency
    p1_10_quo = np.floor((b2/10.))
    p21_mod = np.mod(b2, 10.)
    p1_5_quo = np.floor((p21_mod/5.))
    p22_mod = np.mod(p21_mod, 5.)
    p1_2_quo = np.floor((p22_mod/2.))
    p1_1_quo = np.mod(p22_mod, 2.)
    np.disp('SHAKE METER Status: Scanning Multiple Frequencies.....')
    for i in np.arange(1., (p1_10_quo)+1):
        fprintf(s, 'a')
        
    for i in np.arange(1., (p1_5_quo)+1):
        fprintf(s, 'b')
        
    for i in np.arange(1., (p1_2_quo)+1):
        fprintf(s, 'c')
        
    for i in np.arange(1., (p1_1_quo)+1):
        fprintf(s, 'd')
        


tic
for autono in np.arange(1., (noOfIteration)+1):
    Signal_Generation_Experiment()
    vid.Timeout = 12.
    
    set(vid, 'FramesPerTrigger', 35.)
   
    nFrames = 35.
    vid.FramesPerTrigger = nFrames

    src = getselectedsource(vid)
    
    set(src, 'FrameRate', '30')
    
    preview(vid)

    start(vid)
   
    frames = getdata(vid)
   
    np.delete(vid)
    clear(vid)
   
    nFrames = matcompat.size(frames, 4.)
    first_frame = frames[:,:,:,0]
    rect = np.array(np.hstack((910., 420., 90., 50.)))
    first_region = imcrop(first_frame, rect)

    frame_regions = matcompat.repmat(np.uint8(0.), np.array(np.hstack((matcompat.size(first_region), nFrames))))
    for count in np.arange(1., (nFrames)+1):
        frame_regions[:,:,:,int(count)-1] = imcrop(frames[:,:,:,int(count)-1], rect)
        
   
    seg_pend = false(np.array(np.hstack((matcompat.size(first_region, 1.), matcompat.size(first_region, 2.), nFrames))))
    centroids = np.zeros(nFrames, 2.)
    se_disk = strel('disk', 3.)
    for count in np.arange(1., (nFrames)+1):
        fr = frame_regions[:,:,:,int(count)-1]
       
        gfr = rgb2gray(fr)
        gfr = imcomplement(gfr)
       
        bw = im2bw(gfr, .65)
        #% threshold is determined experimentally 
        bw = imopen(bw, se_disk)
        bw = imclearborder(bw)
        seg_pend[:,:,int(count)-1] = bw
        
        
    #%% Step 5: Find the Center of the Segmented Pendulum in Each Frame
    #% You can see that the shape of the pendulum varied in different frames.
    #% This is not a serious issue because you just need its center. You will
    #% use the pendulum centers to find the length of the pendulum.
    #%%
    #% Use |regionprops| to calculate the center of the pendulum.
    pend_centers = np.zeros((nFrames-2.), 2.)
    for count in np.arange(1., (nFrames-2.)+1):
        property = regionprops(seg_pend[:,:,int(count)-1], 'Centroid')
        if length(property) == 1.:
            pend_centers[int(count)-1,:] = property.Centroid
        elif count<=1.:
            pend_centers[int(count)-1,:] = np.array(np.hstack((0., 0.)))
            
        else:
            pend_centers[int(count)-1,:] = np.array(np.hstack((0., 0.)))
            
        
        
    #%%Vibration Measured based on the velocity information
    #% [VelocityMarker] = velocityPendCenter(pend_centers);
    x = np.transpose(pend_centers[:,0])
    y = np.transpose(pend_centers[:,1])

    [d1, X1, Y1, f1, NFFT1, d2, g1] = FFT_MultiFrequency_update(x, y)
    #%%%Plot single-sided amplitude spectrum.
    plt.figure(1.)
    
    plt.stem(f1, (2.*Y1[0:NFFT1/2.+1.]))
    plt.title('Speaker Y axis frequency plot')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('|Y(f)|')
    pause(2.)
    dimension1 = matcompat.size(d1)
    for r1 in np.arange(1., (dimension1[0,1])+1):
        x_fft[int(z)-1,int(r1)-1] = d1[0,int(r1)-1]
        
    dimension2 = matcompat.size(d2)
    for r2 in np.arange(1., (dimension2[0,1])+1):
        y_fft[int(z)-1,int(r2)-1] = d2[0,int(r2)-1]
        
    #%        disp(strobeFrequency(z,:))  
    #%        disp(prime1(:,z))
    #%        disp(x_fft)
    #%        disp(y_fft)
    #%        disp(g1)
    if g1 > 2.:
        #%            disp('take this value')
    strobe[int(z2)-1,:] = prime1[:,int(z)-1]
    [size_yfft_r, size_yfft_c] = matcompat.size(y_fft)
    remainder[int(z2)-1,0:size_yfft_c] = y_fft[int(0)-1,:]
    z2 = z2+1.
    
    
    fclose(s)
    s = serial('/dev/ttyACM0')
    fopen(s)
    z = z+1.
    #%          disp('Increase the Strobe Frequency by next prime number');
    b2 = prime1[:,int(z)-1]-prime1[:,int((z-1.))-1]
    p1_10_quo = np.floor((b2/10.))
    p21_mod = np.mod(b2, 10.)
    p1_5_quo = np.floor((p21_mod/5.))
    p22_mod = np.mod(p21_mod, 5.)
    p1_2_quo = np.floor((p22_mod/2.))
    p1_1_quo = np.mod(p22_mod, 2.)
    #%                     disp('strobe frequency increasing .....')
    for i in np.arange(1., (p1_10_quo)+1):
        fprintf(s, 'a')
        
    for i in np.arange(1., (p1_5_quo)+1):
        fprintf(s, 'b')
        
    for i in np.arange(1., (p1_2_quo)+1):
        fprintf(s, 'c')
        
    for i in np.arange(1., (p1_1_quo)+1):
        fprintf(s, 'd')
        
       
fclose(s)
[frequency1, frequency2] = calculateTripleFrequency(remainder, strobe)
counting(frequency1, frequency2)
toc
fclose(s)
