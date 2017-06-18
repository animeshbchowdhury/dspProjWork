clear all
clc
warning off;

% frequency = input('please provide the frequency value');

% disp('PLEASE RESET THE BROAD')
%  pause(10);
%  s = serial('/dev/ttyUSB0');
z = 1;
z2 = 1;
strobeFrequency = 50; 
% prime1=[61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257];
% prime1=[61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281 283 293 307 311 313 317 331 337 347 349 353 359 367 373 379 383 389 397 401 409];
prime1=[251 257 263 269 271 277 281 283 293 307 311 313 317 331 337 347 349 353 359 367 373 379 383 389 397 401 409 419 421 431 433 439 443 449 457 461 463 467 479 487 491 499 503 509 521 523 541 547 557 563 569 571 577 587 593 599 601];
noOfIteration = 44;
pAll = prime1(1,1:noOfIteration);

iteration_primeNo = 0;
nsampl = 4;
nfreq = 1;
strobeText = fopen('strobe_file1.txt','w');
fprintf(strobeText,'%f\n',nfreq);
fprintf(strobeText,'%f\n',nsampl);


s = serial('/dev/ttyACM0');
fopen(s);
if(strobeFrequency ~= prime1(1,z))
         b2 = prime1(1,z) - strobeFrequency;
         p1_10_quo = floor(b2./10);
         p21_mod = mod(b2,10);
         p1_5_quo = floor(p21_mod./5);
         p22_mod = mod(p21_mod,5);
         p1_2_quo = floor(p22_mod./2);
         p1_1_quo = mod(p22_mod,2);
         disp('SHAKE METER Status: Scanning Multiple Frequencies.....')
         for i = 1 : p1_10_quo
             fprintf(s,'a');
         end
         for i = 1 : p1_5_quo
             fprintf(s,'b');
         end
         for i = 1 : p1_2_quo
             fprintf(s,'c');
         end
         for i = 1 : p1_1_quo
             fprintf(s,'d');
         end
end

tic

for autono = 1:noOfIteration
  Signal_Generation_Experiment();
    
%% Step 1: Acquire Images
% Load the image frames of a pendulum in motion. The frames in the MAT-file
% |pendulum.mat| were acquired using the following functions in the Image
% Acquisition Toolbox.

% Access an image acquisition device (video object).
% vidimage=videoinput('winvideo',1,'RGB24_352x288');
  vid = videoinput('linuxvideo',1,'RGB24_1920x1080');
  vid.Timeout = 12;
% Configure object to capture every fifth frame.
% vidimage.FrameGrabInterval = 5;

%Configure the no of frames to be logged
  set(vid, 'FramesPerTrigger',35);
  
% Configure the number of frames to be logged.
  nFrames=35;
  vid.FramesPerTrigger = nFrames;

% Access the device's video source.
  src=getselectedsource(vid);

%Configure the to provide 30 frames per second
  set(src, 'FrameRate', '30');

% Open a live preview window. Focus camera onto a moving pendulum.
  preview(vid);

%   pause(1);
% Initiate the acquisition.
  start(vid);

% Wait for data logging to finish before retrieving the data.
% wait(vidimage, 10);

% Extract frames from memory.
  frames = getdata(vid);

% Clean up. Delete and clear associated variables.
  delete(vid)
  clear vid

%   load MAT-file
% 
%     load SpeakerStrobeIteration10;
    

%% Step 2: Explore Sequence with IMPLAY
% Run the following command to explore the image sequence in |implay|.

%    implay(frames);

%% Step 3: Select Region where Pendulum is Swinging
% You can see that the pendulum is swinging in the upper half of each frame
% in the image series.  Create a new series of frames that contains only
% the region where the pendulum is swinging.
%
% To crop a series of frames using |imcrop|, first perform |imcrop| on one
% frame and store its output image. Then use the previous output's size to
% create a series of frame regions.  For convenience, use the |rect| that
% was loaded by |pendulum.mat| in |imcrop|.

nFrames = size(frames,4);
first_frame = frames(:,:,:,1);
rect = [910 420 90 50];
first_region = imcrop(first_frame,rect);
% imshow(first_region);
frame_regions = repmat(uint8(0), [size(first_region) nFrames]);
for count = 1:nFrames
  frame_regions(:,:,:,count) = imcrop(frames(:,:,:,count),rect);
end
%  imshow(frames(:,:,:,1))

%%
%  imshow(frame_regions(:,:,:,1));

%% Step 4: Segment the Pendulum in Each Frame
% Notice that the pendulum is much darker than the background.  You can
% segment the pendulum in each frame by converting the frame to grayscale,
% thresholding it using |im2bw|, and removing background structures using
% |imopen| and |imclearborder|.

% initialize array to contain the segmented pendulum frames.
seg_pend = false([size(first_region,1) size(first_region,2) nFrames]);
centroids = zeros(nFrames,2);
se_disk = strel('disk',3);

for count = 1:nFrames
    fr = frame_regions(:,:,:,count);
%     imshow(fr)
%     pause(0.2)
    
    gfr = rgb2gray(fr);
    gfr = imcomplement(gfr);
%     imshow(gfr)
%     pause(0.2)
    
    bw = im2bw(gfr,.65);  % threshold is determined experimentally 
    bw = imopen(bw,se_disk);
    bw = imclearborder(bw);
    seg_pend(:,:,count) = bw;
%     imshow(bw)
%     pause(0.2)
end

%% Step 5: Find the Center of the Segmented Pendulum in Each Frame
% You can see that the shape of the pendulum varied in different frames.
% This is not a serious issue because you just need its center. You will
% use the pendulum centers to find the length of the pendulum.

%%
% Use |regionprops| to calculate the center of the pendulum.
pend_centers = zeros(nFrames-2,2);
for count = 1:nFrames-2
    property = regionprops(seg_pend(:,:,count), 'Centroid');
    if(length(property) == 1)
    pend_centers(count,:) = property.Centroid;
    elseif(count <= 1)
        pend_centers(count,:) = [0 0];
    else
%         pend_centers(count,:) = pend_centers(count-1,:);
          pend_centers(count,:) = [0 0];
    end
end

%%Vibration Measured based on the velocity information
% [VelocityMarker] = velocityPendCenter(pend_centers);

x = transpose(pend_centers(:,1));
y = transpose(pend_centers(:,2));

%%
% Display pendulum centers using |plot|.

% x = transpose(pend_centers(:,1));
% y = transpose(pend_centers(:,2));
% plot(x);
% plot(y);
% figure
% plot(x,y,'m.')
% axis ij
% axis equal
% hold on;
% xlabel('x');
% ylabel('y');
% title('pendulum centers');

% disp('Computing FFT for region 1');
[d1,X1,Y1,f1,NFFT1,d2,g1] = FFT_MultiFrequency_update(x,y);
 
%%%Plot single-sided amplitude spectrum.
figure(1)
% subplot(2,1,1)
% stem(f1,2*(X1(1:NFFT1/2+1))) 
% title('Speaker X axis frequency plot')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
%                 
% subplot(2,1,2)
stem(f1,2*(Y1(1:NFFT1/2+1))) 
title('Speaker Y axis frequency plot')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

      
       pause(2);
       dimension1 = size(d1);
       for r1 = 1:dimension1(1,2)
           x_fft(z,r1) = d1(1,r1);
       end
       
       dimension2 = size(d2);
       for r2 = 1:dimension2(1,2)
           y_fft(z,r2) = d2(1,r2);
       end
       
       
       
%        disp(strobeFrequency(z,:))  
%        disp(prime1(:,z))
%        disp(x_fft)
%        disp(y_fft)
%        disp(g1)
       
       if(g1 > 2)
%            disp('take this value')
           strobe(z2,:) = prime1(:,z);
           [size_yfft_r,size_yfft_c] = size(y_fft);
           remainder(z2,1:size_yfft_c) = y_fft(end,:);
           z2 = z2 + 1;
       end
%        if(g1 >= 4)
%            
%            disp('TAKE THIS VALUE');
%            disp(prime1(:,z))
%            disp(y_fft(end,1:2))
%            iteration_primeNo = iteration_primeNo + 1;
%            sampling(1,iteration_primeNo) = prime1(:,z);
%            remainder(1,iteration_primeNo) = round(y_fft(end,1));
%                      
% %                  fprintf(strobeText,'%f ',prime1(:,z));
% %            
% %                  fprintf(strobeText,'\n');
% % %              for i=1:nfreq
% % %                  for j=1:nsampl
% %                     fprintf(strobeText,'%8.2f ',round(y_fft(end,1:1)));
% % %                  end
% %                  fprintf(strobeText,'\t');
%         
%              if(iteration_primeNo == nsampl+1)
%                  for j=2:nsampl+1
%                     fprintf(strobeText,'%f ',sampling(j));
%                  end
%                  fprintf(strobeText,'\n');
%                  for i=1:nfreq
%                      for j=2:nsampl+1
%                         fprintf(strobeText,'%8.2f ',remainder(i,j));
%                      end
%                      fprintf(strobeText,'\n');
%                  end
%                  fclose(strobeText);
%                  fclose(s);
%                  clc
%                  clf
%                  mfreq_simulate5(frequency);
%                  freq_all = mfreq_solve13();
%                  clc
%                  sprintf('The Vibration Frequency of the speaker is %d Hz', freq_all)
%                  toc
%                break;
%              end
%                      
%        end
       
%        if(z>3)
%            if(y_fft(z,1) > 0 && y_fft(z,1) < 1)
%                disp('inside decision loop')
%                a11 = ceil(y_fft(z-1,1) - y_fft(z,1))
%                pause(1)
%                a12 = round(y_fft(z-2,1) - y_fft(z-1,1))
%                pause(1)
%                a13 = round(y_fft(z-3,1) - y_fft(z-2,1))
%                if((a11 == a12) || (a12 == a13) || (a11 == a13))
%                    [w] = calculateFrequency(y_fft,z)
%                    pause(1)
%                    if(w == 1)
%                        disp('the vibrating frequecy of the system is')
%                        strobeFrequency(z,:) * w + y_fft(z,1)
%                        fclose(s);
%                        break;
%                    else
%                        disp('Increase the strobe frequecy to a certain limit')
%                        b1 = strobeFrequency(z,:) * w - 10;
%                        b2 = b1 - strobeFrequency(z,:);
%                        strobeFrequency(z,:) = strobeFrequency(z-2,:) + b2;
%                        p1_10_quo = floor(b2./10);
%                        p21_mod = mod(b2,10);
%                        p1_5_quo = floor(p21_mod./5);
%                        p22_mod = mod(p21_mod,5);
%                        p1_2_quo = floor(p22_mod./2);
%                        p1_1_quo = mod(p22_mod,2);
%                        disp('strobe frequency increasing .....')
%                        for i = 1 : p1_10_quo
%                            fprintf(s,'a');
%                        end
%                        for i = 1 : p1_5_quo
%                            fprintf(s,'b');
%                        end
%                        for i = 1 : p1_2_quo
%                            fprintf(s,'c');
%                        end
%                        for i = 1 : p1_1_quo
%                            fprintf(s,'d');
%                        end
%                    end
% 
%                end
% 
%            end
%        end

        %%%%%%%%Save frequency data into file%%%%%%%%%%%%%%
%         fName = 'strobeIteration';
%         Iteration = num2str(strobeFrequency(z,:));
%         filename = strcat(fName,Iteration);
%         save(filename, 'd1', 'd2', 'X1','Y1','f1','NFFT1');
        
         fclose(s);
         s = serial('/dev/ttyACM0');
         fopen(s);
         
         z = z+1;
%          disp('Increase the Strobe Frequency by next prime number');
         b2 = prime1(:,z) - prime1(:,z-1);
         p1_10_quo = floor(b2./10);
         p21_mod = mod(b2,10);
         p1_5_quo = floor(p21_mod./5);
         p22_mod = mod(p21_mod,5);
         p1_2_quo = floor(p22_mod./2);
         p1_1_quo = mod(p22_mod,2);
%                     disp('strobe frequency increasing .....')
                    
                       for i = 1 : p1_10_quo
                           fprintf(s,'a');
                       end
                       
                       for i = 1 : p1_5_quo
                           fprintf(s,'b');
                       end
                       
                       for i = 1 : p1_2_quo
                           fprintf(s,'c');
                       end
                       
                       for i = 1 : p1_1_quo
                           fprintf(s,'d');
                       end
%          fprintf(s, 'd');
%          strobe = fscanf(s,'%e',20);
%          strobeFrequency(z,:) = strobeFrequency(z-1,:) + 1;
%          fName = 'SpeakerStrobeIteration';
%          Iteration = num2str(z);
%          filename = strcat(fName,Iteration);
%          save(filename);
%          


end
fclose(s);
[frequency1,frequency2] =  calculateTripleFrequency(remainder,strobe);
counting(frequency1,frequency2);
%  save('frequencyValue.mat');
toc

fclose(s);
