function [frequencies,sampling]=mfreq_simulate2(frequency)


nfreq=4;
nsampl=33;
sampling_option=0;
cam_fps=15;
prime1=[61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 ];
% prime1=[61 73 79 83 97];
%prime1=[11 13 17 19 27 29 31 37 43 47 53 61 67 71 73 79 83 89 93 101 107 121 133 137];
s=fopen('strobe_file.txt','w');
freq=rand(1,nfreq);
%  frequencies=100+ceil(freq*900);
 frequencies=[70 100 170 230];
 if (sampling_option==1)
 for i=1:nsampl
     stri=['give ', num2str(i),'th frequency of strobe'];
     disp(stri);
     sampling(i)=input('');
 end
 else
     sampling(1:nsampl)=prime1(1:nsampl);
 end
 sampling;
 
 for j=1:nsampl
     for i=1:nfreq
         if mod(floor(min(mod(frequencies(i),sampling(j)), sampling(j)-mod(frequencies(i),sampling(j)))/cam_fps),2)==0
            remainder(i,j)=mod(min(mod(frequencies(i),sampling(j)), sampling(j)-mod(frequencies(i),sampling(j))),cam_fps);
         else
            remainder(i,j)=15 - mod(min(mod(frequencies(i),sampling(j)), sampling(j)-mod(frequencies(i),sampling(j))),cam_fps);
         end
     end
    remainder2(:,j)=sort(remainder(:,j));    
 end
 %remainder
 %remainder2
 fprintf(s,'%f\n',nfreq);
 fprintf(s,'%f\n',nsampl);
 for j=1:nsampl
     fprintf(s,'%f ',sampling(j));
 end
 fprintf(s,'\n');
 for i=1:nfreq
     for j=1:nsampl
        fprintf(s,'%8.2f ',remainder2(i,j));
     end
     fprintf(s,'\n');
 end
 fclose(s);
 %status=0;
