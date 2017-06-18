function [freq_all,freq_all_all]=mfreq_solve12()

s=fopen('strobe_file.txt','r');
nfreq=fscanf(s,'%f',1);
nsampl=fscanf(s,'%f',1);
for j=1:nsampl
    frequencies(j)=fscanf(s,'%f',1);
end

for i=1:nfreq
    for j=1:nsampl
    sampling(i,j)=fscanf(s,'%f',1);
    end
end
fclose(s);


nfreq_orig=nfreq;
nsampl_orig=nsampl;

for i1=1:nfreq_orig
    il_all=1;
    one_frequency=-1;
    while (one_frequency <0)
           if (nfreq<8)
               nsampl=min(nfreq*3, nsampl); 
           end
           margin=nsampl_orig-nsampl;
       
         if margin > 0
             rand_start=floor(rand(1,1)*(margin-1))+1;
          else
             rand_start=1;
         end
        [one_frequency,probable_mod1]=find_frequency(nfreq,nsampl,nsampl_orig,frequencies,sampling,sampling(1:nfreq,rand_start:rand_start-1+nsampl));
    end
    num_freq=numel(one_frequency);
      
    for i_num=1:num_freq  
    if one_frequency(i_num) >0
        freq_all(i1)=one_frequency(i_num);
        freq_all_all1(il_all)=one_frequency(i_num);
        il_all=il_all+1;
    end;
    end;
    sort(freq_all)
    %pause;
    %probable_mod;
  


  
    for j=1:nsampl_orig
        count=0;
      
        for k=1:nfreq-1       
              if (count== 1) || (sampling(k,j)==probable_mod1(j))
                   count=1;
                   sampling(k,j)=sampling(k+1,j);
              end
        end       
    end
 %sampling;
 nfreq=nfreq-1;
 
     %sampling
     %pause;
 
end
%freq_all;
%nfreq;

% for j=1:nsampl_orig
%      for i=1:nfreq_orig
%     remainder(i,j)=mod(freq_all(i),frequencies(j));
%      end
%     remainder2(:,j)=sort(remainder(:,j));  
% end


 sort(freq_all);
 freq_all_all=unique(sort(freq_all_all1));
end



function [one_frequency, probable_mod1]=find_frequency(nfreq,nsampl,nsampl_orig,frequencies,sampling_all,sampling)
fps=15;
eta=floor(nsampl/nfreq);
for j=1:nsampl
    test_set(j)=sampling(1,j);
end
if (eta ==1)
    eta=2;
end
fields_to_check=zeros(eta);

% fields_to_check=[1 2];
% end;

% fields_to_check=fields_to_check-1;
while numel(unique(fields_to_check)) < eta
    fields_to_check= floor(rand(1,eta)*nsampl)+1;
end
%eta
%fields_to_check
%pause;
search_limit=1;
for k=1:eta
    search_limit=search_limit*frequencies(fields_to_check(k));
end
if (search_limit > 2000)
    search_limit=2000;
end;

for k=1:eta
    m=2;
    field(k,2)=test_set(fields_to_check(k));
    %pause
    num_fps(k)=floor(frequencies(fields_to_check(k))/2/fps);
    %frequencies(fields_to_check(k))
    %num_fps(k)
    %pause;
    for  k_fps=1:num_fps(k)+3 
        if (mod(k_fps,2) == 0 )
            tempfkm=(k_fps-1)*fps+(fps-field(k,2));
            if tempfkm < frequencies(fields_to_check(k))/2
                field(k,m)=tempfkm;
            end
            
        else
            tempfkm=(k_fps-1)*fps+field(k,2);
            if tempfkm < frequencies(fields_to_check(k))/2
                field(k,m)=tempfkm;
            end
        end
        %tempfkm
        %pause
        mm1=1;
        mm2=0;
        while(-tempfkm+mm1*frequencies(fields_to_check(k)) < search_limit)
            field(k,m+1)=-tempfkm+mm1*frequencies(fields_to_check(k));
            mm1=mm1+1;
            m=m+1;
        end
        while(tempfkm+mm2*frequencies(fields_to_check(k)) < search_limit)
            field(k,m+1)=tempfkm+mm2*frequencies(fields_to_check(k));
            mm2=mm2+1;
            m=m+1;
        end
    end
    field(k,1)=m-2;
end
%field
%search_limit
%pause;
   field_common1(1,1)=0;   
for i=1:eta-1
        cardinality=numel(intersect(field(i,2:field(i,1)+1),field(i+1,2:field(i+1,1)+1)));
    field_common1(1:cardinality)=intersect(field(i,2:field(i,1)+1),field(i+1,2:field(i+1,1)+1));
    field_common1(:);
    if (numel(field_common1)>0)
        field(i+1,1)=numel(field_common1(:));
        field(i+1,2:numel(field_common1(:))+1)=field_common1(:);
    else
        break;
    end
end
    if numel(field_common1)< 1
        field_common=0;
    else
        field_common=field_common1(:);
      
   end;
   %field_common
   %pause;
   %frequencies
   fps;
   %pause
   num_of_mod=numel(field_common);
 
   
   for i_nom=1:num_of_mod
       for i_sampln=1:nsampl_orig
        if mod(floor(min(mod(field_common(i_nom),frequencies(i_sampln)), frequencies(i_sampln)-mod(field_common(i_nom),frequencies(i_sampln)))/fps),2)==0
            probable_mod(i_nom,i_sampln)=mod(min(mod(field_common(i_nom),frequencies(i_sampln)), frequencies(i_sampln)-mod(field_common(i_nom),frequencies(i_sampln))),fps);
         else
            probable_mod(i_nom,i_sampln)=15 - mod(min(mod(field_common(i_nom),frequencies(i_sampln)), frequencies(i_sampln)-mod(field_common(i_nom),frequencies(i_sampln))),fps);
        end
       end
   end
   
   
   %probable_mod
   %pause;
   [pmr,pmc]=size(probable_mod);
 
    probable_mod1=zeros(pmc);
   %frequencies
common=zeros(num_of_mod);
%fps
field_common;
%frequencies
%probable_mod
for i_nom=1:num_of_mod
for j=1:nsampl_orig
     common(i_nom)=common(i_nom)+numel(intersect(probable_mod(i_nom,j),sampling_all(:,j)));
end
common(i_nom);
if (common(i_nom) >=  1*nsampl_orig )
    one_frequency(i_nom)=field_common(i_nom);
    field_common(i_nom);
    probable_mod1(1:pmc)=probable_mod(i_nom,1:pmc);
%     disp('The Vibration Frequency of the machine:::::::::::::::');
    %one_frequency
    %pause;
else
 %   [pmr,pmc]=size(probable_mod);
    one_frequency(i_nom)=-1;
    %disp('hi')
  %  probable_mod1=zeros(pmc);
end
end
end
