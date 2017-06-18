
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def mfreq_solve12():

    # Local Variables: il_all, num_freq, sampling, nsampl, one_frequency, freq_all, nfreq, nsampl_orig, freq_all_all1, rand_start, i1, frequencies, probable_mod1, freq_all_all, count, nfreq_orig, i, k, j, i_num, s, margin
    # Function calls: sort, rand, unique, floor, mfreq_solve12, fclose, find_frequency, fscanf, min, numel, fopen
    s = fopen('strobe_file.txt', 'r')
    nfreq = fscanf(s, '%f', 1.)
    nsampl = fscanf(s, '%f', 1.)
    for j in np.arange(1., (nsampl)+1):
        frequencies[int(j)-1] = fscanf(s, '%f', 1.)
        
    for i in np.arange(1., (nfreq)+1):
        for j in np.arange(1., (nsampl)+1):
            sampling[int(i)-1,int(j)-1] = fscanf(s, '%f', 1.)
            
        
    fclose(s)
    nfreq_orig = nfreq
    nsampl_orig = nsampl
    for i1 in np.arange(1., (nfreq_orig)+1):
        il_all = 1.
        one_frequency = -1.
        while one_frequency<0.:
            if nfreq<8.:
                nsampl = matcompat.max((nfreq*3.), nsampl)
            
            
            margin = nsampl_orig-nsampl
            if margin > 0.:
                rand_start = np.floor(np.dot(np.random.rand(1., 1.), margin-1.))+1.
            else:
                rand_start = 1.
                
            
            [one_frequency, probable_mod1] = find_frequency(nfreq, nsampl, nsampl_orig, frequencies, sampling, sampling[0:nfreq,int(rand_start)-1:rand_start-1.+nsampl])
            
        num_freq = numel(one_frequency)
        for i_num in np.arange(1., (num_freq)+1):
            if one_frequency[int(i_num)-1] > 0.:
                freq_all[int(i1)-1] = one_frequency[int(i_num)-1]
                freq_all_all1[int(il_all)-1] = one_frequency[int(i_num)-1]
                il_all = il_all+1.
            
            
            
        np.sort(freq_all)
        #%pause;
        #%probable_mod;
        for j in np.arange(1., (nsampl_orig)+1):
            count = 0.
            for k in np.arange(1., (nfreq-1.)+1):
                if count == 1. or sampling[int(k)-1,int(j)-1] == probable_mod1[int(j)-1]:
                    count = 1.
                    sampling[int(k)-1,int(j)-1] = sampling[int((k+1.))-1,int(j)-1]
                
                
                
            
        #%sampling;
        nfreq = nfreq-1.
        #%sampling
        #%pause;
        
    #%freq_all;
    #%nfreq;
    #% for j=1:nsampl_orig
    #%      for i=1:nfreq_orig
    #%     remainder(i,j)=mod(freq_all(i),frequencies(j));
    #%      end
    #%     remainder2(:,j)=sort(remainder(:,j));  
    #% end
    np.sort(freq_all)
    freq_all_all = np.unique(np.sort(freq_all_all1))
    return [freq_all, freq_all_all]
def find_frequency(nfreq, nsampl, nsampl_orig, frequencies, sampling_all, sampling):

    # Local Variables: i_nom, pmr, mm1, sampling, nsampl, one_frequency, num_of_mod, pmc, nfreq, nsampl_orig, cardinality, test_set, i_sampln, probable_mod, tempfkm, field_common1, field, fields_to_check, fps, field_common, frequencies, probable_mod1, sampling_all, mm2, num_fps, i, k, j, m, eta, search_limit, common, k_fps
    # Function calls: size, rand, intersect, floor, min, find_frequency, zeros, numel, unique, mod
    fps = 15.
    eta = np.floor(matdiv(nsampl, nfreq))
    for j in np.arange(1., (nsampl)+1):
        test_set[int(j)-1] = sampling[0,int(j)-1]
        
    if eta == 1.:
        eta = 2.
    
    
    fields_to_check = np.zeros(eta)
    #% fields_to_check=[1 2];
    #% end;
    #% fields_to_check=fields_to_check-1;
    while numel(np.unique(fields_to_check))<eta:
        fields_to_check = np.floor(np.dot(np.random.rand(1., eta), nsampl))+1.
        
    #%eta
    #%fields_to_check
    #%pause;
    search_limit = 1.
    for k in np.arange(1., (eta)+1):
        search_limit = np.dot(search_limit, frequencies[int(fields_to_check[int(k)-1])-1])
        
    if search_limit > 2000.:
        search_limit = 2000.
    
    
    for k in np.arange(1., (eta)+1):
        m = 2.
        field[int(k)-1,1] = test_set[int(fields_to_check[int(k)-1])-1]
        #%pause
        num_fps[int(k)-1] = np.floor(matdiv(frequencies[int(fields_to_check[int(k)-1])-1]/2., fps))
        #%frequencies(fields_to_check(k))
        #%num_fps(k)
        #%pause;
        for k_fps in np.arange(1., (num_fps[int(k)-1]+3.)+1):
            if np.mod(k_fps, 2.) == 0.:
                tempfkm = np.dot(k_fps-1., fps)+fps-field[int(k)-1,1]
                if tempfkm<frequencies[int(fields_to_check[int(k)-1])-1]/2.:
                    field[int(k)-1,int(m)-1] = tempfkm
                
                
            else:
                tempfkm = np.dot(k_fps-1., fps)+field[int(k)-1,1]
                if tempfkm<frequencies[int(fields_to_check[int(k)-1])-1]/2.:
                    field[int(k)-1,int(m)-1] = tempfkm
                
                
                
            
            #%tempfkm
            #%pause
            mm1 = 1.
            mm2 = 0.
            while -tempfkm+np.dot(mm1, frequencies[int(fields_to_check[int(k)-1])-1])<search_limit:
                field[int(k)-1,int((m+1.))-1] = -tempfkm+np.dot(mm1, frequencies[int(fields_to_check[int(k)-1])-1])
                mm1 = mm1+1.
                m = m+1.
                
            while tempfkm+np.dot(mm2, frequencies[int(fields_to_check[int(k)-1])-1])<search_limit:
                field[int(k)-1,int((m+1.))-1] = tempfkm+np.dot(mm2, frequencies[int(fields_to_check[int(k)-1])-1])
                mm2 = mm2+1.
                m = m+1.
                
            
        field[int(k)-1,0] = m-2.
        
    #%field
    #%search_limit
    #%pause;
    field_common1[0,0] = 0.
    for i in np.arange(1., (eta-1.)+1):
        cardinality = numel(intersect(field[int(i)-1,1:field[int(i)-1,0]+1.], field[int((i+1.))-1,1:field[int((i+1.))-1,0]+1.]))
        field_common1[0:cardinality] = intersect(field[int(i)-1,1:field[int(i)-1,0]+1.], field[int((i+1.))-1,1:field[int((i+1.))-1,0]+1.])
        field_common1.flatten(1)
        if numel(field_common1) > 0.:
            field[int((i+1.))-1,0] = numel(field_common1.flatten(1))
            field[int((i+1.))-1,1:numel[field_common1[:]]+1.] = field_common1.flatten(1)
        else:
            break
            
        
        
    if numel(field_common1)<1.:
        field_common = 0.
    else:
        field_common = field_common1.flatten(1)
        
    
    #%field_common
    #%pause;
    #%frequencies
    fps
    #%pause
    num_of_mod = numel(field_common)
    for i_nom in np.arange(1., (num_of_mod)+1):
        for i_sampln in np.arange(1., (nsampl_orig)+1):
            if np.mod(np.floor(matdiv(matcompat.max(np.mod(field_common[int(i_nom)-1], frequencies[int(i_sampln)-1]), (frequencies[int(i_sampln)-1]-np.mod(field_common[int(i_nom)-1], frequencies[int(i_sampln)-1]))), fps)), 2.) == 0.:
                probable_mod[int(i_nom)-1,int(i_sampln)-1] = np.mod(matcompat.max(np.mod(field_common[int(i_nom)-1], frequencies[int(i_sampln)-1]), (frequencies[int(i_sampln)-1]-np.mod(field_common[int(i_nom)-1], frequencies[int(i_sampln)-1]))), fps)
            else:
                probable_mod[int(i_nom)-1,int(i_sampln)-1] = 15.-np.mod(matcompat.max(np.mod(field_common[int(i_nom)-1], frequencies[int(i_sampln)-1]), (frequencies[int(i_sampln)-1]-np.mod(field_common[int(i_nom)-1], frequencies[int(i_sampln)-1]))), fps)
                
            
            
        
    #%probable_mod
    #%pause;
    [pmr, pmc] = matcompat.size(probable_mod)
    probable_mod1 = np.zeros(pmc)
    #%frequencies
    common = np.zeros(num_of_mod)
    #%fps
    field_common
    #%frequencies
    #%probable_mod
    for i_nom in np.arange(1., (num_of_mod)+1):
        for j in np.arange(1., (nsampl_orig)+1):
            common[int(i_nom)-1] = common[int(i_nom)-1]+numel(intersect(probable_mod[int(i_nom)-1,int(j)-1], sampling_all[:,int(j)-1]))
            
        common[int(i_nom)-1]
        if common[int(i_nom)-1] >= 1.*nsampl_orig:
            one_frequency[int(i_nom)-1] = field_common[int(i_nom)-1]
            field_common[int(i_nom)-1]
            probable_mod1[0:pmc] = probable_mod[int(i_nom)-1,0:pmc]
            #%     disp('The Vibration Frequency of the machine:::::::::::::::');
            #%one_frequency
            #%pause;
        else:
            #%   [pmr,pmc]=size(probable_mod);
            one_frequency[int(i_nom)-1] = -1.
            #%disp('hi')
            #%  probable_mod1=zeros(pmc);
            
        
        
    return [one_frequency, probable_mod1]