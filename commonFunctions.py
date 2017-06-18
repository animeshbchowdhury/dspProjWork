import numpy as np
import numpy.mod as nmod
import numpy.abs as nabs
import numpy.arange as narange
import numpy.round as nround
import numpy.zeros as nzeros
import numpy.mean as nmean
import numpy.dot as ndot
import numpy.floot as nfloor
import matplotlib
import scipy as scp
from math import *

def velocityPendCenter(pend_centers):
    time = 1. / 30.
    [r, c] = np.shape(pend_centers)
    VelocityMarker = nzeros([r,c],dtype=float)
    for i in narange(2., (r) + 1):
        VelocityMarker[int(i) - 1, int((c - 1.)) - 1] = (
            pend_centers[int(i) - 1, int((c - 1.)) - 1] - pend_centers[int((i - 1.)) - 1, int((c - 1.)) - 1])/ time
        VelocityMarker[int(i) - 1, int(c) - 1] = (
            pend_centers[int(i) - 1, int(c) - 1] - pend_centers[int((i - 1.)) - 1, int(c) - 1])/ time

    return VelocityMarker


def calculateFrequency(d11, z):
    colSize = d11.shape[1]
    diff = nzeros([3,colSize],dtype=float)
    diff[0,:] = nabs(nround((d11[int((z-3.))-1]-d11[int((z-2.))-1])))
    diff[1,:] = nabs(nround((d11[int((z-2.))-1]-d11[int((z-1.))-1])))
    diff[2,:] = nabs(nround((d11[int((z-1.))-1]-d11[int(z)-1])))
    w = scp.stats.mode(diff)
    return w


def calculateTripleFrequency(remainder, strobe):
    a = np.array([1., 2., 3.])
    remainder1 = remainder.copy()
    remainder2 = remainder1.copy() - 0.6
    strobeLen = strobe.shape[1]
    aLen = a.shape[1]
    fre_1 = nzeros([strobeLen,aLen],dtype=float)
    for i in narange(1., (strobeLen) + 1):
        for j in narange(1., (aLen) + 1):
            fre_1[int(i) - 1, int(j) - 1] = ndot(strobe[int(i) - 1], a[int(j) - 1])

    t = np.shape(fre_1)
    iter = 1.
    remainder2Len = remainder2.shape[1]
    freq_plus = nzeros([t[0,0]*t[0,1],remainder2Len])
    freq_minus = nzeros([t[0,0]*t[0,1],remainder2Len])
    freq_plus_shift_negative1 = nzeros([t[0,0]*t[0,1],remainder2Len])
    freq_plus_shift_negative2 = nzeros([t[0,0]*t[0,1],remainder2Len])

    for i in narange(1., (t[0, 0]) + 1):
        for j in narange(1., (t[0, 1]) + 1):
            fre_11 = fre_1[int(i) - 1, int(j) - 1].copy()
            remainder_1 = remainder2[int(i) - 1, :].copy()
            for k in narange(1., remainder2Len + 1):
                if remainder_1[0, int(k) - 1] > 0.:
                    freq_plus[int(iter) - 1, int(k) - 1] = fre_11 + remainder_1[0, int(k) - 1]
                    freq_minus[int(iter) - 1, int(k) - 1] = fre_11 - remainder_1[0, int(k) - 1]
                    freq_plus_shift_negative1[int(iter) - 1, int(k) - 1] = fre_11 + 30. - remainder_1[0, int(k) - 1]
                    freq_plus_shift_negative2[int(iter) - 1, int(k) - 1] = fre_11 - 30. - remainder_1[0, int(k) - 1]

            iter = iter + 1.

    est_fre = np.vstack((freq_plus, freq_minus, freq_plus_shift_negative1, freq_plus_shift_negative2))
    frequency1 = nround(est_fre)
    frequency2 = nfloor(est_fre)
    return [frequency1, frequency2]


def findPeakPyVersion(dataSample):
    peakPoints = []
    for i in narange(1.,dataSample.size()):
        if((dataSample[0,int(i-1)] < dataSample[0,int(i)]) and (dataSample[0,int(i)] < dataSample[0,int(i+1)])):
            peakPoints.append(dataSample[int(i)])
    return np.array([peakPoints])


def checkValidity_MultiFrquency(X1, Y1, f1, NFFT1):
    d1_list = []
    d2_list = []
    g2_list = []

    peak1 = np.sort(findPeakPyVersion((2. * X1[0:int(NFFT1 / 2. + 1.)])), 'descend')
    peak1Len = peak1.shape[1]
    chk_val1 = nzeros([peak1Len,1])
    index1 = nzeros([peak1Len,1])
    for i in narange(1., (peak1Len + 1)):
        chk_val1[int(i) - 1, :] = peak1[int(i) - 1]/ (nmean(peak1[int(i):peak1Len]))
        if chk_val1[int(i) - 1, :] > 2.:
            index1[int(i) - 1, :] = np.nonzero(X1[0:NFFT1 / 2. + 1.] == peak1[:, int(i) - 1] / 2.)
            d1_list.append(f1[index1[int(i) - 1, :]])

    mean1 = nmean(peak1[0:peak1Len])
    k = 1.
    peak2 = np.sort(findPeakPyVersion((2. * Y1[0:NFFT1 / 2. + 1.])), 'descend')
    peak2Len = peak2.shape[1]
    chk_val = nzeros([peak2Len, 1])
    index2 = nzeros([peak2Len, 1])
    for i in narange(1., peak2Len):
        chk_val[int(i) - 1, :] = peak2[0, int(i) - 1]/nmean(peak2[int(i):peak2Len])
        if chk_val[int(i) - 1, :] > 3.:
            g2_list.append(chk_val[int(i) - 1, :])
            index2[int(i) - 1, :] = np.nonzero((Y1[0:NFFT1 / 2. + 1.] == peak2[:, int(i) - 1] / 2.))
            d2_list.append(f1[index2[int(i) - 1, :]])

    d2 = np.array([d2_list])
    d1 = np.array([d1_list])
    g2 = np.array([g2_list])

    if g2 >= 0.:
        g1 = nmean(g2)
    else:
        g1 = nmean(g2)

    return [g1, d1, d2]


def pyNextPow2(num):
    return pow(2, ceil(log(num)/log(2))) # Returns next higher power of two

def FFT_MultiFrequency_update(s1, s2):
    Fs = np.array([[31.]])
    #% Sampling frequency
    T = 1./Fs
    #% Sample time
    L = 512.
    #% Length of signal
    t = ndot(narange(0., L), T)

    NFFT1 = float(pow(2,pyNextPow2(L)))
    #% Next power of 2 from length of y
    s1[0,:] = s1[0,:]-nmean(s1[0,:])
    s2[0,:] = s2[0,:]-nmean(s2[0,:])
    X1 = nabs(np.fft(s1, NFFT1)/ L)
    Y1 = nabs(np.fft(s2, NFFT1)/ L)
    f1 = ndot(Fs/2., np.linspace(0., 1., (NFFT1/2.+1.)))

    [g1, d1, d2] = checkValidity_MultiFrquency(X1, Y1, f1, NFFT1)
    return [d1, X1, Y1, f1, NFFT1, d2, g1]


def mfreq_simulate2(frequency):
    nfreq = 4.
    nsampl = 33.
    sampling_option = 0.
    cam_fps = 15.
    prime1 = np.array([61.,67., 71., 73., 79., 83., 89., 97., 101., 103., 107.,
                        109., 113., 127., 131., 137., 139., 149., 151., 157., 163.,
                        167., 173., 179., 181., 191., 193., 197., 199., 211., 223., 227., 229.])
    sFile = open('strobe_file.txt','w+')
    freq = np.random.rand(1., nfreq)
    frequencies = np.array([70., 100., 170., 230.])
    if sampling_option == 1.:
        sampling = nzeros([nsampl])
        for i in narange(1., (nsampl) + 1):
            stri = "give"+str(i)+"th frequency of strobe"
            print stri
            sampling[int(i) - 1] = float(raw_input())

    else:
        sampling = prime1.copy()

    print sampling
    remainder = nzeros([int(nfreq),int(nsampl)])
    remainder2 = remainder.copy()
    for j in narange(1., (nsampl) + 1):
        for i in narange(1., (nfreq) + 1):
            if nmod(nfloor((np.minimum(nmod(frequencies[int(i) - 1], sampling[int(j) - 1]), (
                sampling[int(j) - 1] - nmod(frequencies[int(i) - 1], sampling[int(j) - 1])))/cam_fps)), 2.) == 0.:
                remainder[int(i) - 1, int(j) - 1] = nmod(
                    np.minimum(nmod(frequencies[int(i) - 1], sampling[int(j) - 1]),
                                  (sampling[int(j) - 1] - np.mod(frequencies[int(i) - 1], sampling[int(j) - 1]))),
                    cam_fps)
            else:
                remainder[int(i) - 1, int(j) - 1] = 15. - nmod(
                    np.minimum(nmod(frequencies[int(i) - 1], sampling[int(j) - 1]),
                                  (sampling[int(j) - 1] - np.mod(frequencies[int(i) - 1], sampling[int(j) - 1]))),
                    cam_fps)

        remainder2[:, int(j) - 1] = np.sort(remainder[:, int(j) - 1])


    sFile.write("%f\n" % nfreq)
    sFile.write("%f\n" % nsampl)

    for j in narange(1., (nsampl) + 1):
        sFile.write("%f " % sampling[int(j) - 1])

    sFile.write("\n")
    for i in narange(1., (nfreq) + 1):
        for j in np.arange(1., (nsampl) + 1):
            sFile.write('%8.2f '% remainder2[int(i) - 1, int(j) - 1])
        sFile.write('\n')

    sFile.close()
    return [frequencies, sampling]


def mfreq_solve12():
    s = open('strobe_file.txt', 'r')
    nfreq = fscanf(s, '%f', 1.)
    nsampl = fscanf(s, '%f', 1.)
    for j in np.arange(1., (nsampl) + 1):
        frequencies[int(j) - 1] = fscanf(s, '%f', 1.)

    for i in np.arange(1., (nfreq) + 1):
        for j in np.arange(1., (nsampl) + 1):
            sampling[int(i) - 1, int(j) - 1] = fscanf(s, '%f', 1.)

    fclose(s)
    nfreq_orig = nfreq
    nsampl_orig = nsampl
    for i1 in np.arange(1., (nfreq_orig) + 1):
        il_all = 1.
        one_frequency = -1.
        while one_frequency < 0.:
            if nfreq < 8.:
                nsampl = matcompat.max((nfreq * 3.), nsampl)

            margin = nsampl_orig - nsampl
            if margin > 0.:
                rand_start = np.floor(np.dot(np.random.rand(1., 1.), margin - 1.)) + 1.
            else:
                rand_start = 1.

            [one_frequency, probable_mod1] = find_frequency(nfreq, nsampl, nsampl_orig, frequencies, sampling,
                                                            sampling[0:nfreq,
                                                            int(rand_start) - 1:rand_start - 1. + nsampl])

        num_freq = numel(one_frequency)
        for i_num in np.arange(1., (num_freq) + 1):
            if one_frequency[int(i_num) - 1] > 0.:
                freq_all[int(i1) - 1] = one_frequency[int(i_num) - 1]
                freq_all_all1[int(il_all) - 1] = one_frequency[int(i_num) - 1]
                il_all = il_all + 1.

        np.sort(freq_all)
        # %pause;
        # %probable_mod;
        for j in np.arange(1., (nsampl_orig) + 1):
            count = 0.
            for k in np.arange(1., (nfreq - 1.) + 1):
                if count == 1. or sampling[int(k) - 1, int(j) - 1] == probable_mod1[int(j) - 1]:
                    count = 1.
                    sampling[int(k) - 1, int(j) - 1] = sampling[int((k + 1.)) - 1, int(j) - 1]

        # %sampling;
        nfreq = nfreq - 1.
        # %sampling
        # %pause;

    # %freq_all;
    # %nfreq;
    # % for j=1:nsampl_orig
    # %      for i=1:nfreq_orig
    # %     remainder(i,j)=mod(freq_all(i),frequencies(j));
    # %      end
    # %     remainder2(:,j)=sort(remainder(:,j));
    # % end
    np.sort(freq_all)
    freq_all_all = np.unique(np.sort(freq_all_all1))
    return [freq_all, freq_all_all]


def find_frequency(nfreq, nsampl, nsampl_orig, frequencies, sampling_all, sampling):
    # Local Variables: i_nom, pmr, mm1, sampling, nsampl, one_frequency, num_of_mod, pmc, nfreq, nsampl_orig, cardinality, test_set, i_sampln, probable_mod, tempfkm, field_common1, field, fields_to_check, fps, field_common, frequencies, probable_mod1, sampling_all, mm2, num_fps, i, k, j, m, eta, search_limit, common, k_fps
    # Function calls: size, rand, intersect, floor, min, find_frequency, zeros, numel, unique, mod
    fps = 15.
    eta = np.floor(matdiv(nsampl, nfreq))
    for j in np.arange(1., (nsampl) + 1):
        test_set[int(j) - 1] = sampling[0, int(j) - 1]

    if eta == 1.:
        eta = 2.

    fields_to_check = np.zeros(eta)
    # % fields_to_check=[1 2];
    # % end;
    # % fields_to_check=fields_to_check-1;
    while numel(np.unique(fields_to_check)) < eta:
        fields_to_check = np.floor(np.dot(np.random.rand(1., eta), nsampl)) + 1.

    # %eta
    # %fields_to_check
    # %pause;
    search_limit = 1.
    for k in np.arange(1., (eta) + 1):
        search_limit = np.dot(search_limit, frequencies[int(fields_to_check[int(k) - 1]) - 1])

    if search_limit > 2000.:
        search_limit = 2000.

    for k in np.arange(1., (eta) + 1):
        m = 2.
        field[int(k) - 1, 1] = test_set[int(fields_to_check[int(k) - 1]) - 1]
        # %pause
        num_fps[int(k) - 1] = np.floor(matdiv(frequencies[int(fields_to_check[int(k) - 1]) - 1] / 2., fps))
        # %frequencies(fields_to_check(k))
        # %num_fps(k)
        # %pause;
        for k_fps in np.arange(1., (num_fps[int(k) - 1] + 3.) + 1):
            if np.mod(k_fps, 2.) == 0.:
                tempfkm = np.dot(k_fps - 1., fps) + fps - field[int(k) - 1, 1]
                if tempfkm < frequencies[int(fields_to_check[int(k) - 1]) - 1] / 2.:
                    field[int(k) - 1, int(m) - 1] = tempfkm


            else:
                tempfkm = np.dot(k_fps - 1., fps) + field[int(k) - 1, 1]
                if tempfkm < frequencies[int(fields_to_check[int(k) - 1]) - 1] / 2.:
                    field[int(k) - 1, int(m) - 1] = tempfkm

            # %tempfkm
            # %pause
            mm1 = 1.
            mm2 = 0.
            while -tempfkm + np.dot(mm1, frequencies[int(fields_to_check[int(k) - 1]) - 1]) < search_limit:
                field[int(k) - 1, int((m + 1.)) - 1] = -tempfkm + np.dot(mm1, frequencies[
                    int(fields_to_check[int(k) - 1]) - 1])
                mm1 = mm1 + 1.
                m = m + 1.

            while tempfkm + np.dot(mm2, frequencies[int(fields_to_check[int(k) - 1]) - 1]) < search_limit:
                field[int(k) - 1, int((m + 1.)) - 1] = tempfkm + np.dot(mm2, frequencies[
                    int(fields_to_check[int(k) - 1]) - 1])
                mm2 = mm2 + 1.
                m = m + 1.

        field[int(k) - 1, 0] = m - 2.

    # %field
    # %search_limit
    # %pause;
    field_common1[0, 0] = 0.
    for i in np.arange(1., (eta - 1.) + 1):
        cardinality = numel(intersect(field[int(i) - 1, 1:field[int(i) - 1, 0] + 1.],
                                      field[int((i + 1.)) - 1, 1:field[int((i + 1.)) - 1, 0] + 1.]))
        field_common1[0:cardinality] = intersect(field[int(i) - 1, 1:field[int(i) - 1, 0] + 1.],
                                                 field[int((i + 1.)) - 1, 1:field[int((i + 1.)) - 1, 0] + 1.])
        field_common1.flatten(1)
        if numel(field_common1) > 0.:
            field[int((i + 1.)) - 1, 0] = numel(field_common1.flatten(1))
            field[int((i + 1.)) - 1, 1:numel[field_common1[:]] + 1.] = field_common1.flatten(1)
        else:
            break

    if numel(field_common1) < 1.:
        field_common = 0.
    else:
        field_common = field_common1.flatten(1)

    # %field_common
    # %pause;
    # %frequencies
    fps
    # %pause
    num_of_mod = numel(field_common)
    for i_nom in np.arange(1., (num_of_mod) + 1):
        for i_sampln in np.arange(1., (nsampl_orig) + 1):
            if np.mod(np.floor(matdiv(
                    matcompat.max(np.mod(field_common[int(i_nom) - 1], frequencies[int(i_sampln) - 1]), (
                        frequencies[int(i_sampln) - 1] - np.mod(field_common[int(i_nom) - 1],
                                                                frequencies[int(i_sampln) - 1]))), fps)), 2.) == 0.:
                probable_mod[int(i_nom) - 1, int(i_sampln) - 1] = np.mod(
                    matcompat.max(np.mod(field_common[int(i_nom) - 1], frequencies[int(i_sampln) - 1]), (
                    frequencies[int(i_sampln) - 1] - np.mod(field_common[int(i_nom) - 1],
                                                            frequencies[int(i_sampln) - 1]))), fps)
            else:
                probable_mod[int(i_nom) - 1, int(i_sampln) - 1] = 15. - np.mod(
                    matcompat.max(np.mod(field_common[int(i_nom) - 1], frequencies[int(i_sampln) - 1]), (
                    frequencies[int(i_sampln) - 1] - np.mod(field_common[int(i_nom) - 1],
                                                            frequencies[int(i_sampln) - 1]))), fps)

    # %probable_mod
    # %pause;
    [pmr, pmc] = matcompat.size(probable_mod)
    probable_mod1 = np.zeros(pmc)
    # %frequencies
    common = np.zeros(num_of_mod)
    # %fps
    field_common
    # %frequencies
    # %probable_mod
    for i_nom in np.arange(1., (num_of_mod) + 1):
        for j in np.arange(1., (nsampl_orig) + 1):
            common[int(i_nom) - 1] = common[int(i_nom) - 1] + numel(
                intersect(probable_mod[int(i_nom) - 1, int(j) - 1], sampling_all[:, int(j) - 1]))

        common[int(i_nom) - 1]
        if common[int(i_nom) - 1] >= 1. * nsampl_orig:
            one_frequency[int(i_nom) - 1] = field_common[int(i_nom) - 1]
            field_common[int(i_nom) - 1]
            probable_mod1[0:pmc] = probable_mod[int(i_nom) - 1, 0:pmc]
            # %     disp('The Vibration Frequency of the machine:::::::::::::::');
            # %one_frequency
            # %pause;
        else:
            # %   [pmr,pmc]=size(probable_mod);
            one_frequency[int(i_nom) - 1] = -1.
            # %disp('hi')
            # %  probable_mod1=zeros(pmc);

    return [one_frequency, probable_mod1]