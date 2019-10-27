#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

@author: jlvilas
"""
import numpy as np
import utilFuntions as uf


def estimateGaborFRC(half1, half2, sigma, freq, threshold, HPF = True):
    assert half1.shape == half2.shape, 'Halves has not the same dimensions'

    xdim = half1.shape[0]
    ydim = half1.shape[1]
    dim = min(xdim, ydim)/2

    F1 = half1
    F2 = np.conj(half2)

    FSC = np.zeros(np.int(dim))

    for i in range(0, np.int(dim)):
        if HPF is True:
            ring = uf.createGaussianHPF(dim, sigma, i)
        else:
            ring = uf.createGaussian(dim, sigma, i)
        f1 = np.multiply(F1, ring)
        f2 = np.multiply(F2, ring)
        num = np.sum( np.real( np.multiply( f1, f2) ) )
        den = np.linalg.norm(f1) * np.linalg.norm(f2)
        FSC[i] = num/den
    #     if i % 50 ==0:
    #         representImg(ring, 'gabor', False)
    idx = FSC <= threshold
    idx = freq[idx]
    resolution = 1 / idx[0]
    # uf.representCurve(freq, FSC, 'FSC', True)

    # print('resolution = ', resolution)

    return FSC, resolution


def estimateFRC(half1, half2, freq, threshold, HPF = True):
    assert half1.shape == half2.shape, 'Halves has not the same dimensions'

    xdim = half1.shape[0]
    ydim = half1.shape[1]
    dim = int(np.floor(min(xdim, ydim)/2))

    F1 = half1
    F2 = np.conj(half2)

    FSC = np.zeros(dim)
    for i in range(0, dim):

        if HPF is True:
            ring = uf.createRingHPF(xdim, ydim, i)
        else:
            ring = uf.createRing(xdim, ydim, i)
        f1 = F1[ring]
        f2 = F2[ring]
        num = np.sum( np.real( np.multiply( f1, f2) ) )
        den = np.linalg.norm(f1) * np.linalg.norm(f2)
        FSC[i] = num/den

    idx = FSC <= threshold
    idx = freq[idx]
    resolution = 1/idx[0]
    # print('resolution = ', resolution)

    return FSC, resolution


def estimateLocalFRC(half1, half2, freq, threshold, HPF = True):
    assert half1.shape == half2.shape, 'Halves has not the same dimensions'

    xdim = half1.shape[0]
    ydim = half1.shape[1]
    dim = int(np.floor(min(xdim, ydim) / 2))

    F1 = half1
    F2 = np.conj(half2)

    freq_aux = np.fft.fftshift(freq)
    print(freq.shape)
    freq_x = np.matlib.repmat(freq_aux, len(freq_aux), 1 )
    freq_y = np.transpose(np.fliplr(freq_x))

    sigma = 10
    filterSize = min(xdim, ydim)/2

    for i in range(0, ydim):
        for k in range(0, xdim):
            exp_term = np.exp(-(i*freq_y + k*freq_x))




    # FSC = np.zeros(dim)
    # for i in range(0, dim):
    #
    #     if HPF is True:
    #         ring = uf.createRingHPF(xdim, ydim, i)
    #     else:
    #         ring = uf.createRing(xdim, ydim, i)
    #     f1 = F1[ring]
    #     f2 = F2[ring]
    #     num = np.sum(np.real(np.multiply(f1, f2)))
    #     den = np.linalg.norm(f1) * np.linalg.norm(f2)
    #     FSC[i] = num / den
    #
    # idx = FSC <= threshold
    # idx = freq[idx]
    # resolution = 1 / idx[0]

    return freq




