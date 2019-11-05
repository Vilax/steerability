#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

@author: jlvilas
"""
import numpy as np
import util.utilFuntions as uf


def FRC(half1, half2, sampling, fourier = True, threshold = 0.143):
    # It estimates the FRC curves from the corresponding fourier transfroms half1, half2.
    # If fourier = True, it means the input half1, and half2 are images, instead
    # of their respective Fourier transforms. Thus, FRC will compute their corresponding
    # Fourier transfroms

    assert half1.shape == half2.shape, 'Halves has not the same dimensions'
    if fourier is False:
        half1 = np.fft.fft2(half1)
        half1 = np.fft.fftshift(half1)
        half2 = np.fft.fft2(half2)
        half2 = np.fft.fftshift(half2)

    freq = np.fft.fftfreq(half1.shape[1], d=sampling)
    freq = freq[freq > 0]

    xdim = half1.shape[0]
    ydim = half1.shape[1]
    dim = int(np.floor(min(xdim, ydim) / 2))

    FSCcurve = np.zeros(dim)
    for i in range(0, dim):
        ring = uf.createRing(xdim, ydim, i)
        f1 = half1[ring]
        f2 = np.conj(half2[ring])
        num = np.sum(np.real(np.multiply(f1, f2)))
        den = np.linalg.norm(f1) * np.linalg.norm(f2)
        FSCcurve[i] = num / den
    # uf.representCurve(range(0, len(FSCcurve)), FSCcurve, 'FSC', True)
    idx = FSCcurve <= threshold
    idx = freq[idx]
    resolution = 1 / idx[0]
    # print('resolution = ', resolution)
    return FSCcurve, resolution


def localFRC(half1, half2, sampling, fourier = True, threshold = 0.143):

    assert half1.shape == half2.shape, 'Halves has not the same dimensions'

    if fourier is False:
        half1 = np.fft.fft2(half1)
        # half1 = np.fft.fftshift(half1)
        half2 = np.fft.fft2(half2)
        # half2 = np.fft.fftshift(half2)

    freq = np.fft.fftfreq(half1.shape[1], d=sampling)
    freq = freq[freq > 0]

    xdim = half1.shape[0]
    ydim = half1.shape[1]
    dim = int(np.floor(min(xdim, ydim) / 2))

    # Gabor filter
    center = [int(xdim/2), int(ydim/2)]
    Y, X = np.ogrid[:xdim, :ydim]
    r = (X - center[0])**2 + (Y-center[1])**2

    for i in range(0, dim):
        print("freq = ", 1/(freq[i]+1e-38))
        sigma = 4*sampling / freq[i]
        # gauss = np.zeros(half1.shape)
        # gauss[256,256] = 1
        gauss = np.exp( -r**2/(2*sigma**2) )
        # uf.representImg(gauss,'Gauss', False)
        # ring = uf.createRing(xdim, ydim, i)
        aux = np.multiply( np.sqrt(r), np.cos(np.arctan2(Y, X ) ))
        gauss = np.multiply(gauss, np.exp(-1j * i * aux))
        gauss = np.fft.fft2(gauss)
        # f1 = np.fft.ifft2(np.multiply(np.multiply(half1, gauss), ring) )
        f1 = np.multiply(half1, gauss)
        f1 = np.fft.ifftshift(f1)
        # f2 = np.fft.ifft2(np.multiply(np.multiply(half2, gauss), ring) )
        f2 = np.multiply(half2, gauss)
        f2 = np.conj(f2)
        f2 = np.fft.ifftshift(f2)

        num = np.real(np.multiply(f1, f2))
        den = np.linalg.norm(f1) * np.linalg.norm(f2)
        print(num/den)

        FSC = np.real(np.fft.ifft2(np.true_divide(num, den)))
        print(FSC)
        uf.representImg(FSC,'FSC', True)

    resolution = 2


    # uf.representCurve(range(0, len(FSCcurve)), FSCcurve, 'FSC', True)
    return resolution






def estimateGaborFRC(half1, half2, sigma, freq, threshold, HPF=True):
    assert half1.shape == half2.shape, 'Halves has not the same dimensions'

    xdim = half1.shape[0]
    ydim = half1.shape[1]
    dim = min(xdim, ydim) / 2

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
        num = np.sum(np.real(np.multiply(f1, f2)))
        den = np.linalg.norm(f1) * np.linalg.norm(f2)
        FSC[i] = num / den
    #     if i % 50 ==0:
    #         representImg(ring, 'gabor', False)

    idx = FSC <= threshold
    idx = freq[idx]
    resolution = 1 / idx[0]

    return FSC, resolution


def gaborFunctionFourier(sigma, omega_x, omega_y, omega_0):
    gabor = (omega_x - omega_0[0]) ** 2 + (omega_y - omega_0[1]) ** 2
    gabor = np.exp(-0.5 * sigma ** 2 * gabor)

    return gabor


def estimateFSC(half1, half2, freq, threshold, HPF=True):
    assert half1.shape == half2.shape, 'Halves has not the same dimensions'

    dims = len(half1.shape)
    if dims == 2:
        xdim = half1.shape[0]
        ydim = half1.shape[1]
        dim = int(np.floor(min(xdim, ydim) / 2))
    if dims == 3:
        xdim = half1.shape[0]
        ydim = half1.shape[1]
        zdim = half1.shape[2]
        dim = int(np.floor(min(xdim, ydim, zdim) / 2))
    else:
        print('Dimension Error: the FSc cannot be calculated')

    F1 = half1
    F2 = np.conj(half2)

    FSC = np.zeros(dim)
    for i in range(0, dim):

        if HPF is True:
            ring = uf.createShellHPF(xdim, ydim, i)
        else:
            ring = uf.createShell(xdim, ydim, i)
        f1 = F1[ring]
        f2 = F2[ring]
        num = np.sum(np.real(np.multiply(f1, f2)))
        den = np.linalg.norm(f1) * np.linalg.norm(f2)
        FSC[i] = num / den
    idx = (FSC <= threshold)
    idx = freq[idx]
    resolution = 1 / idx[0]
    print('resolution = ', resolution)

    return FSC, resolution
