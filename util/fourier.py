#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

@author: jlvilas
"""
import numpy as np
import util.utilFuntions as uf




def lowPassFilter(img, Digfreq, sampling, fourier = False):
    # This fucntion computes the low pass filtered image of an input image "img"
    # The parameter Digfreq defines the filtering frequency
    # Sampling is the sampling rate or also called pixel size
    # The boolean parameter fourier defines is the input image is real or a Fourier
    # Transform

    if fourier is False:
        lowpassImg = np.fft.fft2(img)
        lowpassImg = np.fft.fftshift(lowpassImg)
    else:
        lowpassImg = img

    freq = np.fft.fftfreq(lowpassImg.shape[1], d=sampling)
    freq = freq[freq > 0]

    dims = img.shape

    # Low pass filtering the original image
    idx = np.arange(0, len(freq), 1)
    idx = idx[freq <= Digfreq]
    cutoff = idx[-1]

    # Filtering the image
    mask = uf.create_circular_mask(dims[0], dims[1], radius=cutoff)
    lowpassImg = np.multiply(mask, lowpassImg)
    lowpassImg = np.real(np.fft.ifft2(np.fft.ifftshift(lowpassImg)))

    return lowpassImg


def highPassFilter(img, Digfreq, sampling, fourier=False):
    # This fucntion computes the high pass filtered image of an input image "img"
    # The parameter Digfreq defines the filtering frequency
    # Sampling is the sampling rate or also called pixel size
    # The boolean parameter fourier defines is the input image is real or a Fourier
    # Transform

    if fourier is False:
        lowpassImg = np.fft.fft2(img)
        lowpassImg = np.fft.fftshift(lowpassImg)
    else:
        lowpassImg = img

    freq = np.fft.fftfreq(lowpassImg.shape[1], d=sampling)
    freq = freq[freq > 0]

    dims = img.shape

    # Low pass filtering the original image
    idx = np.arange(0, len(freq), 1)
    idx = idx[freq >= Digfreq]
    cutoff = idx[-1]

    # Filtering the image
    mask1 = uf.create_circular_mask(dims[0], dims[1], radius=cutoff)
    mask2 = uf.create_circular_mask(dims[0], dims[1])
    lowpassImg = np.multiply(mask1 ^ mask2, lowpassImg)
    lowpassImg = np.real(np.fft.ifft2(np.fft.ifftshift(lowpassImg)))

    return lowpassImg
