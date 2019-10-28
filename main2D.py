#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

# @author: jlvilas
"""
import util.testFunctions as TstFun
import numpy as np
import util.utilFuntions as uf
import util.resolution as res
import scipy as sp
from scipy import fftpack, interpolate
import matplotlib.pyplot as plt
import util.steering2D as ste


# Defining image parameters
xdim = 512
ydim = 512
wavelength = 5
angle = 0
mu, sigma = 0, 0

# Reading Image:

# Defining a fringe pattern
img1 = TstFun.define_sinusoidal_pattern(wavelength, xdim, ydim, 1, 0.0)
img2 = TstFun.define_sinusoidal_pattern(wavelength, xdim, ydim, 1, 0.0)

# adding noise to the fringe pattern
imgNoise1 = TstFun.add_gaussian_noise(img1, mu, sigma)
imgNoise2 = TstFun.add_gaussian_noise(img2, mu, sigma)

# uf.representImg(imgNoise1, 'original Image Img1', False)
# uf.representImg(imgNoise2, 'original Image Img2', False)

# Ensuring odd dimensions for the fft
imgNoise1 = uf.paddingImageIfIsOdd(imgNoise1)
imgNoise2 = uf.paddingImageIfIsOdd(imgNoise2)
uf.representImg(imgNoise1, 'original Image Img1', True)

"""
# Fourier transform of the image
imgNoise_fft1 = np.fft.fft2(imgNoise1)
imgNoise_fft1 = np.fft.fftshift(imgNoise_fft1)
imgNoise_fft2 = np.fft.fft2(imgNoise2)
imgNoise_fft2 = np.fft.fftshift(imgNoise_fft2)

# uf.representImg(np.log(np.abs(imgNoise_fft2)**2), 'original FFT Img2', False)

timestep = 1
FreqCompRows = np.fft.fftfreq(imgNoise1.shape[0], d=timestep)
FreqCompCols = np.fft.fftfreq(imgNoise1.shape[1], d=timestep)
freq = FreqCompCols[FreqCompCols>0]


# Defining the directional filter
#Filter parameters
r0 = 150
sigma = 5.0
cap = 2.0
N = 4
filtSize = xdim
print("filtSize", filtSize)

direction = np.pi/10
dirFilt = ste.directionalFilter2D(N, cap, r0, sigma, direction, filtSize)

# # Aplying a filter to the image
fft_filt1 = np.multiply(dirFilt, imgNoise_fft1)
fft_filt2 = np.multiply(dirFilt, imgNoise_fft2)
FSC, resolution = res.estimateLocalFRC(imgNoise_fft1, imgNoise_fft2, FreqCompCols, 0.143, True)
print('Normal HPF resolution = ', resolution)
FSC, resolution = res.estimateGaborFRC(imgNoise_fft1, imgNoise_fft2, 2, freq, 0.143, True)
print('Gabor HPF resolution = ', resolution)
FSC, resolution = res.estimateFRC(imgNoise_fft1, imgNoise_fft2, freq, 0.143, False)
print('Normal resolution = ', resolution)
uf.representCurve(freq, FSC, 'FSCnormal', False)
FSC, resolution = res.estimateGaborFRC(imgNoise_fft1, imgNoise_fft2, 2, freq, 0.143, False)
print('Gabor resolution = ', resolution)
uf.representCurve(freq, FSC, 'FSCGabor', True)


# count = 0
# for direc in np.arange(0,np.pi/2, np.pi/18.0):
#     count = count + 1
#     print('direction = ',  direc*180/np.pi, 'number = ', count)
#     dirFilt = ste.directionalFilter2D(N, cap, r0, sigma, direc, filtSize)
#     uf.representImg(dirFilt, "filter"+str(count), False)
#
#     # # Aplying a filter to the image
#     fft_filt1 = np.multiply(dirFilt, imgNoise_fft1)
#     fft_filt2 = np.multiply(dirFilt, imgNoise_fft2)
#
#     # imReal = np.real(np.fft.ifft2(np.fft.ifftshift(fft_filt1)))
#     # uf.representImg(imReal, 'Img Real', False)
#
#     FSC, resolution = res.estimateFSC(fft_filt1, fft_filt2, freq, 0.143, HPF = True)
#     # plt.figure()
#     # plt.plot(freq, FSC)
#     # plt.title("FSC dir"+str(count))
#
# plt.show()
"""