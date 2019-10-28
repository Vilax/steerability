# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:24:03 2019

@author: tommy
"""

import numpy as np
import math
import scipy.linalg as la
import matplotlib.pyplot as plt
# from scipy.spatial.transform import Rotation as R
from numpy import linalg as npla
# from scipy.spatial.transform import Rotation as R
# from scipy import ndimage as ndimg
import sys
import util.steering3D as steer
import util.testFunctions as Tsf
import util.resolution as res
import util.filter3d as filt
import util.polynomial3d as poly3d
import mrcfile
import util.utilFuntions as uf

sys.path.insert(0, '../util')

# testing steering a filter derived from a gaussian of the radius times an even
# polynomial

N = 4
cap_size = np.pi / 9.0
r0 = 20
sigma = 5
filtSize = 200

direction = [0, 0, 1]

wavelength = 10
angle = 1
xdim, ydim, zdim = 201, 201, 201
direction = [1, 1, 0]

mu, sigma = 0, 0.5
import time

filter_range = np.arange(-filtSize, filtSize, 1)
X, Y, Z = np.meshgrid(filter_range, filter_range, filter_range)
start = time.time()
a = uf.s1(X)
end = time.time()
print(end - start)
start = time.time()
b = uf.s2(X)
end = time.time()
print(end - start)
start = time.time()
c = uf.s3(X)
end = time.time()
print(end - start)

start = time





# Defining a fringe pattern
vol1 = Tsf.define_sinusoidal_pattern(wavelength, xdim, ydim, zdim, direction)
vol2 = Tsf.define_sinusoidal_pattern(wavelength, xdim, ydim, zdim, direction)

# adding noise to the fringe pattern
imgNoise1 = Tsf.add_gaussian_noise(vol1, mu, sigma)
imgNoise2 = Tsf.add_gaussian_noise(vol2, mu, sigma)

# # Ensuring odd dimensions for the fft
imgNoise1 = uf.paddingImageIfIsOdd(imgNoise1)
imgNoise2 = uf.paddingImageIfIsOdd(imgNoise2)

# Fourier transform of the image
imgNoise_fft1 = np.fft.fftn(imgNoise1)
imgNoise_fft1 = np.fft.fftshift(imgNoise_fft1)
imgNoise_fft2 = np.fft.fftn(imgNoise2)
imgNoise_fft2 = np.fft.fftshift(imgNoise_fft2)

uf.representImg((imgNoise1[:, :, 100]), 'fringes', True)
filtSize = imgNoise1.shape[0]
dirFilt = steer.directionalFilter3D(N, cap_size, r0, sigma, direction, filtSize)
uf.representImg((dirFilt[:, :, 100]), 'Directional Filter', True)

# Aplying a filter to the image
fft_filt1 = np.multiply(dirFilt, imgNoise_fft1)
fft_filt2 = np.multiply(dirFilt, imgNoise_fft2)

timestep = 1
FreqCompRows = np.fft.fftfreq(imgNoise1.shape[0], d=timestep)
FreqCompCols = np.fft.fftfreq(imgNoise1.shape[1], d=timestep)
freq = FreqCompCols[FreqCompCols > 0]
# print(1/freq)
FSC, resolution = res.estimateFSC(imgNoise_fft1, imgNoise_fft2, freq, 0.143, True)
