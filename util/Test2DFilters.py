#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

# @author: jlvilas
"""
import numpy as np
import scipy as sp
from scipy import fftpack, interpolate
import matplotlib.pyplot as plt
import steer2D as ste
import random
from skimage.transform import rotate
import mrcfile

xdim = 1000
ydim = 1000
wavelength = 20
angle = 0
mu = 0
sigma = 0.5

# Defining the directional filter
N = 6
capSize = np.pi/4.0

Theta = capSize

f, u, bCos, theta = ste.steer2dGeneral(Theta, N)

normVals = ste.computeNormalizationConstant2D(N)

u = np.true_divide(u, normVals)
p = np.flipud(u) #its the same vector with an scale factor

aux = np.sum(p)
p = p/aux

#TODO: change x=200
xmin = -200.0
xmax = 200.0
ymin = -200.0
ymax = 200.0

bigImgLims = [xmin, xmax, ymin, ymax]
# filt = ste.poly2filter(p, bigImgLims, 0.04, 0)

r0 = 20
sigma = 5.0
filt = ste.makeSteerFilt(100, f, theta, r0, sigma)

nonzeroCoefficients = np.arange(-N, N+1, 2)
nonzeroCoeffBool = ste.index2boolean(nonzeroCoefficients, N)
M = ste.validateNonzeroCoefficients(nonzeroCoeffBool, N)

steerBasisMat = ste.computeSteerBasis2D(N, nonzeroCoeffBool)
Phi = np.linalg.inv(steerBasisMat)

steerFiltVol = ste.makeSteerBasis(filt, M)

steerFiltVol2 = np.float32(steerFiltVol)
fileToSave = mrcfile.new('tmp.mrc', overwrite = True)
fileToSave.set_data(steerFiltVol2)

testThetaArr = [np.pi/5.5, 3.7*np.pi/5.0, 4.7*np.pi/4.0, 5*np.pi/3.0]
ntestTheta = len(testThetaArr)
for itheta in [0]:#1:ntestTheta:
    # extract desired rotation
    thetaTest = testThetaArr[itheta]
    thetaTestDeg = thetaTest * 180 / np.pi

    # compute manually with imrotate
    bruteForceRotatedFilt = np.transpose(rotate(filt, -thetaTestDeg))

    # compute e ^ itheta powers
    dirVals = ste.computeSteerDirection2D(N, thetaTest, nonzeroCoeffBool)

    k = np.real(np.dot(Phi, dirVals))

    # # compute steered filter
    steerCoeff = np.reshape(k, [1, 1, len(k)])
    steeredFilt = np.sum((np.multiply(steerFiltVol, steerCoeff)), axis=2)

    plt.figure()
    plt.imshow(bruteForceRotatedFilt)
    plt.colorbar()
    plt.figure()
    plt.imshow(steeredFilt)
    plt.colorbar()
    plt.show()

    approxError = (bruteForceRotatedFilt - steeredFilt)
    tmpDenom = bruteForceRotatedFilt
    tmpDenom[tmpDenom == 0] = 1


# visualizing steering function
step = 0.1
thetaRange = np.arange(0,3*np.pi+step,step)
nrotations = len(thetaRange)
steerCoefficientsArray = np.zeros(thetaRange.shape)
steerCoefficientsArray2 = np.zeros(thetaRange.shape)
for irotation in range(0,nrotations):
    theta = thetaRange[irotation]
    dirVals = ste.computeSteerDirection2D(N, theta, nonzeroCoeffBool)
    print(Phi.shape)
    print(dirVals.shape)
    k = np.real(np.dot(Phi,dirVals))
    print(k.shape)
    steerCoefficientsArray[irotation] = k[0]
    steerCoefficientsArray2[irotation] = k[1]


plt.figure()
plt.plot(thetaRange[:], steerCoefficientsArray[:])
plt.show()
plt.figure()
plt.plot(thetaRange[:], steerCoefficientsArray2[:])
plt.show()
