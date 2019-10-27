#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

@author: jlvilas
"""
import numpy as np
import numpy.matlib as npm
import util.utilFuntions as uf
import random
import scipy as sp
from scipy import linalg, interpolate, misc
import skimage as sk
from skimage.transform import rotate


def directionalFilter2D(N, cap, r0, sigma, direction, filtSize):
    # This function gives and image with the directional filter
    # N represent the order of the polynomial of the filter
    # The angle of the filter is defined as pi/cap
    # sigma, si the standard deviation of the gaussian which defines the frequencies
    # direction in radians

    Theta = np.pi / cap

    f, u, bCos, theta = steer2dGeneral(Theta, N)

    normVals = computeNormalizationConstant2D(N)

    u = np.true_divide(u, normVals)
    p = np.flipud(u)  # its the same vector with an scale factor

    aux = np.sum(p)
    p = p / aux

    filtSize = int(filtSize)/2 # in pixels

    filt = makeSteerFilt(filtSize, f, theta, r0, sigma)

    nonzeroCoefficients = np.arange(-N, N + 1, 2)
    nonzeroCoeffBool = index2boolean(nonzeroCoefficients, N)
    M = validateNonzeroCoefficients(nonzeroCoeffBool, N)

    steerBasisMat = computeSteerBasis2D(N, nonzeroCoeffBool)
    Phi = np.linalg.inv(steerBasisMat)

    steerFiltVol = makeSteerBasis(filt, M)

    # compute e ^ itheta powers
    dirVals = computeSteerDirection2D(N, direction, nonzeroCoeffBool)

    k = np.real(np.dot(Phi, dirVals))

    # # compute steered filter
    steerCoeff = np.reshape(k, [1, 1, len(k)])
    steeredFilt = np.abs(np.sum((np.multiply(steerFiltVol, steerCoeff)), axis=2))
    steeredFilt = steeredFilt/np.max(steeredFilt)

    angleCritic, fillingValue = estimateFilterWidth(steeredFilt, direction)
    steeredFilt = maskrippling(steeredFilt, direction, filtSize, angleCritic, fillingValue)
    mask = uf.create_circular_mask(steeredFilt.shape[0], steeredFilt.shape[0])

    steeredFilt = np.multiply(steeredFilt, mask)

    return steeredFilt


def estimateFilterWidth(filter, direction):
    # Center of the image
    center_dir = int(np.floor(0.5 * filter.shape[0]))

    r = (center_dir * 2/3) * np.transpose([np.sin(direction), -np.cos(direction)])

    last_idx_x = np.int(center_dir + r[0])
    last_idx_y = np.int(center_dir + r[1])

    angleCritic = 0

    # valuetest = np.array([])
    ran = np.arange(0, np.pi / 2, np.pi / 180)
    lastAngle = ran[0]
    for theta in ran:
        rotMatrix = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
        r_rotated = np.dot(rotMatrix, r)

        idx_x = np.int(center_dir + r_rotated[0])
        idx_y = np.int(center_dir + r_rotated[1])

        # valuetest = np.append(valuetest, [filter[idx_x, idx_y]])

        if filter[idx_x, idx_y] > filter[last_idx_x, last_idx_y]:
            angleCritic = lastAngle
            value = filter[last_idx_x, last_idx_y]
            break
        lastAngle = theta
        last_idx_x = idx_x
        last_idx_y = idx_y
    return angleCritic, value


def maskrippling(steeredFilt, direction, filtSize, angleCritic, value):
    # Directional filters based on steerability usualy present a rippling
    # this function mask that rippling, resulting in a monotonic and
    # smooth directional filter.

    nn = np.arange(-filtSize, filtSize, 1)
    x, y = np.meshgrid(nn, nn)

    # x = np.true_divide(x, np.sqrt(x ** 2 + y ** 2))

    angleCone = -(np.arctan2(y, x))
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.imshow(angleCone * 180 / np.pi)
    # plt.title('angles')
    # plt.show()

    idx1 = np.abs(angleCone - direction) > angleCritic
    idx2 = np.flipud(np.fliplr(idx1))

    # plt.figure()
    # plt.imshow(idx1)
    # plt.title('idx1')
    # plt.show()

    idx = np.logical_and(idx1, idx2)
    steeredFilt[idx] = value

    return steeredFilt


def steer2dGeneral(Theta, N, *args):
    # N must be maximum bandwidth allowed, nonzeroBool should have length
    # N + 1. If nonzeroBool not provided, will take every other power to
    # be included, starting with the MAXIMUM POWER (fills in every other entry
    # with 1, starting with the rightmost entry)

    varargin = args
    nargin = 2 + len(varargin)
    nonzeroBool = []
    if nargin < 3:
        nonzeroBool = np.zeros((N + 1))
        id = N - np.arange(0, N+2, 2)
        nonzeroBool[id] = 1
    else:
        nonzeroBool = args

    bCos, sqrtSin, theta, dtheta = getBasis(N, 400, nonzeroBool)

    #Get the Gram matrices
    G1 = makeGramMatrix(bCos, theta, dtheta,Theta)
    G2 = makeGramMatrix(bCos, theta, dtheta, np.pi/2)
    # Generalized  eigenvector problem
    d, v = sp.linalg.eig(G1, G2, left=True, right=False)
    min_idx = np.argmin(d)

    u = v[:, min_idx]/np.linalg.norm(v[:, min_idx])
    ut = np.transpose(u)
    constraint = np.dot(np.dot(u, G2), ut)
    aux = np.dot(np.dot(u, G1), ut)
    obj = aux / constraint

    # Make the steerable function and plot it
    f = makeSteerableFunction(v[:, min_idx], bCos)
    f = f / f[0]

    return f, u, bCos, theta


def getBasis(n, nSamples, nonzeroBool):
    #
    dtheta = np.pi / nSamples
    theta = np.arange(0, np.pi+dtheta, dtheta)
    bCos = np.zeros((len(theta), n+1))
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    sinTheta[-1] = 0
    sqrtSin = np.sqrt(sinTheta)
    for k in range(0, n+1):
        bCos[:, k] = np.power(cosTheta, k)
    bCosNorm = np.transpose( np.sqrt( np.diag( np.dot(np.transpose(bCos), bCos)*dtheta) ) )
    bCosNorm = np.matlib.repmat(bCosNorm, len(theta), 1)
    bCos = np.true_divide(bCos, bCosNorm)
    # ind = np.arange(0,len(nonzeroBool))
    aux = nonzeroBool == 1
    bCos = bCos[:, aux]
    return bCos, sqrtSin, theta, dtheta


def makeGramMatrix(bCos, theta, dtheta, thetaE):
    G = 0
    index = theta <= thetaE
    bCos = bCos[index,:]
    G = np.dot(np.transpose(bCos), bCos)*dtheta

    return G


def makeSteerableFunction(v, bCos):
    aux = bCos.shape
    # f = np.zeros( (aux[0], 1) )
    f = np.zeros((aux[0]))
    for k in range(0, len(v)):
        f = f + v[k] * bCos[:, k]
    return f


def computeNormalizationConstant2D(N):
    # % COMPUTENORMALIZATIONCONSTANT2D    Compute normalization constants for 2d
    # % s1 basis
    # % normVals=COMPUTENORMALIZATIONCONSTANT2D(N) computes for up
    # % to bandwidth N the normalization constants for the basis vectors 1, x^2,
    # % x^4,... x^N by saving the integral from 0 to pi/2 of cos^{2m} xdx where m
    # % ranges from 0 to N, where N=2M is assumed even.
    # % result will be an array normVals where normVals(k) = 1/(int_0^pi/2
    # % cos^{4(k-1)}x dx), that is, writing alpha_k = normVals(k+1), that is,
    # % adjusting for 1-indexing in matlab, alpha_k is the
    # % constant s.t. each alpha_k is the coefficient to x^{2k}, for k=0,..,M,
    # % such that the int_0^pi/2 (alpha_k x^2k)^2 dx is 1.

    assert N >= 1, 'Error: N is lesser than 1, N<1'

    N = validateBandwidth(N)
    M = int(N/2)

    integralVals = np.zeros([M+1])
    # % initialize first steps
    integralVals[0] = np.pi/2.0

    # % go through array
    # note that in general int_0^pi/2 cos^m xdx = (m-1)/m * int_0^pi/2
    # % cos^{m-2} x dx = (m-1)(m-3)/m(m-2) int_0^pi/2 cos^{m-4}xdx
    for ival in range(1, M+1):
        m = 4.0*ival # retrieve power
        scaling = (m-1)*(m-3)/(m*(m-2))
        integralVals[ival] = scaling*integralVals[ival-1]

    normVals = np.sqrt(integralVals)
    return normVals


def validateBandwidth(L):
    # perform validation for bandwidth / order L. First ensures it is even, and
    # returns error if modified version is less than 0.
    if (L % 2) != 0:
        L = L - 1

    assert (L >= 0), 'L invalid'

    return L


def makeSteerFilt(filtSize, f, theta, r0, sigmaDenom):
    print(filtSize)
    nn = np.arange(-filtSize, filtSize, 1)
    # nn = np.arange(-filtSize, filtSize+1, 1) #previously
    x, y = np.meshgrid(nn, nn)
    print(x.shape)
    r = np.sqrt(x**2 + y**2)
    sigma = filtSize/sigmaDenom
    g = np.ones(x.shape)
    # g = np.exp(-((r-r0)**2) / (2*sigma**2))
    angle = np.arctan2(x, y)
    angle = np.multiply(angle, (angle >= 0)) - np.multiply(angle, (angle < 0))
    angleSize = angle.size
    a = np.reshape(np.transpose(angle), angleSize)
    values_fun = sp.interpolate.interp1d(theta, f, fill_value="extrapolate")  # interp1(theta,f,a)
    values = values_fun(a)

    filt = np.multiply(g, np.reshape(values, g.shape))

    return filt


def index2boolean(indexArray, N):
    # first verify positive integers
    M = 2 * N + 1

    assert (max(abs(indexArray[:])) <= N), 'Error: max(abs(indexArray(:))) <= N'

    id_ = indexArray + N

    boolArray = np.zeros(M)
    boolArray[id_] = 1
    return boolArray


def validateNonzeroCoefficients(nonzeroCoefficients, N):

    # assert(nonzeroCoefficients.size == 2*N+1)

    # firstHalf = nonzeroCoefficients[0:N]
    # secondHalf = nonzeroCoefficients(N+2:end)

    # assert(isequal(firstHalf(:), flipud(secondHalf(:))))

    numNonzero = np.sum(nonzeroCoefficients)

    return numNonzero


def computeSteerBasis2D(N, nonzeroCoeff):
    # % assumes nonzeroCoeff is a vector of bools
    numNonzero = validateNonzeroCoefficients(nonzeroCoeff, N)
    powersBase = np.transpose(np.arange(-N,N+1,1))
    powerSpec = powersBase[nonzeroCoeff == 1]
    steerThetaBase = getSteerAngles2D(numNonzero) #np.matlib.repmat(steerTheta, numNonzero, 1)
    steerThetaBase = np.matlib.repmat(steerThetaBase, int(numNonzero), 1)
    aux = np.multiply(np.reshape(powerSpec, (len(powerSpec), 1)), steerThetaBase)
    steerBasisVals = np.exp(1j*aux)

    return steerBasisVals


def getSteerAngles2D(M):
    step = np.pi/M
    angles = np.arange(0, (M-1)*step + step, step)

    return angles


def makeSteerBasis(filt, M):
    angles = getSteerAngles2D(M)
    anglesDeg = angles * 180 / np.pi

    dims = filt.shape
    nrows = dims[0]
    ncols = dims[1]
    nfilts = int(M)
    steerFiltVol = np.zeros((nrows, ncols, nfilts), dtype=float)
    for id in range(0, nfilts):
        angleRotate = anglesDeg[id]
        thisFilt = rotate(filt, -angleRotate)
        steerFiltVol[:,:, id] = np.transpose(thisFilt)

    return steerFiltVol


def computeSteerDirection2D(N, Theta, *args):
    varargin = args
    nargin = 2 + len(varargin)
    if nargin < 3:
        nonzeroCoeff = np.ones((2 * N + 1, 1))
    else:
        nonzeroCoeff = args[0]
    # numNonzero = validateNonzeroCoefficients(nonzeroCoeff, N)
    powersBase = np.arange(-N, N+1, 1)
    aux = nonzeroCoeff.astype(int) == 1
    powerSpec = powersBase[aux]
    dirVals = np.exp(1j * (np.multiply(Theta, powerSpec) ) )

    return dirVals

