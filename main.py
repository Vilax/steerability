#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

# @author: jlvilas
"""
import testFunctions as TstFun
import numpy as np
import matplotlib.pyplot as plt
import tommyFunctions as tomFun
import random
import _tkinter


# xdim = 1000
# ydim = 1000
# wavelength = 20
# angle = 0
# mu = 0
# sigma = 0.5
#
# # Defining a fringe pattern
# img = TstFun.define_sinusoidal_pattern(wavelength, angle, xdim, ydim)
#
# # adding noise to the fringe pattern
# imgNoise = TstFun.add_gaussian_noise(img, mu, sigma)
#
# # imgplot = plt.imshow(imgNoise)
# # plt.show(imgplot)

#######################
capdenomArr = [3, 4, 5, 6, 8, 10, 20]


#colorMat = tomFun.distinguishable_colors(len(capdenomArr))

L = 4
allRoots4 = []

for capdenom in capdenomArr:
    print(capdenom)
    theta = np.pi/capdenom
    print("theta %f", theta*180/np.pi)
    a = tomFun.concentratedPolyPoleCaps(theta, L)
    print('checkpoint 1')
    print(np.polyval([a],1.))
    p = np.true_divide(a, np.sum(a))  #According with Tommy np.true_divide(a, np.polyval(a,1.))
    print('checkpoint 2')
    # aa = np.concatenate( (p, np.zeros(p.shape) ))

    bb = [p, np.asmatrix(np.zeros(p.shape)) ]

    aux = np.transpose(bb)

    print('checkpoint 3')
    print(np.size(p))
    p = np.reshape(aux, 2*np.size(p) )
    p = np.squeeze(p)
    print('checkpoint 4')
    p = p[:-1]
    rts = np.roots(p)

    # cc = [allRoots4, rts]
    # allRoots4 = np.concatenate(allRoots4, rts)
    # #color = colorMat(capId, :)
    realpart = np.real(rts)
    imagpart = np.imag(rts)
    print("real", realpart, "imag", imagpart)
    # plt.plot(rts, '.')
    # plt.show()

#TODO: Add legend (Matlab legend;)

