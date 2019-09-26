#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

# @author: jlvilas
"""
import testFunctions as TstFun
import numpy as np
import matplotlib.pyplot as plt
import random
import _tkinter

xdim = 1000
ydim = 1000
wavelength = 20
angle = 0
mu = 0
sigma = 0.5

# Defining a fringe pattern
img = TstFun.define_sinusoidal_pattern(wavelength, angle, xdim, ydim)

# adding noise to the fringe pattern
imgNoise = TstFun.add_gaussian_noise(img, mu, sigma)

imgplot = plt.imshow(imgNoise)
plt.show(imgplot)

