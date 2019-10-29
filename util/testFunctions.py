#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

@author: jlvilas
"""
import numpy as np
import random
import scipy as sp


def define_sinusoidal_pattern(wavelength, xdim, ydim, zdim = 1, direction = None):
    # This functions create an image/volume with dimensions "xdim" rows and "ydim"
    # columns (for iamges) and "zdim" slices for volumes. The choice of create an image
    # or a volume is given by zdim.
    # The image is a sinusoidal wave with wavelength defined by "wavelength".
    # The propagation direction of the wave is along the direction defined
    # by the parameter "direction". In the case of defining an images, direction represent
    # the angle of the fringe pattern measured in counter clockwise respect the x-axis.
    # In the case of volumes, direction must be a vector.

    # Initialize the image
    wave_vector = 2*np.pi/wavelength
    x = np.arange(0, xdim, 1)
    y = np.arange(0, ydim, 1)

    if zdim == 1:
        print('----')
        xx, yy = np.meshgrid(x, y)
        img = np.sin(wave_vector * (xx * np.cos(direction) - yy * np.sin(direction)))
    else:
        z = np.arange(0, zdim, 1)
        xx, yy, zz = np.meshgrid(x, y, z)
        img = np.sin(wave_vector * (xx * direction[0] + yy * direction[1] + zz * direction[2]))

    return img


def define_chirp_pattern(wavelength, angle, xdim, ydim):
    # Initialize the image

    x = np.arange(0, xdim, 1)
    y = np.arange(0, ydim, 1)
    xx, yy = np.meshgrid(x, y)
    sqxx = np.sqrt(xx)
    idx = sqxx<=10
    wave_vector = 2* np.pi / ( sqxx + 1e-38)
    wave_vector[idx] = 2* np.pi / wavelength
    img = np.sin(wave_vector * (xx * np.cos(angle) - yy * np.sin(angle)))

    return img


def add_gaussian_noise(img, mu, sigma):
    # This functions takes an input image "img", and adds to each pixels
    # a random value that follows a gaussian distribution with mean "mu" and
    # standard deviation "sigma". Thus, the output image NoisyImg = img + Noise

    NoisyImg = np.copy(img)
    NoisyImg = NoisyImg + np.random.normal(mu, sigma, img.shape)

    return NoisyImg