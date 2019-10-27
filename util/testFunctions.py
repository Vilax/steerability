#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

@author: jlvilas
"""
import numpy as np
import random
import scipy as sp


def define_sinusoidal_pattern(wavelength, angle, xdim, ydim):
    # This functions create an image with dimensions "xdim" rows and "ydim"
    # columns. The image is a sinusoidal wave with wavelength, wavelength.
    # The propagation direction of the wave is along the direction defined
    # by the angle "angle" measured in counter clockwise respect the x-axis

    # Initialize the image
    img = np.zeros((xdim, ydim))
    wave_vector = 2*np.pi/wavelength
    x = np.arange(0, xdim, 1)
    y = np.arange(0, ydim, 1)
    xx, yy = np.meshgrid(x, y)
    img = np.sin(wave_vector * (xx * np.cos(angle) - yy * np.sin(angle)))

    return img


def define_chirp_pattern(wavelength, angle, xdim, ydim):
    # Initialize the image
    img = np.zeros((xdim, ydim))

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

    # for i in range(0, NoisyImg.shape[0]):
    #     for j in range(0, NoisyImg.shape[0]):
    #         NoisyImg[i, j] += random.gauss(mu, sigma)
    return NoisyImg