#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

@author: jlvilas
"""
import numpy as np


def paddingImageIfIsOdd(img):
    # The input image/volume is padded to have odd dimensions along all axis

    dims = len(img.shape)
    if dims == 2:
        if (img.shape[0] % 2) == 0:
            imgPadded = np.pad(img, ((0, 1), (0, 0)), 'edge')
        else:
            imgPadded = img
        if (img.shape[1] % 2) == 0:
            imgPadded = np.pad(imgPadded, ((0, 0), (0, 1)), 'edge')
    if dims == 3:
        if (img.shape[0] % 2) == 0:
            imgPadded = np.pad(img, ((0, 1), (0, 0), (0, 0)), 'edge')
            # imgPadded = np.pad(img, ((0, 1), (0, 0), (0, 0)), 'constant', constant_values=(0,))
            # imgPadded[-1, :-1] = lastrow
        else:
            imgPadded = img

        if (img.shape[1] % 2) == 0:
            # lastcol = imgPadded[:, -1]
            # lastrow = imgPadded[-1, :]
            imgPadded = np.pad(imgPadded, ((0, 0), (0, 1), (0, 0)), 'edge')
            # imgPadded = np.pad(imgPadded, ((0, 0), (0, 1), (0, 0)), 'constant', constant_values=(0,))
            # imgPadded[:-1, -1] = lastcol
            # imgPadded[-1, -1] = 0.5 * (lastrow[-1] + lastcol[-1])

        if (img.shape[2] % 2) == 0:
            imgPadded = np.pad(imgPadded, ((0, 0), (0, 0), (0, 1)), 'edge')
            # lastcol = imgPadded[:, -1]
            # lastrow = imgPadded[-1, :]
            # imgPadded = np.pad(imgPadded, ((0, 0), (0, 0), (0, 1)), 'constant', constant_values=(0,))
            # imgPadded[:-1, -1] = lastcol
            # imgPadded[-1, -1] = 0.5 * (lastrow[-1] + lastcol[-1])
        print("imgPadded = ", imgPadded.shape)
    return imgPadded


def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask


def createRing(height, width, rad):
    mask1 = create_circular_mask(height, width, radius=rad)
    mask2 = create_circular_mask(height, width, radius=rad+1)

    ring = mask1 ^ mask2

    return ring


def createRingHPF(height, width, rad):
    dim = min(height, width)
    mask1 = create_circular_mask(height, width, radius=rad)
    mask2 = create_circular_mask(height, width, radius=dim)

    ring = mask1 ^ mask2

    return ring


def createGaussianHPF(filterSize, sigma, r0):
    nn = np.arange(-filterSize, filterSize, 1)
    x, y = np.meshgrid(nn, nn)
    r = np.sqrt(x ** 2 + y ** 2)
    sigma = 0.368*r0
    gaussianImg = np.exp(-((r-r0)**2) / (2*sigma**2))
    gaussianImg[r>r0] = 1.0
    # representImg(gaussianImg, 'gauss', True)

    return gaussianImg


def createGaussian(filterSize, sigma, r0):
    nn = np.arange(-filterSize, filterSize, 1)
    x, y = np.meshgrid(nn, nn)
    r = np.sqrt(x ** 2 + y ** 2)

    gaussianImg = np.exp(-((r-r0)**2) / (2*sigma**2))

    return gaussianImg


def create_spherical_mask(xdim, ydim, center=None, radius=None):

    if center is None: # use the middle of the image
        center = [int(xdim/2), int(ydim/2)]#, int(zdim/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], center[2], xdim-center[0], ydim-center[1]) #, zdim-center[2])
    print(center)
    Y, X = np.ogrid[:xdim, :ydim]#, :zdim]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2 )#+)# (Z - center[2])**2 )

    mask = dist_from_center <= radius

    return mask

def createShell(xdim, ydim, zdim, rad):
    mask1 = create_spherical_mask(xdim, ydim, zdim, radius=rad)
    mask2 = create_spherical_mask(xdim, ydim, zdim, radius=rad+1)

    shell = mask1 ^ mask2

    return shell


def createShellHPF(xdim, ydim, zdim, rad):

    dim = min(xdim, ydim, zdim)

    mask1 = create_spherical_mask(xdim, ydim, zdim, radius=rad)
    mask2 = create_spherical_mask(xdim, ydim, zdim, radius=dim)

    ring = mask1 ^ mask2

    return ring


def createGaussian3DHPF(filterSize, sigma, r0):
    nn = np.arange(-filterSize, filterSize, 1)

    x, y, z = np.meshgrid(nn, nn, nn)
    r = np.sqrt(x ** 2 + y ** 2, z ** 2)

    gaussianImg = np.exp(-((r-r0)**2) / (2*sigma**2))
    gaussianImg[r>r0] = 1.0
    # representImg(gaussianImg, 'gauss', True)

    return gaussianImg


def createGaussian3D(filterSize, sigma, r0):
    nn = np.arange(-filterSize, filterSize, 1)
    x, y, z = np.meshgrid(nn, nn, nn)
    r = np.sqrt(x ** 2 + y ** 2, z ** 2)

    gaussianImg = np.exp(-((r-r0)**2) / (2*sigma**2))

    return gaussianImg


def representImg(img, title, blockplot = False):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.imshow(img)
    plt.colorbar()
    plt.title(title)
    plt.show(block = blockplot)

    pass


def representCurve(x, y, title, blockplot = False):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(x, y)
    plt.title(title)
    plt.show(block = blockplot)

    pass
