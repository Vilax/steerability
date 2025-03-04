# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 16:03:22 2019

@author: tommy
"""
import numpy as np


def make_filter(filtSize, r0, sigma_denom, f, phi):
    filter_range = np.arange(-filtSize, filtSize, 1)
    # filter_range = np.array(np.arange(-n,n+1))
    print('filter_range.shape = ', filter_range.shape)
    X, Y, Z = np.meshgrid(filter_range, filter_range, filter_range)
    R = np.sqrt(X**2 + Y**2 + Z**2)
    sigma = filtSize / sigma_denom
    
    # g = np.exp(-np.power(R-r0,2) / (2*(np.power(sigma,2))))
    g = np.ones(X.shape)
    import time
    start = time.time()
    u, v, w = project_to_sphere(X, Y, Z, R)
    end = time.time()
    print("Time", end-start)
    
    angles = np.arccos(w)
    a = np.reshape(angles, (np.size(angles), 1))
    phi = np.reshape(phi, np.size(f))
    values = np.interp( np.ndarray.flatten(a), phi, \
                       np.ndarray.flatten(np.array(f)))
    spherical_vals = np.reshape(values, angles.shape)
    
    filt = np.multiply(g, spherical_vals)
    return filt

    
def project_to_sphere(X, Y, Z, R):
    SMALL_CONSTANT = 1e-8
    R = R+SMALL_CONSTANT
    
    u = np.divide(X, R)
    v = np.divide(Y, R)
    w = np.divide(Z, R)
    
    return u, v, w