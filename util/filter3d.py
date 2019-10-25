# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 16:03:22 2019

@author: tommy
"""
import numpy as np

def make_filter(n, r0, sigma_denom, f, phi):
    
    filter_range = np.array(np.arange(-n,n+1))
    X, Y, Z = np.meshgrid(filter_range, filter_range, filter_range)
    R = np.sqrt(np.power(X,2)+np.power(Y,2)+np.power(Z,2))
    sigma = n / sigma_denom
    
    g = np.exp(-np.power(R-r0,2) / (2*(np.power(sigma,2))))
    
    u, v, w = project_to_sphere(X, Y, Z)
    
    angles = np.arccos(w)
    a = np.reshape(angles, (np.size(angles), 1))
    phi = np.reshape(phi, np.size(f))
    values = np.interp( np.ndarray.flatten(a), phi, \
                       np.ndarray.flatten(np.array(f)))
    spherical_vals = np.reshape(values, angles.shape)
    
    filt = np.multiply(g, spherical_vals)
    return filt
    
    
    
def project_to_sphere(X, Y, Z):
    SMALL_CONSTANT = 1e-8
    r = np.sqrt(np.power(X,2)+np.power(Y,2)+np.power(Z,2))+SMALL_CONSTANT
    
    u = np.divide(X, r)
    v = np.divide(Y, r)
    w = np.divide(Z, r)
    
    return u, v, w