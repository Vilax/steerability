# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 16:03:22 2019

@author: tommy
"""
import numpy as np
from numpy import linalg as npla
import util.steering3D as ste


def make_filter(n, r0, sigma_denom, f, phi):
    filter_range = np.arange(-n, n+1, 1)
    # filter_range = np.array(np.arange(-n,n+1))
    print('filter_range.shape = ', filter_range.shape)
    X, Y, Z = np.meshgrid(filter_range, filter_range, filter_range)
    R = np.sqrt(X**2 + Y**2 + Z**2)
    sigma = n / sigma_denom
    
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

def makeOrientedFilter3d(n, f, phi, orientation, radial_type, *args):
    filter_range = np.arange(-n, n+1, 1)
    X, Y, Z = np.meshgrid(filter_range, filter_range, filter_range)
    R = np.sqrt(X**2 + Y**2 + Z**2)
    
    north_pole = [0,0,1]
    orientation = orientation/npla.norm(orientation)
    rotMat = ste.get_rotation_matrix(orientation, north_pole)
    
    num_coordinates = np.size(R)
    coordinates = np.concatenate((X.reshape(num_coordinates,1),\
                            Y.reshape(num_coordinates,1),\
                            Z.reshape(num_coordinates,1)), axis = 1)
    
    rotCoordinates = np.dot(rotMat,coordinates.transpose()).transpose()
    rotCoordinates = np.array(rotCoordinates)
    xrot = np.reshape(rotCoordinates[:,0], X.shape)
    yrot = np.reshape(rotCoordinates[:,1], Y.shape)
    zrot = np.reshape(rotCoordinates[:,2], Z.shape)
    u, v, w = project_to_sphere(X, Y, Z, R)
    
    angles = np.arccos(w)
    a = np.reshape(angles, (np.size(angles), 1))
    phi = np.reshape(phi, np.size(f))
    values = np.interp( np.ndarray.flatten(a), phi, \
                       np.ndarray.flatten(np.array(f)))
    spherical_vals = np.reshape(values, angles.shape)
    
    params = args
    g = get_radial_function(n, radial_type, params)
    filt = np.multiply(g, spherical_vals)
    return filt
    
def get_radial_function(n,radial_type, params):
    params = params[0]
    if radial_type is 'gaussian':
        r0 = params[0]
        sigma = params[1]
        g = makeGaussianRadial3d(n,r0,sigma)
    elif radial_type is 'spline':
        pass
    return g

def makeGaussianRadial3d(n,r0,sigma):
    filter_range = np.arange(-n, n+1, 1)
    X, Y, Z = np.meshgrid(filter_range, filter_range, filter_range)
    R = np.sqrt(X**2 + Y**2 + Z**2)  
    g = np.exp(-np.power(R-r0,2) / (2*(np.power(sigma,2))))
    return g

def project_to_sphere(X, Y, Z, R):
    SMALL_CONSTANT = 1e-8
    R = R+SMALL_CONSTANT
    
    u = np.divide(X, R)
    v = np.divide(Y, R)
    w = np.divide(Z, R)
    
    return u, v, w