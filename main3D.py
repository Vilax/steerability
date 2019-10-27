# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:24:03 2019

@author: tommy
"""

import numpy as np
import math
import scipy.linalg as la
import matplotlib.pyplot as plt
# from scipy.spatial.transform import Rotation as R
from numpy import linalg as npla
# from scipy.spatial.transform import Rotation as R
# from scipy import ndimage as ndimg
import sys
import util.steering3D as steer
import util.filter3d as filt
import util.polynomial3d as poly3d
import mrcfile
sys.path.insert(0, '../util')

# testing steering a filter derived from a gaussian of the radius times an even
# polynomial

N = 4
cap_size = math.pi/3.5

f, v, bCos, phi = poly3d.steerable_polynomial3d(cap_size, N)

n = 100
r0 = 20
sigma_denom = 5
orig_filt = filt.make_filter(n, r0, sigma_denom, f, phi)

# view this filter
file = mrcfile.new('orig_filt.mrc', overwrite = True)
file.set_data(orig_filt.astype('float32'))
#with mrcfile.new('tmp.mrc') as mrc:
#    mrc.set_data(orig_filt)

direction_cosines = np.loadtxt('3d_direction_cosines_N4.txt')
steer_basis_hypervol = steer.make_steer_basis3d(orig_filt, direction_cosines)
num_basis_filters = steer_basis_hypervol.shape[3]

steer_basis_angles = steer.compute_steer_basis_angles3d(N, direction_cosines)
Phi = npla.inv(steer_basis_angles)

sym_test_axes = np.array([[1, 2, 3],[3,4,5]])
nAxes = sym_test_axes.shape[0]

for iaxis, axis_sym in enumerate(sym_test_axes):
    axis_sym = axis_sym / npla.norm(axis_sym)
    alpha, beta, gamma = axis_sym
    
    direction_angles = steer.compute_direction_angle_powers(N, axis_sym)
    
    k = Phi * direction_angles
    k = np.ndarray.flatten(np.array(k))
    
    steered_filter = np.dot(steer_basis_hypervol, k)
    brute_force_filter = steer.rotate_filter3d(orig_filt, axis_sym)
    
    diff = steered_filter - brute_force_filter
    normalized_diff = npla.norm(diff) / npla.norm(steered_filter)
    print(normalized_diff)
    assert(normalized_diff < 0.005)
    
    
    