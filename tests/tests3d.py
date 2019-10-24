# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 12:33:23 2019

@author: tommy
"""

import unittest
import sys
import pandas as pd
import numpy as np
import math
sys.path.insert(0, 'E:\\steerability\\util')
sys.path.insert(0, 'E:\\steerability\\scripts')
import steering3D as steer
import filter3d as filt
import polynomial3d as poly3d
import scipy.io as spio
from numpy import linalg as npla

class TestPolynomialMethods(unittest.TestCase):
    N = 4
    cap_size = math.pi/3    
    f, v, bCos, phi = poly3d.steerable_polynomial3d(cap_size, N)
    def test_steerable_polynomial3d(self):
        mat = spio.loadmat('f_N4.mat')
        f_matlab = mat['f']
        results = np.array(self.f)
        expected = f_matlab
        np.testing.assert_array_almost_equal(expected, results, decimal = 4)
    
    def test_make_default_nonzeroBool_raises_assertion(self): 
        # makes sure error is raised if 0 is passed
        self.assertRaises(AssertionError, poly3d.make_default_nonzeroBool, 0)

    def test_make_default_nonzeroBool_returns_correct_even_length(self):
        N = 4
        expected = [1, 0, 1, 0, 1]
        results = poly3d.make_default_nonzeroBool(N)
        self.assertEqual(expected, results)
    
    def test_make_default_nonzeroBool_returns_correct_even_length(self):
        N = 5
        expected = [0, 1, 0, 1, 0, 1]
        results = poly3d.make_default_nonzeroBool(N)
        self.assertEqual(expected, results)        
    def test_get_basis_returns_correct_phi(self):
        N = 4
        nonzeroBool = poly3d.make_default_nonzeroBool(N)
        num_samples = 4
        bCos, sqrtSin, phi, dphi = poly3d.get_basis(N, num_samples, nonzeroBool)
        expected = np.array([0, math.pi/3, 2*math.pi/3, math.pi])      
        np.testing.assert_array_equal(expected, phi)    
    def test_get_basis_returns_correct_bcos(self):
        results = np.array(self.bCos)
        mat = spio.loadmat('bCos_N4.mat')
        expected = mat['bCos']
        np.testing.assert_allclose(expected, results, rtol=0.05)      
    
class TestMakeFilterMethods(unittest.TestCase):
    N = 4
    cap_size = math.pi/3
    f, v, bCos, phi = poly3d.steerable_polynomial3d(cap_size, N)
    n = 100
    r0 = 20
    sigma_denom = 5
    orig_filt = filt.make_filter(n, r0, sigma_denom, f, phi)
    
    def test_makes_correct_filter(self):
        mat = spio.loadmat('orig_filt_matlab.mat')
        filt_matlab = mat['filt']
        expected = filt_matlab
        results = self.orig_filt
        np.testing.assert_array_almost_equal(expected, results)
    
    
class TestSteeringMethods(unittest.TestCase): 
    direction_cosines = np.loadtxt('3d_direction_cosines_N4.txt')
    N = 4
    def test_steer_basis_angles_makes_correct_matrix(self):
        steer_basis_angles = steer.compute_steer_basis_angles3d(self.N, \
                                                        self.direction_cosines)
        # test determinant
        results = npla.det(steer_basis_angles)
        expected = -0.0000014387
        self.assertAlmostEqual(expected, results)    

    
if __name__ == '__main__':
    unittest.main()