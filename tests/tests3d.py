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
sys.path.insert(0, '../util')
import steering3D as steer

class TestPolynomialMethods(unittest.TestCase):
    def test_steerable_polynomial3d(self):
        pass
    
    def test_make_default_nonzeroBool_raises_assertion(self): 
        # makes sure error is raised if 0 is passed
        self.assertRaises(AssertionError, steer.make_default_nonzeroBool, 0)

    def test_make_default_nonzeroBool_returns_correct_even_length(self):
        N = 4
        expected = [1, 0, 1, 0, 1]
        results = steer.make_default_nonzeroBool(N)
        self.assertEqual(expected, results)
    
    def test_make_default_nonzeroBool_returns_correct_even_length(self):
        N = 5
        expected = [0, 1, 0, 1, 0, 1]
        results = steer.make_default_nonzeroBool(N)
        self.assertEqual(expected, results)        
    def test_get_basis_returns_correct_phi(self):
        N = 4
        nonzeroBool = steer.make_default_nonzeroBool(N)
        num_samples = 4
        bCos, sqrtSin, phi, dphi = steer.get_basis(N, num_samples, nonzeroBool)
        expected = np.matrix([0, math.pi/3, 2*math.pi/3, math.pi])      
        np.testing.assert_array_equal(expected, phi)                           
    
class TestMakeFilterMethods(unittest.TestCase):
    def test_something(self):
        pass
class TestSteeringMethods(unittest.TestCase):    
    def test_something(self):
        pass
    
if __name__ == '__main__':
    unittest.main()