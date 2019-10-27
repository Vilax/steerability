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
import util.utilFuntions as uf
sys.path.insert(0, '../util')

# testing steering a filter derived from a gaussian of the radius times an even
# polynomial

N = 4
cap_size = np.pi/3.5
r0 = 20
sigma = 5
filtSize = 200

direction = [0, 0, 1]

Dirfilt = steer.directionalFilter3D(N, cap_size, r0, sigma, direction, filtSize)
