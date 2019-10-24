import numpy as np
import math
import scipy as sp
from scipy.spatial.transform import Rotation as R
from numpy import linalg as la
from scipy.spatial.transform import Rotation as R
from scipy import ndimage as ndimg


def compute_steer_basis_angles3d(N, direction_cosines):

    direction_cosines = np.matrix(direction_cosines)
    steer_angles = compute_direction_angle_powers(N, direction_cosines)

    return steer_angles

def compute_direction_angle_powers(N, directions):
    M = int((N+1)*(N+2)/2)
    directions = np.matrix(directions)
    
    powers = direction_cosine_powers3d(N)
    alpha_powers = powers[:,0]
    beta_powers = powers[:,1]
    gamma_powers = powers[:,2]  
    
    steer_alpha_base = np.tile(directions[:,0].transpose(), (M,1))
    steer_beta_base = np.tile(directions[:,1].transpose(), (M,1))
    steer_gamma_base = np.tile(directions[:,2].transpose(), (M,1))

    steer_alpha = np.power(steer_alpha_base, alpha_powers)
    steer_beta = np.power(steer_beta_base, beta_powers)
    steer_gamma = np.power(steer_gamma_base, gamma_powers)

    directions_powers = np.multiply(np.multiply(steer_alpha,steer_beta), \
                                    steer_gamma)

    return directions_powers    

def direction_cosine_powers3d(N):
    M = int((N+1)*(N+2)/2)
    powers = np.zeros((M, 3))
    counter = 0
    for ialpha in np.arange(N,-1,-1):
      remainder = N - ialpha
      for ibeta in np.arange(remainder,-1,-1):
        igamma = remainder - ibeta
        powers[counter,:] = [ialpha, ibeta, igamma]
        counter += 1
    return powers


def make_steer_basis3d(filt, steer_direction_array):
    nfilt = steer_direction_array.shape[0]
    nrows, ncols, nframes = filt.shape

    steer_filt_hypervolume = np.zeros((nrows,ncols,nframes,nfilt))

    for ifilt in np.arange(nfilt):
        steer_direction = steer_direction_array[ifilt,:]
        steer_filt = rotate_filter3d(filt, steer_direction)
        steer_filt_hypervolume[:,:,:,ifilt] = steer_filt
        
    return steer_filt_hypervolume
        
def rotate_filter3d(filt, target_direction, start_direction = (0, 0, 1)):
    
    target_direction = np.array(target_direction)
    start_direction = np.array(start_direction)
    
    assert(la.norm(start_direction) > 0)
    start_direction = start_direction / la.norm(start_direction)
    assert(la.norm(target_direction) > 0)
    target_direction = target_direction / la.norm(target_direction)
    
    rotation_matrix = get_rotation_matrix(start_direction,\
                                                    target_direction)
    
    z1, x1, z2 = get_euler_angles(rotation_matrix)
    # apply python rotations
    rotated_filt = ndimg.rotate(filt, z1, [0,1], reshape = False)
    rotated_filt = ndimg.rotate(rotated_filt, x1, [1,2], reshape = False)
    rotated_filt = ndimg.rotate(rotated_filt, z2, [0,1], reshape = False)
    
    return rotated_filt

def get_rotation_matrix(start_direction, target_direction):
    # first write in axis angle form
    angle = np.arccos(np.dot(target_direction, start_direction))
    w = np.cross(start_direction, target_direction)
    axis = w / la.norm(w)
    u1, u2, u3 = axis
    
    cosTheta = np.cos(angle)
    sinTheta = np.sin(angle)
    
    rotation_matrix = np.matrix([[cosTheta + (u1**2)*(1-cosTheta), u1*u2*\
                               (1-cosTheta)-u3*sinTheta, u1*u3*(1-cosTheta) + \
                               u2*sinTheta],\
                                 [u2*u1*(1-cosTheta)+u3*sinTheta, cosTheta + \
                                  (u2**2)*(1-cosTheta), u2*u3*(1-cosTheta)-u1\
                                  *sinTheta],\
                                 [u3*u1*(1-cosTheta)-u2*sinTheta, u3*u2*\
                                  (1-cosTheta) + u1*sinTheta, cosTheta + (u3**2) \
                                  *(1-cosTheta)]])
    
    return rotation_matrix

def get_euler_angles(rotation_matrix):  
    r = R.from_dcm(np.array(rotation_matrix))
    z2, x1, z1 = r.as_euler('zxz', degrees=True)
    return z1, x1, z2
    
    