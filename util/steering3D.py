import numpy as np
import math
import scipy as sp
from scipy.spatial.transform import Rotation as R
from numpy import linalg as la
from scipy import ndimage as ndimg
import util.polynomial3d as poly3d
import util.filter3d as filt
import util.utilFuntions as uf


def directionalFilter3D(N, cap, r0, sigma, direction, filtSize):
    # This function gives volume with the directional filter
    # N represent the order of the polynomial of the filter
    # The angle of the filter is defined as pi/cap
    # sigma, si the standard deviation of the gaussian which defines the frequencies
    # direction in radians
    import time

    f, v, bCos, phi = poly3d.steerable_polynomial3d(cap, N)
    filtSize = int(filtSize) / 2  # in pixels
    start = time.time()
    orig_filt = filt.make_filter(filtSize, r0, sigma, f, phi)
    end = time.time()
    print('Total Time ', end - start)
    img = orig_filt[:, 100, :]
    uf.representImg(img, 'Filter', True)

    direction_cosines = np.loadtxt('3d_direction_cosines_N4.txt')
    start = time.time()
    steer_basis_hypervol = make_steer_basis3d(orig_filt, direction_cosines)
    end = time.time()
    print('Total Time ', end - start)
    num_basis_filters = steer_basis_hypervol.shape[3]

    start = time.time()
    steer_basis_angles = compute_steer_basis_angles3d(N, direction_cosines)
    end = time.time()
    print('Total Time ', end - start)
    start = time.time()
    Phi = np.linalg.inv(steer_basis_angles)
    end = time.time()
    print('Total Time ', end - start)
    sym_test_axes = np.array([[1, 2, 3], [3, 4, 5]])
    nAxes = sym_test_axes.shape[0]
    iaxis = sym_test_axes[0, :]
    axis_sym = sym_test_axes[1, :]
    axis_sym = axis_sym / np.linalg.norm(axis_sym)
    alpha, beta, gamma = axis_sym

    axis_sym = direction/np.linalg.norm(axis_sym)
    start = time.time()
    direction_angles = compute_direction_angle_powers(N, axis_sym)
    end = time.time()
    print('Total Time ', end - start)
    k = Phi * direction_angles
    k = np.ndarray.flatten(np.array(k))

    steeredFilt = np.abs(np.dot(steer_basis_hypervol, k))
    steeredFilt = steeredFilt / np.max(steeredFilt)

    uf.representImg(steeredFilt[:,:,100], 'testDir', True)

    angleCritic, fillingValue = estimateFilterWidth3D(steeredFilt, direction)
    print("angleCritic =", angleCritic*180/np.pi, " ", fillingValue)
    steeredFilt = maskrippling3D(steeredFilt, direction, filtSize, angleCritic, fillingValue)
    mask = uf.create_spherical_mask(steeredFilt.shape[0], steeredFilt.shape[1], steeredFilt.shape[2])

    steeredFilt = np.multiply(steeredFilt, mask)

    return steeredFilt


def rotationAroundAxis(axisVector, angle):
    # This function implements the Rodrigues formula
    ux = axisVector[0]
    uy = axisVector[1]
    uz = axisVector[2]

    m11 = np.cos(angle) + (1 - np.cos(angle)) * ux ** 2
    m12 = ux*uy*(1 - np.cos(angle)) - uz * np.sin(angle)
    m13 = ux*uz*(1 - np.cos(angle)) + uy * np.sin(angle)
    m21 = ux*uy*(1 - np.cos(angle)) + uz * np.sin(angle)
    m22 = np.cos(angle) + (1 - np.cos(angle)) * uy ** 2
    m23 = uy*uz*(1 - np.cos(angle)) - ux * np.sin(angle)
    m31 = ux*uz*(1 - np.cos(angle)) - uy * np.sin(angle)
    m32 = uy*uz*(1 - np.cos(angle)) + ux * np.sin(angle)
    m33 = np.cos(angle) + (1 - np.cos(angle)) * uz ** 2

    rotation = [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]

    return rotation


def maskrippling3D(steeredFilt, direction, filtSize, angleCritic, value):
    # Directional filters based on steerability usualy present a rippling
    # this function mask that rippling, resulting in a monotonic and
    # smooth directional filter.

    direction = direction/np.linalg.norm(direction)
    nn = np.arange(-filtSize, filtSize, 1)
    x, y, z = np.meshgrid(nn, nn, nn)
    r = np.sqrt(x**2 + y**2 + z**2)

    x = np.arccos(np.true_divide( np.abs(np.multiply(y, direction[0]) + np.multiply(x, direction[1]) + np.multiply(z, direction[2])), r))

    # uf.representImg(x[:,:,100],'angles', True)
    idx = x > angleCritic

    steeredFilt[idx] = value

    return steeredFilt


def estimateFilterWidth3D(filter, direction):
    # Center of the image
    center_dir = int(np.floor(0.5 * filter.shape[0]))

    ux = direction[0]
    uy = direction[1]
    uz = direction[2]
    direction = direction / np.linalg.norm(direction)

    ind = np.argmax([np.arccos(ux), np.arccos(uy), np.arccos(uz)])

    eyeMat = np.eye(3)
    print(eyeMat)
    axisVector = np.cross(eyeMat[:,ind], direction)
    axisVector = axisVector/ np.linalg.norm(direction)

    r = (center_dir * 2/3) * direction

    last_idx_x = np.int(center_dir + r[0])
    last_idx_y = np.int(center_dir + r[1])
    last_idx_z = np.int(center_dir + r[2])
    print("Init idx_x ", last_idx_x, "  idx_y ", last_idx_y, "  idx_z ", last_idx_z)

    angleCritic = 0
    value = 0
    valuetest = np.array([])
    ran = np.arange(0, np.pi / 2, np.pi / 180)
    lastAngle = ran[0]
    for theta in ran:
        rotMatrix = rotationAroundAxis(axisVector, theta)
        r_rotated = np.dot(rotMatrix, r)

        idx_x = np.int(center_dir + r_rotated[0])
        idx_y = np.int(center_dir + r_rotated[1])
        idx_z = np.int(center_dir + r_rotated[2])
        # print("idx_x ", idx_x, "  idx_y ", idx_y, "  idx_z ", idx_z)

        valuetest = np.append(valuetest, [filter[idx_x, idx_y, idx_z]])

        if filter[idx_x, idx_y, idx_z] > filter[last_idx_x, last_idx_y, last_idx_z]:
            angleCritic = lastAngle
            value = filter[last_idx_x, last_idx_y, last_idx_z]
            break
        lastAngle = theta
        last_idx_x = idx_x
        last_idx_y = idx_y
        last_idx_z = idx_z

    uf.representCurve(np.arange(0,len(valuetest),1), valuetest, 'filter_size', True)
    return angleCritic, value


def compute_steer_basis_angles3d(N, direction_cosines):

    direction_cosines = np.matrix(direction_cosines)
    steer_angles = compute_direction_angle_powers(N, direction_cosines)

    return steer_angles


def compute_direction_angle_powers(N, directions):
    M = int((N+1)*(N+2)/2)
    directions = np.matrix(directions)
    
    powers = np.matrix(direction_cosine_powers3d(N))
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
    
    