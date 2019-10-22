import numpy as np

def compute_steer_basis_angles3d(N, direction_cosines):
    M = (N+1)*(N+2)/2
    direction_cosines = np.matrix(direction_cosines)

    powers = direction_cosine_powers3d(N)
    alpha_powers = powers[:,0]
    beta_powers = powers[:,1]
    gamma_powers = powers[:,2]

    steer_alpha_base = np.tile(direction_cosines[:,0].transpose(), (M,1))
    steer_beta_base = np.tile(direction_cosines[:,1].transpose(), (M,1))
    steer_gamma_base = np.tile(direction_cosines[:,2].transpose(), (M,1))

    steer_alpha = np.powers(steer_alpha_base, alpha_powers)
    steer_beta = np.powers(steer_beta_base, beta_powers)
    steer_gamma = np.powers(steer_gamma_base, gamma_powers)

    steer_angles = np.multiply(np.multiply(steer_alpha, steer_beta), \
                                                        steer_gamma)

    return steer_angles

def direction_cosine_powers3d(N):
    M = (N+1)*(N+2)/2
    powers = np.zeros((M, 3))
    counter = 0
    for ialpha in range(N+1):
      remainder = N - ialpha
      for ibeta in range(remainder+1):
        igamma = N - ialpha - ibeta
        powers[counter, 0] = ialpha
        powers[counter, 1] = ibeta
        powers[counter, 2] = igamma
        counter += 1
    return powers


def make_steer_basis3d(filt, steer_direction_array):
    nfilt = steer_direction_array.shape[0]
    nrows, ncols, nframes = filt.shape

    steer_filt_hypervolume = np.zeros((nrows,ncols,nfilt))

    for ifilt in arange(nfilt):
        steer_direction = steer_direction_array[ifilt,:]
        steer_filt = rotate_filter3d(filt, steer_direction)
        steer_filter_hypervolume[:,:,:,ifilt] = steer_filt

def rotate_filter3d(filt, target_direction, start_direction = (0, 0, 1)):
    northpole = [0, 0, 1]
    
    assert(norm(startDir) > 0)
    startDir = startDir / norm(startDir)
    