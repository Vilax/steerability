import numpy as np
import math
import scipy.linalg as la
import matplotlib.pyplot as plt
from numpy import linalg as npla
#import scipy.sparse.linalg as sla

def steerable_polynomial3d(poleCap, N, nonzeroBool=None):
    
    DEFAULT_SAMPLES = 400
    num_samples = DEFAULT_SAMPLES
    
    if nonzeroBool is None:
        nonzeroBool = make_default_nonzeroBool(N)
    assert(len(nonzeroBool) == N+1)    
    nonzeroBool = np.array(nonzeroBool)        
    
    bCos, sqrtSin, phi, dphi = get_basis(N, num_samples, nonzeroBool)    

    G1 = make_gram_matrix(bCos,sqrtSin,phi,dphi,poleCap)  
    G2 = make_gram_matrix(bCos,sqrtSin,phi,dphi,math.pi/2)
    
    D, V = la.eig(G1, G2)
    eigval_max = get_max_arg_real(D)
    V = np.matrix(V)
    v = V[:,eigval_max]
    v = v / npla.norm(v)
    if v[0] < 0:
        v = v * (-1)
        
    print('eigenvalues are: \n', D)
    print('eigenvector is \n', v)
    constraint = v.transpose() * G2 * v
    obj = (v.transpose() * G1 * v) / constraint
    print('printing constraint and then objective')
    print(constraint)
    print(obj)
    f = make_steerable_function(v, bCos)
    f = f/f[0]
    plt.figure()
    plt.plot(phi, f)
    
    # now plot f
    return f, v, bCos, phi

def get_max_arg_real(a):
    # for array a, gets argument of max value, with complex values treated as
    # negative infinity
    neg_inf = float('-inf')
    a[np.iscomplex(a)] = neg_inf
    return np.argmax(a)
    

def make_steerable_function(v, bCos):
    nVals = bCos.shape[0]
    npowers = bCos.shape[1]
    f = np.matrix(np.zeros((nVals, 1)))
    for j in range(npowers):
        f += bCos[:,j] * v[j]
    return f        

def make_default_nonzeroBool(N):
    assert(N>0)
    nonzeroBool = [0] *(N+1)
    num_default_ones = math.ceil(len(nonzeroBool)/2)
    nonzeroBool[-1::-2] = [1]*num_default_ones
    return nonzeroBool

def get_basis(N, num_samples, nonzeroBool):
    assert(num_samples > 1)
    dphi = math.pi/(num_samples-1)
    phi = np.arange(0, math.pi+dphi,dphi)
    bCos = np.matrix(np.zeros((len(phi),N+1)))
    cosPhi = np.matrix(np.cos(phi)).transpose()
    sqrtSin = np.matrix(np.sqrt(np.sin(phi))).transpose()
    for j in range(N+1):
        bCos[:,j] = np.power(cosPhi, j)
    sinMat = np.tile(sqrtSin, [1, bCos.shape[1]])    
    bCosSin = np.multiply(bCos, sinMat)
    bCosNorm = np.matrix(np.sqrt(np.diagonal(bCosSin.transpose() * \
                                             bCosSin * dphi)))
    bCosNorm = np.tile(bCosNorm, [bCos.shape[0], 1])
    bCos = np.divide(bCos, bCosNorm)
    bCos = bCos[:, nonzeroBool == [1]*len(nonzeroBool)]

    return bCos, sqrtSin, phi, dphi

def make_gram_matrix(bCos, sqrtSin, phi, dphi, phiEnd = math.pi/2):
    bCos = bCos[phi<phiEnd, :]
    sqrtSin = sqrtSin[phi<phiEnd]
    sinMat = np.tile(sqrtSin, [1, bCos.shape[1]])    
    bCos = np.multiply(bCos, sinMat)
    G = bCos.transpose() * bCos * dphi
    return G

