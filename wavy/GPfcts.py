import numpy as np
from numpy.linalg import inv
from numpy.linalg import cholesky, det, lstsq
import scipy as sp

def kernel(X1, X2, l=1.0, sigma_f=1.0):
    '''
    Isotropic squared exponential kernel. Computes
    a covariance matrix from points in X1 and X2.
    Args:
        X1: Array of m points (m x d).
        X2: Array of n points (n x d).
    Returns:
        Covariance matrix (m x n).
    '''
    sqdist = np.sum(X1**2, 1).reshape(-1, 1) \
            + np.sum(X2**2, 1) - 2 * np.dot(X1, X2.T)
    M = sigma_f**2 * np.exp(-0.5 / l**2 * sqdist)
    return M

def posterior_predictive_nigp(X_s, X_train, Y_train,
l=None,sigma_f=None, sigma_y=None, sigma_x=None,
Grad_fmean=None):
    '''
    Computes statistics of the GP posterior predictive distribution
    from m training data X_train and Y_train and n new inputs X_s.
    Args:
        X_s: New input locations (n x d)
        X_train: Training locations (m x d)
        Y_train: Training targets (m x 1)
        l: length scale parameter
        sigma_f: signal variance parameter
        sigma_y: noise paramter on y
        sigma_x: noise parameter on x
    Returns:
        Posterior mean vector (n x d) and covariance matrix (n x n)
    '''
    # Gradient at training points
    Grad_fmean_train = np.interp(
            X_train.ravel(), X_s.ravel(), Grad_fmean.ravel())
    Grad_fmean_train = Grad_fmean_train.reshape(-1,1)
    Grad_fmean_input = Grad_fmean.copy()
    Grad_fmean_input = Grad_fmean_input.reshape(-1,1)
    # constituents of equations for mean and covariance
    K = kernel(X_train, X_train, l, sigma_f) + \
        sigma_y**2 * np.eye(len(X_train)) + \
        np.diag(np.diag(np.dot(Grad_fmean_train,
            np.atleast_2d(Grad_fmean_train.ravel())))) * sigma_x**2
    K_s = kernel(X_train, X_s, l, sigma_f)
    K_ss = kernel(X_s, X_s, l, sigma_f) + \
        sigma_y**2 * np.eye(len(X_s)) + \
        np.diag(np.diag(np.dot(Grad_fmean_input,
            np.atleast_2d(Grad_fmean_input.ravel())))) * sigma_x**2
    ### possible solutions:
    ## 1. classical
    K_inv = inv(K)
    # Equation (4)
    mu_s = K_s.T.dot(K_inv).dot(Y_train)
    # Equation (5)
    cov_s = K_ss - K_s.T.dot(K_inv).dot(K_s)
#    ## 2. use of cholesky decomposition
#    L = sp.linalg.cholesky(K, lower=True)
#    a1 = np.linalg.lstsq(L, Y_train, rcond=None)[0]
#    alpha = np.linalg.lstsq(L.T, a1, rcond=None)[0]
#    # Equation (4) of Krasser
#    mu_s = K_s.T.dot(alpha)
#    # Equation (5) of Krasser
#    v = np.linalg.lstsq(L, K_s, rcond=None)[0]
#    cov_s = K_ss - v.T.dot(v)
    return mu_s, cov_s

def nll_fn_nigp(X_train,Y_train,Grad_fmean,naive=False):
    '''
    Returns a function that computes the negative log marginal
    likelihood for training data X_train and Y_train and given
    noise level.
    Args:
        X_train: training locations (m x d).
        Y_train: training targets (m x 1).
        noise: known noise level of Y_train.
        naive: if True use a naive implementation of Eq. (7), if
               False use a numerically more stable implementation.
    Returns:
        Minimization objective.
    '''
    def nll_naive(theta):
        # Naive implementation of Eq. (7). Works well for the examples
        # in this article but is numerically less stable compared to
        # the implementation in nll_stable below.
        K = kernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + \
            theta[2]**2 * np.eye(len(X_train)) + \
            np.diag(np.diag(np.dot(Grad_fmean,
                np.atleast_2d(Grad_fmean.ravel())))) * theta[3]**2
        F = 0.5 * np.log(det(K)) + \
               0.5 * Y_train.T.dot(inv(K).dot(Y_train)) + \
               0.5 * len(X_train) * np.log(2*np.pi)
        return F.ravel()
    def nll_stable(theta):
        # Numerically more stable implementation of Eq. (7) as described
        # in http://www.gaussianprocess.org/gpml/chapters/RW2.pdf, Section
        # 2.2, Algorithm 2.1.
        K = kernel(X_train, X_train, l=theta[0], sigma_f=theta[1]) + \
            theta[2]**2 * np.eye(len(X_train)) + \
            np.diag(np.diag(np.dot(Grad_fmean,
                np.atleast_2d(Grad_fmean.ravel())))) * theta[3]**2
        #L = cholesky(K)
        L = sp.linalg.cholesky(K,lower=True)
        F = np.sum(np.log(np.diagonal(L))) + \
                  0.5 * Y_train.T.dot(lstsq(
                  L.T, lstsq(
                    L, Y_train,rcond=None
                    )[0],rcond=None)[0]
                  ) + \
                  0.5 * len(X_train) * np.log(2*np.pi)
        return F.ravel()
    if naive:
        return nll_naive
    else:
        return nll_stable
