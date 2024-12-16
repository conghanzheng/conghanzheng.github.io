"""
Preliminaries
"""
##  Package Dependencies
import os
import numpy as np
from scipy.io import loadmat, savemat
from scipy.optimize import minimize
import numba

## Set working directory to the script's directory, for saving results
script_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_path)

"""
Load online data
"""

data = loadmat('https://raw.githubusercontent.com/conghanzheng/conghanzheng.github.io/master/assets/Python/cleandata.mat')
EExt = data['EExt']  # Cumulative mileage
EEit = data['EEit']  # Dummy indicator matrix

"""
Transitions
"""

## Empirical frequencies
interv = np.ceil(EExt / 5000).astype(int)  # same as ceil(EExt/5000)
dfint = np.diff(interv, axis=0)  # difference along the time dimension (row-wise)

## Prob Vector
phi = np.zeros(3)
for i in range(3):
    ind = (dfint == (i)).astype(float)
    phi[i] = np.nanmean(ind)

phi = phi / np.sum(phi)

## Building Transition matrix
trans_mat = np.zeros((90, 90))

for i in range(0, 89):
    trans_mat[i, i] = phi[0]    # phi(1,1) in MATLAB is phi[0] in Python
    trans_mat[i, i+1] = phi[1]  # phi(1,2) in MATLAB is phi[1] in Python

for i in range(0, 88):
    trans_mat[i, i+2] = phi[2]  # phi(1,3) in MATLAB is phi[2] in Python

trans_mat[89, 89] = 1
trans_mat[88, 89] = 1 - phi[0]

## Save intermediate results
savemat('results_nfxp.mat', {'phi': phi, 'trans_mat': trans_mat})

## Parameters and Initial Values
beta = 0.99
thetastart = np.array([10, 2])

"""
    Function: Compute the value function by fixed point iteration

    Inputs:
    - beta: discount factor
    - theta: vector of parameters [theta_r, theta_m]
    - trans_mat: 90x90 transition matrix

    Output:
    - ev: 90x1 vector of expected value function
"""
@numba.njit
def valuefunction(beta, theta, trans_mat):
    
    tol = 1e-6
    ev = np.zeros((90, 1))
    rule = 1.0
    while rule > tol:
        x = np.arange(1, 91).reshape(-1, 1) # state: 1 to 90
        c = 0.001 * theta[1] * x   # maintenance cost
        exp1 = np.exp(-theta[0] + beta * ev[0, 0])   # scalar
        exp0 = np.exp(-c + beta * ev)                # 90x1
        # ev_new = trans_mat * log(exp1 + exp0), matrix multiplication and log inside
        # Carefully: log(exp1+exp0) is elementwise. 
        # Need to do matrix multiplication of trans_mat with the vector log(exp1+exp0)
        
        vals = np.log(exp1 + exp0)  # 90x1
        ev_new = trans_mat @ vals   # 90x1
        
        diff = ev_new - ev
        rule = np.sqrt(np.sum(diff**2))
        ev = ev_new
    return ev

"""
    Function: Negative log-likelihood for Rust model.

    Inputs:
    - EExt: T x N matrix of cumulative mileage states
    - EEit: T x N matrix of engine replacement dummies
    - beta: discount factor
    - theta: vector of parameters [theta_r, theta_m]
    - trans_mat: 90x90 transition matrix

    Output:
    - out: scalar, negative log-likelihood
"""
@numba.njit
def loglike_rust(EExt, EEit, beta, theta, trans_mat):
    
    ## Get value function fixed point
    ev0 = valuefunction(beta, theta, trans_mat)
    x = np.arange(1,91).reshape(-1,1)
    c = 0.001 * theta[1] * x
    v1 = -theta[0] - c[0] + beta * ev0[0, 0]
    v0 = -c + beta * ev0
    payoff_diff = v1 - v0  # 90x1

    ## Transform states
    if np.nanmax(EExt) > 100:
        state = np.ceil(EExt / 5000).astype(int)
    else:
        state = EExt.astype(int)

    ## Compute log-likelihood
    f = np.full(EExt.shape, np.nan)
    T, N = EExt.shape
    for i in range(N):
        st = state[:, i]
        it = EEit[:, i]
        for t in range(T):
            if not np.isnan(st[t]):
                s_ind = int(st[t]) - 1  # zero-based indexing
                p1t = 1.0 - 1.0/(1.0 + np.exp(payoff_diff[s_ind]))
                if it[t] == 1:
                    f[t, i] = np.log(p1t)
                else:
                    f[t, i] = np.log(1 - p1t)

    out = -np.nansum(f)
    return out

"""
    Optimization
"""
## Objective function for scipy minimize
def objective(theta):
    return loglike_rust(EExt, EEit, beta, theta, trans_mat)

## Optimization
options = {
    'disp': True,
    'maxiter': 10**7,   # large maxiter
    'gtol': 1e-20       # similar to TolFun/TolX in MATLAB
}

res = minimize(objective, thetastart, method='BFGS', options=options)

theta_nfxp = res.x
hessian_pooled = res.hess_inv  # approximate Hessian inverse
stdtheta_nfxp = np.sqrt(np.diag(hessian_pooled))

## Save the results to local drive
savemat('results_nfxp.mat', {
    'phi': phi, 
    'trans_mat': trans_mat, 
    'theta_nfxp': theta_nfxp,
    'stdtheta_nfxp': stdtheta_nfxp
}, appendmat=True)