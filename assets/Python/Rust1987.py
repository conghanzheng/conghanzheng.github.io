"""
Preliminaries
"""
## Package Dependencies
# !pip install numpy scipy
import numpy as np
from scipy.optimize import minimize, approx_fprime
from scipy.stats import norm
import logging

logging.basicConfig(level=logging.INFO)

"""
Load online data
"""

url = 'https://raw.githubusercontent.com/conghanzheng/conghanzheng.github.io/master/assets/Python/cleandata.mat'
with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
    urllib.request.urlretrieve(url, tmp_file.name)

data = loadmat(tmp_file.name)

os.unlink(tmp_file.name)

EExt = data['EExt']  ## Cumulative mileage
EEit = data['EEit']  ## Dummy indicator matrix

"""
Transitions
"""

def trans_prob(EExt):
    """
    Function: Calculate Empirical Transition Probabilities
    """
    EExt_clean = np.nan_to_num(EExt, nan=0)
    interv = np.ceil(EExt_clean / 5000).astype(int)
    dfint = np.diff(interv, axis=0)
    
    phi = np.array([np.mean(dfint == i) for i in range(3)])
    phi_sum = np.sum(phi)
    if phi_sum == 0:
        logging.warning("Zero transition probabilities detected")
        return np.array([1/3, 1/3, 1/3])
    return phi / phi_sum

def tran_mat(phi):
    """
    Function: Build transition matrix from probabilities
    """
    F = np.zeros((90, 90))
    np.fill_diagonal(F[:89, :89], phi[0])
    np.fill_diagonal(F[:89, 1:], phi[1])
    np.fill_diagonal(F[:88, 2:], phi[2])
    F[89, 89] = 1
    F[88, 89] = 1 - phi[0]
    return F

"""
NFXP Steps
"""

def value_function(beta, theta, trans_mat):
    """
    Function: Compute Value Function
    """
    ev = np.zeros((90, 1))
    tol, diff = 1e-6, 1.0
    max_iter = 1000
    iter_count = 0
    
    while diff > tol and iter_count < max_iter:
        c = 0.001 * theta[1] * np.arange(1, 91).reshape(-1, 1)
        max_v = np.maximum(-theta[0] + beta * ev[0, 0], -c + beta * ev)
        log_sum = np.log(np.exp(-theta[0] + beta * ev[0, 0] - max_v) + 
                        np.exp(-c + beta * ev - max_v))
        ev_new = trans_mat @ (max_v + log_sum)
        diff = np.linalg.norm(ev_new - ev)
        ev = ev_new
        iter_count += 1
    
    if iter_count == max_iter:
        logging.warning("Value function iteration did not converge")
    
    return ev

def log_likelihood(theta, beta, trans_mat, state, choices):
    """
    Function: Compute Log-Likelihood
    """
    if np.any(theta <= 0):
        return -1e10
    
    try:
        ev = value_function(beta, theta, trans_mat)
        c = 0.001 * theta[1] * np.arange(1, 91).reshape(-1, 1)
        v_diff = -theta[0] - c[0] + beta * ev[0, 0] - (-c + beta * ev)
        
        loglike = 0.0
        valid_obs = 0
        
        for i in range(state.shape[1]):
            for t in range(state.shape[0]):
                if not np.isnan(state[t, i]):
                    s_ind = int(state[t, i]) - 1
                    if 0 <= s_ind < 90:  # Ensure valid state index
                        p1 = 1.0 / (1.0 + np.exp(-v_diff[s_ind]))
                        loglike += np.log(p1 if choices[t, i] == 1 else 1 - p1)
                        valid_obs += 1
        
        logging.info(f"Log-likelihood computed with {valid_obs} valid observations")
        return loglike
    except Exception as e:
        logging.error(f"Error in log-likelihood computation: {str(e)}")
        return -1e10

def estimate_rust(EExt, EEit, beta):
    """
    Function: Main Rust NFXP Estimation
    """
    logging.info("Starting Rust model estimation")
    
    # States
    state = np.ceil(EExt / 5000).astype(int) if np.nanmax(EExt) > 100 else EExt.astype(int)
    
    # Transitions
    phi = trans_prob(EExt)
    logging.info(f"Transition probabilities: {phi}")
    trans_mat = tran_mat(phi)
    
    # Optimization
    objective = lambda theta: -log_likelihood(theta, beta, trans_mat, state, EEit)
    
    try:
        result = minimize(objective, [10.0, 2.0], method='L-BFGS-B',
                         bounds=[(1e-10, None), (1e-10, None)])
        
        n_obs = np.sum(~np.isnan(state))
        
        # Modified Hessian calculation with more stable step size
        step_size = np.maximum(1e-4 * np.abs(result.x), 1e-4)  # Adaptive step size
        
        def compute_hessian(x, func, h):
            n = len(x)
            hessian = np.zeros((n, n))
            
            for i in range(n):
                for j in range(n):
                    x_ij = x.copy()
                    x_ij[i] += h[i]
                    x_ij[j] += h[j]
                    
                    x_i = x.copy()
                    x_i[i] += h[i]
                    
                    x_j = x.copy()
                    x_j[j] += h[j]
                    
                    hessian[i,j] = (func(x_ij) - func(x_i) - func(x_j) + func(x)) / (h[i] * h[j])
            
            return (hessian + hessian.T) / 2  # Ensure symmetry
        
        hessian = compute_hessian(result.x, objective, step_size)
        
        # Add small diagonal perturbation if necessary
        min_eig = np.min(np.linalg.eigvals(hessian))
        if min_eig < 1e-6:
            hessian += np.eye(len(result.x)) * (1e-6 - min_eig)
        
        vcov = np.linalg.inv(hessian) / n_obs
        se = np.sqrt(np.diag(vcov))
        t_stats = result.x / se
        p_vals = 2 * (1 - norm.cdf(np.abs(t_stats)))
        
        # Results
        param_names = ['RC cost (θ_r)', 'MC cost (θ_m)']
        print("\nRust Model Estimation Results:")
        print("=" * 65)
        print(f"{'Parameter':<12} {'Estimate':>10} {'Std.Err':>10} {'t-stat':>10} {'p-value':>10}")
        print("-" * 65)
        
        for i, name in enumerate(param_names):
            print(f"{name:<12} {result.x[i]:10.4f} {se[i]:10.4f} {t_stats[i]:10.4f} {p_vals[i]:10.4f}")
        
        print(f"\nLog-likelihood: {-result.fun:.4f}")
        print(f"Observations: {n_obs}")
        print(f"Converged: {result.success}")
        
        return result.x, se, t_stats, p_vals
        
    except Exception as e:
        logging.error(f"Estimation failed: {str(e)}")
        return None, None, None, None

"""
Optimization
"""
beta = 0.99
estimates, std_errors, t_stats, p_values = estimate_rust(EExt, EEit, beta)