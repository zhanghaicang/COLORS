import sys
import numpy

def calc_lrs(lambda_, M):
    nrow, ncol = M.shape
    dim = np.minimum(nrow, ncol)
    Y = M
    
    L= np.zeros((nrow, ncol)); 
    S= np.zeros((nrow, ncol))
    Z= np.zeros((nrow, ncol))
    
    singular_values = np.linalg.svd(M, full_matrices=True, compute_uv=False)

    norm_two = np.max(singular_values)
    norm_inf = np.max(np.abs(Y)) / lambda_
    dual_norm = np.maximum(norm_two, norm_inf)
    
    d_norm = np.linalg.norm(M, ord=2)
    
    Y /= dual_norm
    
    mu = 1.25 / norm_two
    rho = 1.5
    mu_bar = mu * 1.0e+7
    converged = False
    max_iter = 10
    error_tolerance = 1.0e-7
    total_svd = 0
    sv = 10
    iter = 0
        
    for _ in range(max_iter):     
        
        #update sparse matrix S
        temp_T = M - L + (1.0 / mu) * Y;
        S = np.maximum(temp_T - lambda_ / mu, Z) + np.minimum(temp_T + lambda_ / mu, Z)
        S = np.maximum(S, Z)

        
        U, singular_values, V = np.linalg.svd(M - S + 1.0 / mu * Y, full_matrices=True,compute_uv=True)
        iter = iter+1

        svp = np.int(np.sum(np.greater(singular_values, 1.0 / mu)))
        print(svp, 'test')
        if svp < sv:
            sv = np.minimum(svp + 1, dim);
        else:
            sv = np.minimum(svp + np.int(0.05 * dim + 0.5), dim);

        S_th = np.diag(singular_values[:svp] - 1.0 / mu)

        L = np.matmul(np.matmul(U[:,:svp], S_th), V[:,:svp].transpose())
        L = np.maximum(L, Z)

        total_svd += 1;
        D = M - L - S;
        Y = Y + mu * D;

        mu = np.minimum(mu * rho, mu_bar);
        objective = np.linalg.norm(D, ord=2) / d_norm;
        print(objective, error_tolerance)
        if objective < error_tolerance:
            break;
      
    return S 
