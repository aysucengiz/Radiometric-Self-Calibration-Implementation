import numpy as np

EPSILON_R = 1e-6
EPSILON_F = 1e-6
LIMIT_F = 10
LIMIT_R = 10

class Radiance_Map:
    def __init__(rm,Image, R_est_list):
        rm.R_est_list = R_est_list # check len of R_est_list
        rm.R_est_list_prev = None
        rm.degree = 1
        rm.coeff = []
        rm.coeff_prev = None
        rm.I_max = 0 # recover I_max
        rm.Q = 0
        rm.P = 0
        rm.M = Image


    
    def derivation_coeff(rm,n):
        
        # f(1) = I_max
        # eliminate c_N

        # A: precompute M_pq^n - R_q,q+1 * M_pq^n
        A = []
        for q in range(rm.Q-1):
            A.append([])
            R_qq1 = rm.R_est_list[q]
            for p in range(rm.P):
                A[q].append([])
                M_pq = rm.M[p][q]
                M_pq1 = rm.M[p][q+1]
                for m in range(rm.degree + 1):
                    temp = M_pq**m - R_qq1 * M_pq1**m
                    A[q][p].append(temp)
        
        # Im dying here
        # record the coefficients of all cn
        L_n = []  
        for m in range(rm.degree):
            tot = 0
            for q in range(rm.Q-1):
                for p in range(rm.P):
                    tot += (A[q][p][m]-A[q][p][rm.degree])  * (A[q][p][n]-A[q][p][rm.degree])
            L_n.append(tot)
        
        R_n = 0  
        for q in range(rm.Q-1):
            for p in range(rm.P):
                R_n -= rm.I_max * (A[q][p][n] - A[q][p][rm.degree]) * A[q][p][rm.degree]
                
        return L_n, R_n
            
        
    def update_coeff(rm):
        rm.coeff_prev = rm.coeff.copy() 
        # L C = R
        L = []
        R = []
        for n in range(rm.degree):
            L_n, R_n = rm.derivation_coeff(n)
            L.append(L_n) 
            R.append(R_n)
        
        C = np.linalg.solve(L, R)
        c_N = rm.I_max - np.sum(C)
        rm.coeff = list(C) + [c_N]



    def update_R_est(rm):
        rm.R_est_list_prev = rm.R_est_list.copy() 
        for q in range (rm.Q - 1):
            res = 0
            for p in range(rm.P):
                num = 0
                den = 0
                for n in range(rm.degree+1):
                    num += rm.coeff[n] * rm.M[p][q]**n
                    den += rm.coeff[n] * rm.M[p][q+1]**n
                res += num/den
            rm.R_est_list[q] = res
           
    def Func_converged(rm) -> bool:
        # prev has rm.degree - 1
        if rm.res_prev is None or rm.coeff_prev is None:
            return False
        
        diff = []
        
        for p in range(rm.P):
            for q in range(rm.Q-1):
                res_curr = 0
                res_prev = 0
                for n in range(rm.degree+1):
                    res_curr += rm.coeff[n] * rm.M[p][q] ** n
                for n in range(rm.degree):
                    res_prev += rm.coeff_prev[n] * rm.M[p][q] ** n
                diff.append(abs(res_curr - res_prev))
                
        if max(diff) < EPSILON_F:
            return True
        else:
            return False
        
        
    def new_iteration(rm):
        rm.res_prev = rm.res_curr
        rm.degree += 1

    def R_is_converged(rm):
        if rm.R_est_list_prev is None:
            return False
        
        errors = [abs(rm.R_est_list_prev[i]-rm.R_est_list[i]) for i in range(rm.Q-1)]
        
        rm.R_est_list_prev = rm.R_est_list 
        return max(errors) < EPSILON_R

    # order of images in increasing exposure
    # M: brightness measurement
    # p: pixel index
    # q: exposure index
    def radiance_map_recovery(rm):
        while rm.degree < rm.limit:
            rm.coeff = [0 for i in range(rm.degree+1)]
            while True:
                rm.update_coeff()
                rm.update_R_est()
                if rm.R_is_converged():
                    break
                
            if rm.Func_converged():
                break
                    
            rm.new_iteration()
        return rm.degree, rm.coeff