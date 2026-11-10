import numpy as np

EPSILON = 1e-6

class Radiance_Map:
    def __init__(rm,Image, R_est_list, limit = 10):
        rm.limit = limit
        rm.R_est_list = R_est_list # check len of R_est_list
        rm.degree = 0
        rm.coeff = []
        rm.I_max = 0 # recover I_max
        rm.Q = 0
        rm.P = 0
        rm.M = Image
        rm.res_prev = None
        rm.res_curr = None


    
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
        
        # E = sum_m c_m * A_pqm + I_max (M_pq^N - R_qq1 * M_pq1**N)
        # deriv(E_sum_of_sq,c_n) = sum_p sum_q (E * A_pqn)
        #                        = sum_p sum_q  sum_m c_m * A_qpm * A_qpn
        #                        = sum_m c_m * (sum_q sum_p A_qpm * A_qpn)
        # Im dying here
        # B: sum_p sum_q A_pqm * A_pqn
        # record the coefficients of all cn
        L_n = []  
        for m in range(rm.degree):
            tot = 0
            for q in range(rm.Q-1):
                for p in range(rm.P):
                    tot += (A[q][p][m]-A[q][p][rm.degree])  * A[q][p][n]
            L_n.append(tot)
        
        R_n = 0  
        for q in range(rm.Q-1):
            for p in range(rm.P):
                R_n += - rm.I_max * A[q][p][n] * A[q][p][rm.degree]
                
        return L_n, R_n
            
        
    def update_coeff(rm):
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
        for q in range (rm.Q - 1):
            for p in range(rm.P):
                num = 0
                den = 0
                for n in range(rm.degree):
                    num += rm.coeff[n] * rm.M[p][q]**n
                    den += rm.coeff[n] * rm.M[p][q+1]**n
           
    def is_converged(rm) -> bool:
        # prev has rm.degree - 1
        if rm.res_prev is None:
            return False
        
        rm.res_curr = 0
        
        for p in range(rm.P):
            for q in range(rm.Q):
                for n in range(rm.degree):
                    rm.res_curr += rm.coeff * rm.M[p,q] ** n
                
        if abs(rm.res_prev - rm.res_curr) < EPSILON:
            return True
        else:
            return False
        
        
    def new_iteration(rm):
        rm.res_prev = rm.res_prev
        rm.degree += 1
         

    # order of images in increasing exposure
    # M: brightness measurement
    # p: pixel index
    # q: exposure index
    def radiance_map_recovery(rm):
        while rm.degree < rm.limit:
            rm.coeff.append(0)
            rm.coeff = rm.update_coeff()
            rm.R_est_list = rm.update_R_est()
            
            if rm.is_converged():
                break
                
            rm.new_iteration()
        return rm.degree, rm.coeff