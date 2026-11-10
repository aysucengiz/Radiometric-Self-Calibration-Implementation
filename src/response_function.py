import numpy as np
import math
EPSILON_R = 1e-6
LIMIT_F = 10
LIMIT_R = 10

class Radiance_Map:
    def __init__(rm,Image, R_est_list):
        rm.R_est_list_orig = R_est_list.copy() # check len of R_est_list
        rm.R_est_list = R_est_list.copy() # check len of R_est_list
        rm.R_est_list_prev = None
        rm.degree = 0
        rm.coeff = []
        rm.I_max = 1 
        rm.M =  Image.astype(np.float64) / 255.0
        rm.Q = rm.M.shape[1]
        rm.P = rm.M.shape[0]
        rm.final_degree = 1


    
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
        # L C = R
        L = []
        R = []
        for n in range(rm.degree):
            L_n, R_n = rm.derivation_coeff(n)
            L.append(L_n) 
            R.append(R_n)
        
        try:
            C = np.linalg.solve(L, R)
        except np.linalg.LinAlgError:
            raise ValueError(f"Singular coefficient system at degree {rm.degree}")

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
                if abs(den) < 1e-12:
                    print("Denominator was zero!!")
                else:
                    res += num/den
            rm.R_est_list[q] = res
           
    def compute_error(rm):
        err = 0
        
        for q in range(rm.Q-1):
            for p in range(rm.P):
                tot = 0
                M_pq = rm.M[p][q]
                M_pq1 = rm.M[p][q+1]
                R_qq1 = rm.R_est_list[q]
                for n in range(rm.degree + 1):
                    tot += rm.coeff[n] * (M_pq**n - R_qq1 * M_pq1**n)
                err += tot**2
        return err
        
        
    def new_iteration(rm):
        rm.R_est_list_prev = None

    def R_is_converged(rm):
        if rm.R_est_list_prev is None:
            return False
        
        errors = [abs(rm.R_est_list_prev[i]-rm.R_est_list[i]) for i in range(rm.Q-1)]
        return max(errors) < EPSILON_R

    def choose_best(rm):
        min_ = float('inf')
        curr_degree = -1
        for i in range(1,LIMIT_F+1):
            if rm.errors[i] < min_:
                min_ = rm.errors[i]
                curr_degree = i
        return curr_degree
            

    # order of images in increasing exposure
    # M: brightness measurement
    # p: pixel index
    # q: exposure index
    def response_function(rm):
        rm.errors = [float('inf') for _ in range(LIMIT_F+1)]
        rm.coeff_final = [[] for i in range(LIMIT_F+1)]
        rm.Rs = [[] for i in range(LIMIT_F+1)]
        for deg in range(1,LIMIT_F+1):
            rm.degree = deg
            rm.new_iteration()
            rm.R_est_list = rm.R_est_list_orig.copy()
            #polynomial degree
            rm.coeff = [0 for i in range(rm.degree+1)]
            for _ in range(LIMIT_R):
                # r est coeff computation loop
                rm.update_coeff()
                rm.update_R_est()
                if rm.R_is_converged():
                    rm.errors[rm.degree] = rm.compute_error()
                    rm.coeff_final[rm.degree] = rm.coeff.copy()  
                    rm.Rs[rm.degree] = rm.R_est_list.copy() 
                    break
                   
        rm.final_degree = rm.choose_best()
        if rm.final_degree == -1:
            print("No degree converged")
            return False
        else:
            print(f"Degree {rm.final_degree} has been chosen.")
        
    
import matplotlib.pyplot as plt
def plot_response_function(rm):

    if rm.final_degree == -1 or not rm.coeff_final[rm.final_degree]:
        print("No response function :(")
        return

    coeff = rm.coeff_final[rm.final_degree]
    x = np.linspace(0, 1, 256)

    y = np.zeros_like(x)
    for n, c in enumerate(coeff):
        y += c * x**n

    plt.figure(figsize=(8, 5))
    plt.plot(x, y, label=f"Recovered response (degree {rm.final_degree})")

    plt.xlabel("Pixel intensity (M)")
    plt.ylabel("Radiance (f(M))")
    plt.title("Recovered Camera Response Function")
    plt.xlim(0, 1)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()