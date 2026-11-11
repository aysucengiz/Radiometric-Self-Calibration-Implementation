import numpy as np
import matplotlib.pyplot as plt
import math
EPSILON_R = 1e-2
LIMIT_F = 7
LIMIT_R = 10
EPSILON_M = 0.02

class Radiance_Map:
    def __init__(rm,Image, R_est_list):
        rm.R_est_list_orig = R_est_list.copy() # check len of R_est_list
        rm.R_est_list = R_est_list.copy() # check len of R_est_list
        rm.R_est_list_prev = None
        rm.degree = 0
        rm.coeff = []
        rm.I_max = 1 
        rm.M =  Image.astype(np.float64) / 255.0
        rm.Q = rm.M.shape[0]
        rm.P = rm.M.shape[1]
        rm.final_degree = 1
        rm.M_n = np.stack([rm.M ** n for n in range(LIMIT_F + 1)],axis=0)



    def derivation_coeff(rm, n):

        R = np.asarray(rm.R_est_list)

        A = rm.M_n[:rm.degree + 1, :-1, :] - (
                R[None, :, None] * rm.M_n[:rm.degree + 1, 1:, :]
        )

        A_diff = A - A[rm.degree, :, :][None, :, :]

        valid_pairs = (
                (rm.M[:-1, :] > EPSILON_M) &
                (rm.M[:-1, :] < 1 - EPSILON_M) &
                (rm.M[1:, :] > EPSILON_M) &
                (rm.M[1:, :] < 1 - EPSILON_M)
        )

        L_n = np.sum(
            A_diff[:rm.degree] *
            A_diff[n][None, :, :] *
            valid_pairs[None, :, :],
            axis=(1, 2)
        )

        R_n = -rm.I_max * np.sum(
            A_diff[n] * A[rm.degree] * valid_pairs
        )

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
            print(f"L: {L}")
            print(f"R: {R}")
            C = np.linalg.solve(L, R)
        except np.linalg.LinAlgError:
            raise ValueError(f"Singular coefficient system at degree {rm.degree}")

        c_N = rm.I_max - np.sum(C)
        rm.coeff = list(C) + [c_N]



    def update_R_est(rm):
        rm.R_est_list_prev = rm.R_est_list.copy() 
        for q in range (rm.Q - 1):
            res = 0
            pixels_used = 0
            for p in range(rm.P):
                num = 0
                den = 0
                if rm.M_n[1][q][p] <= EPSILON_M or rm.M_n[1][q][p] >= 1- EPSILON_M:
                    continue
                if rm.M_n[1][q+1][p] <= EPSILON_M or rm.M_n[1][q+1][p] >= 1- EPSILON_M:
                    continue
                for n in range(rm.degree+1):
                    num += rm.coeff[n] * rm.M_n[n][q][p]
                    den += rm.coeff[n] * rm.M_n[n][q+1][p]
                if abs(den) < 1e-12:
                    print("Denominator was zero!!")
                else:
                    res += num/den
                    pixels_used += 1
            rm.R_est_list[q] = res / pixels_used
           
    def compute_error(rm):
        M_n = rm.M_n[:rm.degree + 1]

        R = np.asarray(rm.R_est_list)
        coeff = np.asarray(rm.coeff)

        A = M_n[:, :-1, :] - R[None, :, None] * M_n[:, 1:, :]

        err = np.sum((coeff[:, None, None] * A).sum(axis=0) ** 2)
        return err
        
        
    def new_iteration(rm):
        rm.R_est_list_prev = None

    def R_is_converged(rm):
        if rm.R_est_list_prev is None:
            return False
        
        errors = [abs(rm.R_est_list_prev[i]-rm.R_est_list[i]) for i in range(rm.Q-1)]
        print(f"R is not converged yet. Max error: {max(errors)}. R estimations: {rm.R_est_list}")
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
            print(f"Degree: {deg}")
            rm.degree = deg
            rm.new_iteration()
            rm.R_est_list = rm.R_est_list_orig.copy()
            #polynomial degree
            rm.coeff = [0 for i in range(rm.degree+1)]
            for _ in range(LIMIT_R):
                # r est coeff computation loop
                rm.update_coeff()
                rm.errors[rm.degree] = rm.compute_error()
                rm.coeff_final[rm.degree] = rm.coeff.copy()  
                rm.Rs[rm.degree] = rm.R_est_list.copy()  
                print(f"Error: {rm.errors[rm.degree]}.\n Coeffs: {rm.coeff_final[rm.degree]}")
                break
                """rm.update_R_est()
                if rm.R_is_converged():
                    rm.errors[rm.degree] = rm.compute_error()
                    rm.coeff_final[rm.degree] = rm.coeff.copy()  
                    rm.Rs[rm.degree] = rm.R_est_list.copy() 
                    break"""
                   
        rm.final_degree = rm.choose_best()
        if rm.final_degree == -1:
            print("No degree converged")
            return False
        else:
            print(f"Degree {rm.final_degree} has been chosen.")
        
    
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