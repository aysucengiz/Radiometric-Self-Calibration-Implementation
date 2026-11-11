import numpy as np
import matplotlib.pyplot as plt
import csv
import json
import math
EPSILON_R = 1e-4
LIMIT_F = 5
LIMIT_R = 50
EPSILON_M = 0.02


class Radiance_Map:
    def __init__(rm,Image, R_est_list):
        rm.R_est_list_orig = R_est_list.copy() # check len of R_est_list -> gerek kalmadı buna numpya geçince
        rm.R_est_list = R_est_list.copy() # check len of R_est_list -> bunlar hala kullanılıyor mu
        rm.R_est_list_prev = None
        rm.degree = 0
        rm.coeff = []
        rm.I_max = 1
        rm.M =  Image.astype(np.float64) / 255.0
        rm.M = rm.M[:, np.random.default_rng(0).choice(rm.M.shape[1], min(20000, rm.M.shape[1]), replace=False)] # random pixel subsample
        rm.Q = rm.M.shape[0]
        rm.P = rm.M.shape[1]
        rm.final_degree = 1
        rm.M_n = np.stack([rm.M ** n for n in range(LIMIT_F + 1)],axis=0)

    # numpya geçince birleştirince daha hızlı oldu
    # def derivation_coeff(rm, n):
    #
    #     R = np.asarray(rm.R_est_list)
    #
    #     A = rm.M_n[1:rm.degree + 1, :-1, :] - (
    #             R[None, :, None] * rm.M_n[1:rm.degree + 1, 1:, :]
    #     )
    #
    #     A_diff = A - A[rm.degree - 1, :, :][None, :, :]
    #
    #     valid_pairs = (
    #             (rm.M[:-1, :] > EPSILON_M) &
    #             (rm.M[:-1, :] < 1 - EPSILON_M) &
    #             (rm.M[1:, :] > EPSILON_M) &
    #             (rm.M[1:, :] < 1 - EPSILON_M)
    #     )
    #
    #     L_n = np.sum(
    #         A_diff[:rm.degree - 1] *
    #         A_diff[n][None, :, :] *
    #         valid_pairs[None, :, :],
    #         axis=(1, 2)
    #     )
    #
    #     R_n = -rm.I_max * np.sum(
    #         A_diff[n] * A[rm.degree - 1] * valid_pairs
    #     )
    #
    #     return L_n, R_n


    def update_coeff(rm):
        # L C = R
        R = np.asarray(rm.R_est_list)

        if rm.degree == 1:
            rm.coeff = [0.0, rm.I_max]
            return

        A1 = rm.M_n[1:rm.degree + 1, :-1, :]
        A2 =  R[None, :, None] * rm.M_n[1:rm.degree + 1, 1:, :]
        A = A1 - A2
        A_diff = A[:-1] - A[-1][None, :, :]

        valid_pairs = (
                (rm.M[:-1, :] > EPSILON_M) &
                (rm.M[:-1, :] < 1 - EPSILON_M) &
                (rm.M[1:, :] > EPSILON_M) &
                (rm.M[1:, :] < 1 - EPSILON_M)
        )

        X = A_diff[:, valid_pairs]
        Y = A[-1, valid_pairs]

        L = X @ X.T
        RHS = -rm.I_max * (X @ Y)
        try:
            C = np.linalg.solve(L, RHS)
        except np.linalg.LinAlgError:
            raise ValueError(
                f"Singular coefficient system at degree {rm.degree}"
            )

        c_N = rm.I_max - np.sum(C)
        rm.coeff = [0.0] + list(C) + [c_N]



    def update_R_est(rm):
        rm.R_est_list_prev = rm.R_est_list.copy()

        coeff = np.asarray(rm.coeff)

        # Evaluate the recovered response for every image and pixel
        f_M = np.sum(
            coeff[:, None, None] * rm.M_n[:rm.degree + 1],
            axis=0
        )

        valids = (
                (rm.M > EPSILON_M) &
                (rm.M < 1 - EPSILON_M)
        )

        for q in range(rm.Q - 1):
            valid = valids[q] & valids[q + 1]

            if not np.any(valid):
                continue

            numerator = np.sum(
                f_M[q, valid] * f_M[q + 1, valid]
            )

            denominator = np.sum(
                f_M[q + 1, valid] ** 2
            )

            if denominator < 1e-12:
                continue

            rm.R_est_list[q] = numerator / denominator

        rm.R_est_list[0] = rm.R_est_list_orig[0] # anchor R[0]: removes the f^t / R^t ambiguity

    def compute_error(rm):
        R = np.asarray(rm.R_est_list)

        A = rm.M_n[1:rm.degree + 1, :-1, :] - (
                R[None, :, None] * rm.M_n[1:rm.degree + 1, 1:, :]
        )

        valid_pairs = (
                (rm.M[:-1, :] > EPSILON_M) &
                (rm.M[:-1, :] < 1 - EPSILON_M) &
                (rm.M[1:, :] > EPSILON_M) &
                (rm.M[1:, :] < 1 - EPSILON_M)
        )

        coeff = np.asarray(rm.coeff[1:])

        residual = (coeff[:, None, None] * A).sum(axis=0)

        if not np.any(valid_pairs):
            return float('inf')

        return np.mean(residual[valid_pairs] ** 2)


    def new_iteration(rm):
        rm.R_est_list_prev = None

    def R_is_converged(rm):
        if rm.R_est_list_prev is None:
            return False

        errors = [abs(rm.R_est_list_prev[i]-rm.R_est_list[i]) for i in range(rm.Q-1)]
        if max(errors) < EPSILON_R:
            return True
        # print(f"R is not converged yet. Max error: {max(errors)}. R estimations: {rm.R_est_list}")
        return False

    def choose_best(rm):
        cand = {d: rm.errors[d] for d in range(1, LIMIT_F + 1) if rm.coeff_final[d]}
        if not cand: return -1
        best = min(cand.values())
        return min(d for d, e in cand.items() if e <= best * 1.05)


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
                rm.update_R_est()
                if rm.R_is_converged():
                    break

            rm.update_coeff()
            rm.errors[rm.degree] = rm.compute_error()
            rm.coeff_final[rm.degree] = rm.coeff.copy()
            rm.Rs[rm.degree] = rm.R_est_list.copy()
            print(f"Error: {rm.errors[rm.degree]}.\n Coeffs: {rm.coeff_final[rm.degree]}")
        rm.final_degree = rm.choose_best()

        for i in range(1, LIMIT_F + 1):
            print(f"Degree {i}: {rm.coeff_final[i]}")
            print(f"Degree {i}: {rm.errors[i]}")
        if rm.final_degree == -1:
            print("No degree converged")
            return False
        else:
            print(f"Degree {rm.final_degree} has been chosen.")


    def plot_response_function(rm, save_path):

        if rm.final_degree == -1 or not rm.coeff_final[rm.final_degree]:
            print("No response function :(")
            return

        coeff = rm.coeff_final[rm.final_degree]
        x = np.linspace(0, 1, 256)

        y = np.zeros_like(x)
        for n, c in enumerate(coeff):
            y += c * x**n

        plt.figure(figsize=(8, 5))
        plt.plot(x, y, label=f"Response Function (degree {rm.final_degree})")

        # take inverse of hw pdf
        y_srgb = np.where( x <= 0.0031308, x / 12.92,((x + 0.055) / 1.055) ** 2.4)


        plt.plot(x, y_srgb, label="sRGB", color="orange")
        plt.plot(x, x, "--", label="Linear ", color="gray")
        plt.xlabel("M")
        plt.ylabel("f(M)")
        plt.title("Camera Response Function")
        plt.xlim(0, 1)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")


        plt.close()

    def save_results(rm, filepath, time):

        with open(filepath, "w", newline="") as file:
            writer = csv.writer(file)

            writer.writerow([
                "Degree",
                "Error",
                "Time (seconds)",
                "Coefficients",
                "Exposure Ratios"
            ])

            for d in range(1, LIMIT_F + 1):
                if not rm.coeff_final[d]:
                    continue

                writer.writerow([
                    d,
                    rm.errors[d],
                    time[d] if isinstance(time, dict) else time,
                    json.dumps(rm.coeff_final[d]),
                    json.dumps(rm.Rs[d])
                ])

        print(f"Results saved to: {filepath}")