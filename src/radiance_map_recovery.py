
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


    def update_coeff(rm):
        for i in range(rm.degree):
            rm.coeff[i] = rm.get_coeff(i)
        return rm.coeff
    
    def get_coeff(rm,n):
        
        
        # f(1) = I_max
        if n - 1 is rm.degree - 1:
            return rm.I_max - sum(rm.coeff[:-1])
        else:
            
            for q in range(rm.Q):
                for p in range(rm.P):
                    (rm.M[p,q] **n-1) * n
                    # TODO
            return 


    # order of images in increasing exposure
    # M: brightness measurement
    # p: pixel index
    # q: exposure index
    def radiance_map_recovery(rm):
        converged = False
        while rm.degree < rm.limit or not converged:
            rm.coeff.append(0)
            rm.coeff = rm.update_coeff()
            rm.R_est_list = rm.update_R_est()

            rm.degree += 1
        return rm.degree, rm.coeff