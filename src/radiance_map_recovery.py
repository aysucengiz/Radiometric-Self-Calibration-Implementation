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


    def update_coeff(rm):
        for i in range(rm.degree):
            rm.coeff[i] = rm.get_coeff(i)
        return rm.coeff
    
    def derivation_coeff(rm,n,p,q):
        # (a - b)**2
        # (a - b) * (a - b)'
        M_pq = rm.M[p][q]
        M_pq1 = rm.M[p][q+1]
        R_qq1 = rm.R_est_list[q]
        
        # compute a and b
        a = 0
        b = 0
        for n in range(rm.degree):
            a += rm.coeff[n] * M_pq**n
            b += rm.coeff[n] * M_pq**n
        
        b *= R_qq1
        
        temp1 = a - b
        temp2 = M_pq**n -  R_qq1 * M_pq1**n
        
        # MATRIX
            
            
            
        
    def get_coeff(rm,n):
        
        # f(1) = I_max
        if n - 1 is rm.degree - 1:
            return rm.I_max - sum(rm.coeff[:-1])
        
        else:
            sum = 0
            for q in range(rm.Q):
                for p in range(rm.P):
                    sum += rm.derivation_coeff(n,q,p)
            return sum


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