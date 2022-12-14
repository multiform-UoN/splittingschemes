import numpy as np
import scipy.sparse as sparse
import scipy.linalg as linalg
from scipy import stats

def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = int(np.round(2 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/2$"
    elif N == 2:
        return r"$\pi$"
    elif N % 2 > 0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)

##################################################################################################

def fvm_laplacian_1D(nu, leftBC, rightBC, N, L, kwargs):
    dx = np.divide(L, N)
    fBC = np.zeros(N)
    diag = np.zeros(N)
    lower = np.zeros(N-1)
    upper = np.zeros(N-1)
    
    for k in range(N): #k=0,1,...,N-1
        if k==0: #first cell
            if leftBC['type']=='dirichlet':
                diag[0] = -(3.5*nu(0.0, kwargs) + nu(dx, kwargs))
                upper[0] =  0.5*nu(0.0, kwargs) + nu(dx, kwargs)
                fBC[k] = -3.0*nu(0.0, kwargs)*leftBC.get('value')/np.square(dx) #cambio di segno perchè è termine noto
            if leftBC['type']=='neumann':
                diag[0] = -nu(dx, kwargs)
                upper[0] = nu(dx, kwargs)
                fBC[k] = nu(0.0, kwargs)*leftBC.get('value')/dx #cambio di segno perchè è termine noto
        elif k==N-1: #last cell
            if rightBC['type']=='dirichlet':
                diag[N-1] = -(3.5*nu(L, kwargs) + nu(L-dx, kwargs))
                lower[N-2] =  0.5*nu(L, kwargs) + nu(L-dx, kwargs)
                fBC[k] = -3.0*nu(L, kwargs)*rightBC.get('value')/np.square(dx) #cambio di segno perchè è termine noto
            elif rightBC['type']=='neumann':
                diag[N-1] = -nu(L-dx, kwargs)
                lower[N-2] = nu(L-dx, kwargs)
                fBC[k] = -nu(L, kwargs)*rightBC.get('value')/dx #cambio di segno perchè è termine noto
        else: #k=1,2,3,...,N-2
            diag[k] = -(nu(k*dx, kwargs) + nu((k+1)*dx, kwargs))
            lower[k-1] = nu(k*dx, kwargs)
            upper[k] = nu((k+1)*dx, kwargs)

    return sparse.csc_matrix(sparse.diags([diag/np.square(dx), lower/np.square(dx), upper/np.square(dx)], [0, -1, 1])), fBC

def fvm_reconstruct_1D(sol):
    """second order accuracy reconstructor for cell faces"""
    n = len(sol)
    rec = np.zeros(n+1)
    for i,val in enumerate(rec):
        if i == 0:
            rec[0] = 1.5*sol[0] - 0.5*sol[1]
        elif i==n:
            rec[n] = 1.5*sol[-1] - 0.5*sol[-2]
        else:
            rec[i] = 0.5*(sol[i-1]+sol[i])
    return rec

def slope(residual):
    """compute the slope of the error function in log scale"""
    x = np.arange(len(residual))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,np.log10(residual))
    return slope

def algebraic_residual(A, B, C, D, u, v, f1, f2):
    return np.sqrt(np.square(np.linalg.norm(f1-A@u-B@v)) + np.square(np.linalg.norm(f2-C@u-D@v)))

##################################################################################################

def method_BlockJacobi(A, B, C, D, f1, f2, nit, toll):
    u  = np.zeros(A.shape[0])
    v  = np.zeros(A.shape[0])
    u0 = np.zeros(A.shape[0])
    v0 = np.zeros(A.shape[0])
    res = []
    nitfinal = 0
    for i in range(nit):
        nitfinal += 1
        u = sparse.linalg.spsolve(A, f1 - B@v0)
        v = sparse.linalg.spsolve(D, f2 - C@u0)
        res.append(algebraic_residual(A, B, C, D, u, v, f1, f2))
        u0 = u
        v0 = v
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal

def method_BlockGaussSeidel(A, B, C, D, f1, f2, nit, toll):
    nitfinal = 0
    u = np.zeros(A.shape[0])
    v = np.zeros(A.shape[0])
    res = []
    for i in range(nit):
        nitfinal += 1
        u = sparse.linalg.spsolve(A, f1 - B@v)
        v = sparse.linalg.spsolve(D, f2 - C@u)
        res.append(algebraic_residual(A, B, C, D, u, v, f1, f2))
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal

# def method_BlockSOR(A, B, C, D, f1, f2, nit, alpha, toll):
#     nitfinal = 0
#     u = np.zeros(A.shape[0])
#     v = np.zeros(A.shape[0])
#     res = []
#     for i in range(nit):
#         nitfinal += 1
#         delta_u = sparse.linalg.spsolve(A, f1 - np.dot(A.toarray(), u) - np.dot(B.toarray(), v))
#         u = u + alpha*delta_u
#         delta_v = sparse.linalg.spsolve(D, f2 - np.dot(C.toarray(), u) - np.dot(D.toarray(), v))
#         v = v + delta_v
#         #v = v + alpha*delta_v # this affects the choice of alpha
#         res.append(algebraic_residual(A, B, C, D, u, v, f1, f2))
#         if res[-1] < toll: break
#     return u, v, np.array(res), nitfinal

def method_Lscheme(A, B, C, D, f1, f2, nit, toll, type, L):
    nitfinal = 0
    u = np.zeros(A.shape[0])
    v = np.zeros(A.shape[0])
    res = []
    for i in range(nit):
        nitfinal += 1

        # L on both
        if type=='both':
            u = sparse.linalg.spsolve((L*sparse.eye(A.shape[0]))+A, f1 + L*u - B@v)
            v = sparse.linalg.spsolve((L*sparse.eye(A.shape[0]))+D, f2 + L*v - C@u)
        
        # L on u
        if type=='L_on_u':
            u = sparse.linalg.spsolve((L*sparse.eye(A.shape[0]))+A, f1 - B@v + L*u)
            v = sparse.linalg.spsolve(D, f2 - C@u)

        # L on v
        if type=='L_on_v':
            v = sparse.linalg.spsolve((L*sparse.eye(A.shape[0]))+D, f2 - C@u + L*v)
            u = sparse.linalg.spsolve(A, f1 - B@v)

        res.append(algebraic_residual(A, B, C, D, u, v, f1, f2))
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal

def method_ShurPartialJacobi(A, B, C, D, f1, f2, nit, toll, type):
    nitfinal = 0
    AA    = sparse.csc_matrix(np.diag(A.diagonal(), 0))
    invAA = sparse.csc_matrix(np.diag(1/A.diagonal(), 0))
    DD    = sparse.csc_matrix(np.diag(D.diagonal(), 0))
    invDD = sparse.csc_matrix(np.diag(1/D.diagonal(), 0))
    u = np.zeros(A.shape[0])
    v = np.zeros(A.shape[0])
    res = []
    for i in range(nit):
        nitfinal += 1

        # ALTERNATE
        if type=='alternate':
            u = sparse.linalg.spsolve(A - B@invDD@C, f1 - B@invDD@(f2-(D-DD)@v))
            v = sparse.linalg.spsolve(D - C@invAA@B, f2 - C@invAA@(f1-(A-AA)@u))
        
        # APPROXIMATED SHUR ON u
        if type=='shur_on_u':
            u = sparse.linalg.spsolve(A - B@invDD@C, f1 - B@invDD@(f2-(D-DD)@v))
            v = sparse.linalg.spsolve(D, f2 - C@u)

        # APPROXIMATED SHUR ON v
        if type=='shur_on_v':
            v = sparse.linalg.spsolve(D - C@invAA@B, f2 - C@invAA@(f1-(A-AA)@u))
            u = sparse.linalg.spsolve(A, f1 - B@v)

        res.append(algebraic_residual(A, B, C, D, u, v, f1, f2))
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal

def method_SPJ_stab(A, B, C, D, f1, f2, nit, toll, type):
    nitfinal = 0
    AA    = np.diag(A.diagonal(), 0)
    invAA = np.diag(1/A.diagonal(), 0)
    DD    = np.diag(D.diagonal(), 0)
    invDD = np.diag(1/D.diagonal(), 0)
    u = np.zeros(A.shape[0])
    v = np.zeros(A.shape[0])
    res = []
    for i in range(nit):
        nitfinal += 1

        # ALTERNATE
        if type=='alternate':
            u = sparse.linalg.spsolve(A - B@invDD@C, f1 - B@v - B@invDD@C@u)
            v = sparse.linalg.spsolve(D - C@invAA@B, f2 - C@u - C@invAA@B@v)
        
        # APPROXIMATED SHUR ON u
        if type=='shur_on_u':
            u = sparse.linalg.spsolve(A - B@invDD@C, f1 - B@v - B@invDD@C@u)
            v = sparse.linalg.spsolve(D, f2 - C@u)

        # APPROXIMATED SHUR ON v
        if type=='shur_on_v':
            v = sparse.linalg.spsolve(D - C@invAA@B, f2 - C@u - C@invAA@B@v)
            u = sparse.linalg.spsolve(A, f1 - B@v)

        res.append(algebraic_residual(A, B, C, D, u, v, f1, f2))
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal

def method_ShurDualPartialJacobi(A, B, C, D, f1, f2, nit, toll):
    nitfinal = 0
    AA    = sparse.csc_matrix(np.diag(A.diagonal(), 0))
    invAA = sparse.csc_matrix(np.diag(1/A.diagonal(), 0))
    BB    = sparse.csc_matrix(np.diag(B.diagonal(), 0))
    CC    = sparse.csc_matrix(np.diag(C.diagonal(), 0))
    DD    = sparse.csc_matrix(np.diag(D.diagonal(), 0))
    invDD = sparse.csc_matrix(np.diag(1/D.diagonal(), 0))
    u = np.zeros(A.shape[0])
    v = np.zeros(A.shape[0])
    res = []
    for i in range(nit):
        nitfinal += 1

        # S2PJ alternate
        u = sparse.linalg.spsolve(A - BB@invDD@C, f1 - BB@invDD@f2 - (B-BB@invDD@D)@v)
        v = sparse.linalg.spsolve(D - CC@invAA@B, f2 - CC@invAA@f1 - (C-CC@invAA@A)@u)

        # S2PJ ON u
        # u = sparse.linalg.spsolve(A - BB@invDD@C, f1 - BB@invDD@f2 - (B-BB@invDD@D)@v)
        # v = sparse.linalg.spsolve(D, f2 - C@u)

        # S2PJ ON v
        # v = sparse.linalg.spsolve(D - CC@invAA@B, f2 - CC@invAA@f1 - (C-CC@invAA@A)@u)
        # u = sparse.linalg.spsolve(A, f1 - B@v)

        res.append(algebraic_residual(A, B, C, D, u, v, f1, f2))
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal
