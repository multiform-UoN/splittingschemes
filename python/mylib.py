import numpy as np
import scipy.sparse as sparse
import scipy.linalg as linalg
from scipy import stats

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

def method_BlockJacobi(A, B, C, D, f1, f2, nit, toll):
    u  = np.zeros(A.shape[0])
    v  = np.zeros(A.shape[0])
    u0 = np.zeros(A.shape[0])
    v0 = np.zeros(A.shape[0])
    res = []
    nitfinal = 0
    for i in range(nit):
        nitfinal += 1
        u = sparse.linalg.spsolve(A, f1 - np.dot(B.toarray(), v0))
        v = sparse.linalg.spsolve(D, f2 - np.dot(C.toarray(), u0))
        res.append(
            np.sqrt(
                np.square(np.linalg.norm(f1-np.dot(A.toarray(),u)-np.dot(B.toarray(),v)))
                +
                np.square(np.linalg.norm(f2-np.dot(C.toarray(),u)-np.dot(D.toarray(),v)))
                )
            )
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
        u = sparse.linalg.spsolve(A, f1 - np.dot(B.toarray(), v))
        v = sparse.linalg.spsolve(D, f2 - np.dot(C.toarray(), u))
        res.append(
            np.sqrt(
                np.square(np.linalg.norm(f1-np.dot(A.toarray(),u)-np.dot(B.toarray(),v)))
                +
                np.square(np.linalg.norm(f2-np.dot(C.toarray(),u)-np.dot(D.toarray(),v)))
                )
            )
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal

def method_BlockSOR(A, B, C, D, f1, f2, nit, alpha, toll):
    nitfinal = 0
    u = np.zeros(A.shape[0])
    v = np.zeros(A.shape[0])
    res = []
    for i in range(nit):
        nitfinal += 1
        delta_u = sparse.linalg.spsolve(A, f1 - np.dot(A.toarray(), u) - np.dot(B.toarray(), v))
        u = u + alpha*delta_u
        delta_v = sparse.linalg.spsolve(D, f2 - np.dot(C.toarray(), u) - np.dot(D.toarray(), v))
        v = v + delta_v
        #v = v + alpha*delta_v # this affects the choice of alpha
        res.append(
            np.sqrt(
                np.square(np.linalg.norm(f1-np.dot(A.toarray(),u)-np.dot(B.toarray(),v)))
                +
                np.square(np.linalg.norm(f2-np.dot(C.toarray(),u)-np.dot(D.toarray(),v)))
                )
            )
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal

def method_ShurPartialJacobi(A, B, C, D, f1, f2, nit, toll):
    nitfinal = 0
    DD    = np.diag(D.diagonal(),0)
    invDD = np.diag(1/D.diagonal(),0)
    AA    = np.diag(A.diagonal(),0)
    invAA = np.diag(1/A.diagonal(),0)
    u = np.zeros(A.shape[0])
    v = np.zeros(A.shape[0])
    res = []
    for i in range(nit):
        nitfinal += 1
        u = sparse.linalg.spsolve(
            A - sparse.csc_matrix(np.dot(B.toarray(), np.dot(invDD ,C.toarray()))),
            f1 - np.dot(B.toarray(), np.dot(invDD, f2 - np.dot(D.toarray()-DD, v)))
            )
        # v = sparse.linalg.spsolve(
        #     D - sparse.csc_matrix(np.dot(C.toarray(), np.dot(invAA ,C.toarray()))),
        #     f2 - np.dot(B.toarray(), np.dot(invAA, f1 - np.dot(A.toarray()-AA, u)))
        #     )
        v = sparse.linalg.spsolve(D, f2 - np.dot(C.toarray(),u))
        # u = sparse.linalg.spsolve(A, f1 - np.dot(B.toarray(),v))
        res.append(
            np.sqrt(
                np.square(np.linalg.norm(f1-np.dot(A.toarray(),u)-np.dot(B.toarray(),v)))
                +
                np.square(np.linalg.norm(f2-np.dot(C.toarray(),u)-np.dot(D.toarray(),v)))
                )
            )
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal

def method_ShurDualPartialJacobi(A, B, C, D, f1, f2, nit, toll):
    nitfinal = 0
    DD    = np.diag(D.diagonal(),0)
    invDD = np.diag(1/D.diagonal(),0)
    AA    = np.diag(D.diagonal(),0)
    invAA = np.diag(1/D.diagonal(),0)
    BB    = np.diag(B.diagonal(),0)
    u = np.zeros(A.shape[0])
    v = np.zeros(A.shape[0])
    res = []
    for i in range(nit):
        nitfinal += 1
        u = sparse.linalg.spsolve(
            A - sparse.csc_matrix(np.dot(BB, np.dot(invDD ,C.toarray()))),
            f1 - np.dot(B.toarray(), np.dot(invDD, f2 - np.dot(D.toarray()-DD, v))) + np.dot(B.toarray()-BB,np.dot(invDD,np.dot(C.toarray(),u)))
            )
        v = sparse.linalg.spsolve(D, f2 - np.dot(C.toarray(),u))
        # res.append(np.linalg.norm(sol - np.concatenate((u,v))))
        res.append(
            np.sqrt(
                np.square(np.linalg.norm(f1-np.dot(A.toarray(),u)-np.dot(B.toarray(),v)))
                +
                np.square(np.linalg.norm(f2-np.dot(C.toarray(),u)-np.dot(D.toarray(),v)))
                )
            )
        if res[-1] < toll: break
    return u, v, np.array(res), nitfinal













# def method_ShurPrecond(A, B, C, D, f1, f2, nit, L, sol):
#     #sA = sparse.csc_matrix(D)
#     #sA_iLU = sparse.linalg.spilu(sA)
#     #M = sparse.linalg.LinearOperator((A.shape[0],A.shape[0]), sA_iLU.solve)
#     #invDD = M.toarray()
#     DD    = np.diag(D.diagonal(),0)
#     invDD = np.diag(1/D.diagonal(),0)
#     u = np.zeros(A.shape[0])
#     v = np.zeros(A.shape[0])
#     res = []
#     for i in range(nit):
#         u = sparse.linalg.spsolve(
#             A  - np.dot(B, np.dot(invDD ,C)) , 
#             f1 - np.dot(B, np.dot(invDD, f2 - np.dot(D-DD, v)))
#         )
#         v = sparse.linalg.spsolve(D, f2 - np.dot(C,u))
#         res.append(np.linalg.norm(sol - np.concatenate((u,v))))
#     return u, v, np.array(res)


# def approxInverse(A, n):
#     Id = np.eye(A.shape[0])
#     invA = Id
#     prev = Id-A
#     for i in range(n):
#         invA = invA + prev
#         prev = np.dot(prev, Id-A)
#     return invA

# def approxInverseEpsilon(A, n, epsilon):
#     Id = np.eye(A.shape[0])
#     invA = Id
#     prev = Id-epsilon*A
#     for i in range(n):
#         invA = invA + prev
#         prev = np.dot(prev, Id-epsilon*A)
#     return epsilon*invA

# def method_ShurApproxinv(A, B, C, D, f1, f2, nit, sol, N, epsilon):
#     invA = approxInverseEpsilon(A, N, epsilon)
#     #invA = 2.0*np.diag(1/A.diagonal(),0) - np.dot(np.diag(1/A.diagonal(),0),np.dot(A,np.diag(1/A.diagonal(),0)))
#     u = np.zeros(A.shape[0])
#     v = np.zeros(A.shape[0])
#     res = []
#     for i in range(nit):
#         v = sparse.linalg.spsolve(
#             A - np.dot(C, np.dot(invA ,D)), 
#             f2 - np.dot(B, np.dot(C, np.dot(invA, f1)))
#         )
#         v = sparse.linalg.spsolve(A, f1 - np.dot(B, v))
#         res.append(np.linalg.norm(sol - np.concatenate((u,v))))
#     return u, v, np.array(res)

# def solveApproxInv(A, b, x0, nit, N, epsilon):
#     x = x0
#     AA = epsilon*A
#     operator = np.eye(A.shape[0])
#     prev = np.eye(A.shape[0])
#     for i in range(N):
#         operator = operator + prev*AA
#         prev = prev*AA
    
#     for i in range(nit):
#         x = np.dot(operator, x - epsilon*b)

#     return x
