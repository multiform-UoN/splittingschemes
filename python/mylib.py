import numpy as np
import scipy.sparse as sparse

def laplacian(nu, leftBC, rightBC, N, dx, L, kwargs):
    fBC = np.zeros(N)
    diag = np.zeros(N)
    lower = np.zeros(N-1)
    upper = np.zeros(N-1)
    
    for k in range(N): #k=0,1,...,N-1
        if k==0:
            if leftBC.get('type')=='dirichlet':
                diag[0] = -(3.5*nu(0.0, kwargs) + nu(dx, kwargs))
                upper[0] =  0.5*nu(0.0, kwargs) + nu(dx, kwargs)
                fBC[k] = -3.0*nu(0.0, kwargs)*leftBC.get('value')/np.square(dx) #cambio di segno perchè a dx dell'=
            if leftBC.get('type')=='neumann':
                diag[0] = -nu(dx, kwargs)
                upper[0] = nu(dx, kwargs)
                fBC[k] = nu(0.0, kwargs)*leftBC.get('value')/dx #cambio di segno perchè a dx dell'=
        elif k==N-1:
            if rightBC.get('type')=='dirichlet':
                diag[N-1] = -(3.5*nu(L, kwargs) + nu(L-dx, kwargs))
                lower[N-2] =  0.5*nu(L, kwargs) + nu(L-dx, kwargs)
                fBC[k] = -3.0*nu(L, kwargs)*rightBC.get('value')/np.square(dx) #cambio di segno perchè a dx dell'=
            elif rightBC.get('type')=='neumann':
                diag[N-1] = -nu(L-dx, kwargs)
                lower[N-2] = nu(L-dx, kwargs)
                fBC[k] = -nu(L, kwargs)*rightBC.get('value')/dx #cambio di segno perchè a dx dell'=
        else: #k=1,2,3,...,N-2
            diag[k] = -(nu(k*dx, kwargs) + nu((k+1)*dx, kwargs))
            lower[k-1] = nu(k*dx, kwargs)
            upper[k] = nu((k+1)*dx, kwargs)

    return sparse.diags([diag/np.square(dx), lower/np.square(dx), upper/np.square(dx)], [0, -1, 1]), fBC

# second order accuracy reconstructor for cell faces
def reconstruct(sol):
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
