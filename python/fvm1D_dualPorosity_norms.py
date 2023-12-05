from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from scipy.sparse import diags
import scipy.sparse as sparse
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib import markers
import mylib


def u_sol(x, kwargs):
    return kwargs['N'] + kwargs['O']*np.sin(kwargs['P']*(x-kwargs['Q']))

def dudx_sol(x, kwargs):
    return kwargs['O']*kwargs['P']*np.cos(kwargs['P']*(x-kwargs['Q']))

def d2udx2_sol(x, kwargs):
    return -kwargs['O']*np.square(kwargs['P'])*np.sin(kwargs['P']*(x-kwargs['Q']))


def v_sol(x, kwargs):
    return kwargs['R'] + kwargs['S']*np.exp(-kwargs['T']*x)

def dvdx_sol(x, kwargs):
    return -kwargs['S']*kwargs['T']*np.exp(-kwargs['T']*x)

def d2vdx2_sol(x, kwargs):
    return kwargs['S']*np.square(kwargs['T'])*np.exp(-kwargs['T']*x)


def m_u(x, kwargs):
    return kwargs['D'] + kwargs['E']*np.sin(kwargs['F']*(x-kwargs['G']))

def dm_udx(x, kwargs):
    return kwargs['E']*kwargs['F']*np.cos(kwargs['F']*(x-kwargs['G']))


def m_v(x, kwargs):
    return kwargs['H'] + kwargs['I']*np.sin(kwargs['L']*(x-kwargs['M']))

def dm_vdx(x, kwargs):
    return kwargs['I']*kwargs['L']*np.cos(kwargs['L']*(x-kwargs['M']))


def f_1(x, kwargs):
    return kwargs['beta']*(u_sol(x, kwargs)-v_sol(x, kwargs)) - (dm_udx(x, kwargs)*dudx_sol(x, kwargs) + m_u(x, kwargs)*d2udx2_sol(x, kwargs))

def f_2(x, kwargs):
    return kwargs['beta']*(v_sol(x, kwargs)-u_sol(x, kwargs)) - (dm_vdx(x, kwargs)*dvdx_sol(x, kwargs) + m_v(x, kwargs)*d2vdx2_sol(x, kwargs))



L = np.pi
N = 128
x = np.linspace(0, L, N+1) #faces coordinates
xc = np.linspace(0.5*np.divide(L, N), L-0.5*np.divide(L, N), N) #centercells coordinates
xfine = np.linspace(0, L, 1000)
omega = 1.5 # SOR relaxation factor
nit_max = 100
toll = 1e-8
mesh_size = np.power(2, np.arange(5, 12))

norm_A = []
norm_B = []
norm_C = []
norm_D = []
norm_Z = []


for val in mesh_size:
    kwargs = {
        'beta': 100,

        # m_u and m_v 
        'D': 1e4, #constant
        'E': 5e3, #amplitude
        'F': 2.0, #frequency
        'G': 0.0, #shift

        'H': 1,   #1000.0, #constant
        'I': 0.5, #-999.999, #amplitude
        'L': 4,   #frequency
        'M': 0.0, #shift

        # solutions
        'N': 0,
        'O': 1,
        'P': 2,
        'Q': 0.0,

        'R': 0,
        'S': 1.0,
        'T': 2.0
    }  

    leftBC_u  = {'type':'dirichlet', 'value':u_sol(0.0, kwargs)} #'type':'neumann',
    rightBC_u = {'type':'dirichlet', 'value':u_sol(L,   kwargs)}
    leftBC_v  = {'type':'dirichlet', 'value':v_sol(0.0, kwargs)}
    rightBC_v = {'type':'dirichlet', 'value':v_sol(L,   kwargs)}

    A_laplacian, f1BC = mylib.fvm_laplacian_1D(m_u, leftBC_u, rightBC_u, N, L, kwargs)
    D_laplacian, f2BC = mylib.fvm_laplacian_1D(m_v, leftBC_v, rightBC_v, N, L, kwargs)
    A = sparse.csc_matrix( kwargs['beta']*sparse.eye(N) - A_laplacian)
    B = sparse.csc_matrix(-kwargs['beta']*sparse.eye(N))
    C = sparse.csc_matrix(-kwargs['beta']*sparse.eye(N))
    D = sparse.csc_matrix( kwargs['beta']*sparse.eye(N) - D_laplacian)

    f1 = f_1(xc, kwargs) - f1BC
    f2 = f_2(xc, kwargs) - f2BC

    Z = np.zeros((2*N,2*N))
    Z[:N,:N] = A.toarray()
    Z[:N,N:] = B.toarray()
    Z[N:,:N] = C.toarray()
    Z[N:,N:] = D.toarray()
    Z = sparse.csc_matrix(Z)
    f = np.concatenate((f1,f2))
    
    norm_A.append(linalg.norm(A.toarray(), ord=2))
    norm_B.append(linalg.norm(B.toarray(), ord=2))
    norm_C.append(linalg.norm(C.toarray(), ord=2))
    norm_D.append(linalg.norm(D.toarray(), ord=2))
    norm_Z.append(linalg.norm(Z.toarray(), ord=2))

plt.plot(mesh_size, norm_A, '.-', label="A")
plt.plot(mesh_size, norm_B, '.-', label="B")
plt.plot(mesh_size, norm_C, '.-', label="C")
plt.plot(mesh_size, norm_D, '.-', label="D")
plt.plot(mesh_size, norm_Z, '.-', label="Z")
plt.legend()
plt.grid()
plt.show()
