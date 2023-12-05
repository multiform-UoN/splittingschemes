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
nit_max = 1000
toll = 1e-12

kwargs = {
    'beta': 1e4,

    # m_u and m_v 
    'D': 1e4, #constant
    'E': 5e3, #amplitude
    'F': 2.0, #frequency
    'G': 0.0, #shift

    'H': 1.0, #1000.0, #constant
    'I': 0.5, #-999.999, #amplitude
    'L': 4.0, #frequency
    'M': 0.0, #shift

    # solutions          u(x) = N + O*sin(P*(x-Q))          v(x) = R + S*exp(-T*x)
    'N': 0.0,
    'O': 1.0,
    'P': 2.0,
    'Q': 0.0,

    'R': 0.0,
    'S': 1.0,
    'T': 2.0
}

leftBC_u  = {'type':'dirichlet', 'value':u_sol(0.0, kwargs)}
leftBC_v  = {'type':'dirichlet', 'value':v_sol(0.0, kwargs)}
rightBC_u = {'type':'dirichlet', 'value':u_sol(L, kwargs)}
rightBC_v = {'type':'dirichlet', 'value':v_sol(L, kwargs)}

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

solution = sparse.linalg.spsolve(Z, f)
u = solution[:N]
v = solution[N:]
u_rec = mylib.fvm_reconstruct_1D(u)
v_rec = mylib.fvm_reconstruct_1D(v)

#########################################################################################################

u_BGS,  v_BGS,  res_BGS,  nitfinal_BGS  = mylib.method_BlockGaussSeidel( A, B, C, D, f1, f2, nit_max, toll)
u_SPJa, v_SPJa, res_SPJa, nitfinal_SPJa = mylib.method_ShurPartialJacobi(A, B, C, D, f1, f2, nit_max, toll, type='alternate')
u_SPJu, v_SPJu, res_SPJu, nitfinal_SPJu = mylib.method_ShurPartialJacobi(A, B, C, D, f1, f2, nit_max, toll, type='shur_on_u')
u_SPJv, v_SPJv, res_SPJv, nitfinal_SPJv = mylib.method_ShurPartialJacobi(A, B, C, D, f1, f2, nit_max, toll, type='shur_on_v')

SCALING = 1.0
u_Lu,  v_Lu,  res_Lu,  nitfinal_Lu = mylib.method_Loperators(A, B, C, D, f1, f2, nit_max, toll, type='L_on_u', L=-SCALING*B@sparse.csc_matrix(np.diag(1/D.diagonal(), 0))@C)
u_Lv,  v_Lv,  res_Lv,  nitfinal_Lv = mylib.method_Loperators(A, B, C, D, f1, f2, nit_max, toll, type='L_on_v', L=-SCALING*C@sparse.csc_matrix(np.diag(1/A.diagonal(), 0))@B)

u_2Lu, v_2Lu, res_2Lu, nitfinal_Lu = mylib.method_Loperators(A, B, C, D, f1, f2, nit_max, toll, type='L_on_u', L=-SCALING*sparse.csc_matrix(np.diag(B.diagonal(), 0))@sparse.csc_matrix(np.diag(1/D.diagonal(), 0))@C)
u_2Lv, v_2Lv, res_2Lv, nitfinal_Lv = mylib.method_Loperators(A, B, C, D, f1, f2, nit_max, toll, type='L_on_v', L=-SCALING*sparse.csc_matrix(np.diag(C.diagonal(), 0))@sparse.csc_matrix(np.diag(1/A.diagonal(), 0))@B)


Z = Z.toarray()
mats = [A.toarray(), B.toarray(), C.toarray(), D.toarray()]
mats_names = ['A','B','C','D']

print("cond(Z) =", np.linalg.cond(Z))
print(f"monolithic error = {np.linalg.norm(f-np.dot(Z, solution))}")
print(leftBC_u)
print(rightBC_u)
print(leftBC_v)
print(rightBC_v)

linewidth = 2
markersize = 6
markevery = 1

plt.figure(figsize=(8,7), dpi=100)
nn = 5

plt.semilogy(res_BGS,  'k',  linestyle='dashed',  label=r"$BGS$",  linewidth=linewidth,   markersize=markersize,   markevery=markevery+1)
plt.semilogy(res_SPJa, 'ro', linestyle='solid',   label=r"$SPJa$", linewidth=linewidth,   markersize=markersize,   markevery=markevery+2)
plt.semilogy(res_SPJu, 'bo', linestyle='solid',   label=r"$SPJu$", linewidth=linewidth,   markersize=markersize,   markevery=markevery+3)
plt.semilogy(res_SPJv, 'go', linestyle='solid',   label=r"$SPJv$", linewidth=linewidth,   markersize=markersize,   markevery=markevery+3)
plt.semilogy(res_Lu,   'r+', linestyle='solid',   label=r"$Lu$",   linewidth=linewidth,   markersize=markersize,   markevery=markevery+4)
plt.semilogy(res_Lv,   'b+', linestyle='solid',   label=r"$Lv$",   linewidth=linewidth,   markersize=markersize,   markevery=markevery+5)
plt.semilogy(res_2Lu,  'rx', linestyle='solid',   label=r"$2Lu$",   linewidth=linewidth,   markersize=markersize,   markevery=markevery+4)
plt.semilogy(res_2Lv,  'bx', linestyle='solid',   label=r"$2Lv$",   linewidth=linewidth,   markersize=markersize,   markevery=markevery+5)

plt.legend()
plt.grid()
plt.xlabel('iterations')
plt.ylabel('algebraic residual')
plt.title(r'$\beta = ${:.1e}'.format(kwargs['beta']))
# plt.xlim(right=35)

# plt.savefig('figures/dualporosity_aresidual.pdf')
plt.show()