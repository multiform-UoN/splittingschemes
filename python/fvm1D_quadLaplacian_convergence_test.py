from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from scipy.sparse import diags
import scipy.sparse as sparse
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
# plt.style.use('tableau-colorblind10')
# plt.style.use('default')
import mylib


def u_sol(x, kwargs):
    return np.exp(np.sin(x))

def du_soldx(x, kwargs):
    return np.cos(x)*np.exp(np.sin(x))

def d2u_soldx2(x, kwargs):
    return -np.sin(x)*np.exp(np.sin(x)) + np.square(np.cos(x))*np.exp(np.sin(x))


def v_sol(x, kwargs):
    return - np.square(x) + x - 1.0 

def dv_soldx(x, kwargs):
    return - 2.0*x + 1.0

def d2v_soldx2(x, kwargs):
    return - 2.0 + 0.0*x


def fmuu(x, kwargs):
    return kwargs['A'] + kwargs['B']*np.sin(kwargs['C']*x)

def dfmuudx(x, kwargs):
    return kwargs['B']*kwargs['C']*np.cos(kwargs['C']*x)


def fmuv(x, kwargs):
    return kwargs['D'] + kwargs['E']*x

def dfmuvdx(x, kwargs):
    return kwargs['E'] + 0.0*x


def fmvu(x, kwargs):
    return -fmuv(x, kwargs)

def dfmvudx(x, kwargs):
    return -dfmuvdx(x, kwargs)


def fmvv(x, kwargs):
    return kwargs['F'] + kwargs['G']*np.sin(kwargs['H']*x)

def dfmvvdx(x, kwargs): 
    return kwargs['G']*kwargs['H']*np.cos(kwargs['H']*x)


def f_1(x, kwargs):
    return dfmuudx(x, kwargs)*du_soldx(x, kwargs) + fmuu(x, kwargs)*d2u_soldx2(x, kwargs) + dfmuvdx(x, kwargs)*dv_soldx(x, kwargs) + fmuv(x, kwargs)*d2v_soldx2(x, kwargs)

def f_2(x, kwargs):
    return dfmvudx(x, kwargs)*du_soldx(x, kwargs) + fmvu(x, kwargs)*d2u_soldx2(x, kwargs) + dfmvvdx(x, kwargs)*dv_soldx(x, kwargs) + fmvv(x, kwargs)*d2v_soldx2(x, kwargs)









L = np.pi
N = 128
dx = np.divide(L, N)
x = np.linspace(0, L, N+1)
xc = np.linspace(0.5*dx, L-0.5*dx, N)
xfine = np.linspace(0, L, 1000)
nit_max = 200
toll = 1e-12

kwargs = {
    # muu coefficients
    'A': 1, #const
    'B': 0.5, #amplitude
    'C': 4.0, #frequency
    
    # muv=mvu coefficients
    'D': 1e-1, #const
    'E': 0, #slope
    
    # mvv coefficients
    'F': 1e-1, #const
    'G': 0.5e-1, #amplitude
    'H': 2  #frequency
    }

leftBC_u = {
    'type': 'dirichlet', #'type':'neumann',
    'value': u_sol(0.0, kwargs)}

rightBC_u = {
    'type': 'dirichlet', #'type':'neumann',
    'value': u_sol(L, kwargs)}

leftBC_v = {
    'type': 'dirichlet', #'type':'neumann',
    'value': v_sol(0.0, kwargs)}

rightBC_v = {
    'type': 'dirichlet', #'type':'neumann',
    'value': v_sol(L, kwargs)}

A, fuuBC = mylib.fvm_laplacian_1D(fmuu, leftBC_u, rightBC_u, N, L, kwargs)
B, fuvBC = mylib.fvm_laplacian_1D(fmuv, leftBC_v, rightBC_v, N, L, kwargs)
C, fvuBC = mylib.fvm_laplacian_1D(fmvu, leftBC_u, rightBC_u, N, L, kwargs)
D, fvvBC = mylib.fvm_laplacian_1D(fmvv, leftBC_v, rightBC_v, N, L, kwargs)

f1 = f_1(xc, kwargs) + fuuBC + fuvBC
f2 = f_2(xc, kwargs) + fvuBC + fvvBC

Z = np.zeros((2*N,2*N))
Z[:N,:N] = A.toarray()
Z[:N,N:] = B.toarray()
Z[N:,:N] = C.toarray()
Z[N:,N:] = D.toarray()
Z = sparse.csc_matrix(Z)
f = np.concatenate((f1, f2))
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
u_Lu, v_Lu, res_Lu, nitfinal_Lu = mylib.method_Loperators(A, B, C, D, f1, f2, nit_max, toll, type='L_on_u', L=-SCALING*B@sparse.csc_matrix(np.diag(1/D.diagonal(), 0))@C)
u_Lv, v_Lv, res_Lv, nitfinal_Lv = mylib.method_Loperators(A, B, C, D, f1, f2, nit_max, toll, type='L_on_v', L=-SCALING*C@sparse.csc_matrix(np.diag(1/A.diagonal(), 0))@B)
u_2Lu, v_2Lu, res_2Lu, nitfinal_2Lu = mylib.method_Loperators(A, B, C, D, f1, f2, nit_max, toll, type='L_on_u', L=-SCALING*B@sparse.csc_matrix(np.diag(1/D.diagonal(), 0))@C)
u_2Lv, v_2Lv, res_2Lv, nitfinal_2Lv = mylib.method_Loperators(A, B, C, D, f1, f2, nit_max, toll, type='L_on_v', L=-SCALING*C@sparse.csc_matrix(np.diag(1/A.diagonal(), 0))@B)


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
plt.semilogy(res_Lu,   'r+', linestyle='solid',   label=r"$Lu$",   linewidth=linewidth,   markersize=markersize,   markevery=markevery+4)
plt.semilogy(res_Lv,   'b+', linestyle='solid',   label=r"$Lv$",   linewidth=linewidth,   markersize=markersize,   markevery=markevery+5)

plt.legend()
plt.grid()
plt.xlabel('iterations')
plt.ylabel('algebraic residual')
# plt.title(r'$\beta = ${:.1e}'.format(kwargs['beta']))
# plt.xlim(right=35)

# plt.savefig('figures/dualporosity_aresidual.pdf')
plt.show()