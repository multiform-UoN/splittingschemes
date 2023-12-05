import scipy.sparse as sparse
import scipy.linalg as linalg
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mylib

kwargs = {
    'A': 0,   #const
    'B': 1.0, #amplitude
    'C': 3.0, #frequency
    'D': 0.1, #shift
    
    'E': 3.0, #const
    'F': 2.9, #amplitude
    'G': 9.0, #frequency
    'H': 0.5  #shift
}

L = 2*np.pi
N = 128
dx = np.divide(L, N)
x = np.linspace(0, L, N+1)
xc = np.linspace(0.5*dx, L-0.5*dx, N)
xfine = np.linspace(0, L, 1000)

def alpha(x, kwargs):
    return 1.0 + 0.999*np.sin(8.0*x)

def sol(x, kwargs):
    return kwargs['A'] + kwargs['B']*np.sin(kwargs['C']*(x-kwargs['D']))

def dsoldx(x, kwargs):
    return kwargs['B']*kwargs['C']*np.cos(kwargs['C']*(x-kwargs['D']))

def d2soldx2(x, kwargs):
    return -kwargs['B']*np.square(kwargs['C'])*np.sin(kwargs['C']*(x-kwargs['D']))

def nu(x, kwargs):
    return kwargs['E'] + kwargs['F']*np.sin(kwargs['G']*(x-kwargs['H']))

def dnudx(x, kwargs):
    return kwargs['F']*kwargs['G']*np.cos(kwargs['G']*(x-kwargs['H']))
   
def f(x, kwargs):
    return alpha(x, kwargs)*sol(x, kwargs) + nu(x, kwargs)*d2soldx2(x, kwargs) + dnudx(x, kwargs)*dsoldx(x, kwargs)
    
leftBC = {
    'type':'dirichlet',
    'value':sol(0.0, kwargs) #set the value 
    }

rightBC = {
    'type':'dirichlet', #'type':'neumann',
    'value':sol(L, kwargs)
    }

A_laplacian, fBC = mylib.fvm_laplacian_1D(nu, leftBC, rightBC, N, L, kwargs)
fTOT = f(xc, kwargs) + fBC
A = sparse.csc_matrix(A_laplacian + sparse.diags([alpha(xc, kwargs)], [0]))
u = sparse.linalg.spsolve(A, fTOT)
A = A.toarray()

print(leftBC)
print(rightBC)
print("cond(A) =", np.linalg.cond(A))
print(f"monolithic error = {np.linalg.norm(fTOT-np.dot(A,u))}")

N = np.diag(1.0/np.linalg.norm(A, axis=1))
AA = N@A
ff = N@fTOT
gg = ff/np.linalg.norm(ff)
idxs = np.where(np.abs(AA.T@ff)>1e-4)
print("cond(AA) =", np.linalg.cond(AA))
fig0, axs = plt.subplots(1, 2, figsize=(24,5))
# axsA = axs[0].matshow(A)
# axsA = axs[0].matshow(AA)
axsA = axs[0].matshow(AA[:,idxs[0]])

# axs[1].set_title('fTOT[cell_index]')
# axsf = axs[1].plot(ff, '.-')

axsf = axs[1].plot(np.abs(AA.T@gg), '-')
axs[1].grid()
fig0.colorbar(axsA, ax=axs[0])

u_pre = np.zeros_like(u)
u_pre[idxs[0]] = np.linalg.lstsq(AA[:,idxs[0]], ff)[0]

print(np.linalg.norm(AA@u_pre-ff))#/np.linalg.norm(u))
print(np.linalg.norm(AA@u-ff))


# fig1, ax1 = plt.subplots(figsize=(20,4), dpi=100)
# ax1.set_xlabel(r'$x$')
# ax1.set_ylabel(r'$u(x)$')
# ax1.plot(x, 0*x, '|', color='k', markersize=6.0)
# ax1.plot(xc, 0*xc, '.', color='k', markersize=2.0)
# ax1.plot(xfine, sol(xfine, kwargs), '-k', label=r'$u(x)$')
# ax1.plot(x, mylib.fvm_reconstruct_1D(u), '.r' , label=r'$u_{reconstruct}$', linewidth=0.8)
# ax1.plot(xc, u, '_y', label=r'$\frac{1}{|cell|}\int u\,dx$', markersize=10)
# ax1.legend(bbox_to_anchor=(1.15, 0.8))
# ax1.grid(color='k', alpha=0.3)


# # fig2, ax2 = plt.subplots(1, 2, figsize=(24,3), dpi=100)
# # ax2[0].plot(xfine, nu(xfine, kwargs))
# # ax2[0].set_xlabel(r'$x$')
# # ax2[0].set_ylabel(r'$\nu$')
# # ax2[1].plot(xfine, alpha(xfine, kwargs))
# # ax2[1].set_xlabel(r'$x$')
# # ax2[1].set_ylabel(r'$\alpha$')
# # ax2[0].grid()
# # ax2[1].grid()


plt.show()