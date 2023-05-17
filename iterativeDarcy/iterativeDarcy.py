import numpy as np
from scipy import sparse

def randpd(n):
    # Generate a random positive definite matrix of size n x n

    # Generate a random matrix
    A = np.random.randn(n, n)

    # Construct a positive definite matrix
    A = np.dot(A, A.T)

    # Add a multiple of the identity matrix to ensure positive definiteness
    lambda_min = np.min(np.linalg.eigvals(A)) # Minimum eigenvalue
    if lambda_min <= 0:
        A = A - (2 * lambda_min - 1) * np.eye(n) # Add a multiple of the identity matrix

    return A

def randtridiagpd(n):
    # Generate a random tridiagonal positive definite matrix of size n x n

    # Generate random vectors for the main diagonal and the off-diagonals
    main_diag = np.random.randn(n)
    off_diag = np.random.randn(n-1)

    # Construct the tridiagonal matrix
    A = np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)

    # Add a multiple of the identity matrix to ensure positive definiteness
    lambda_min = np.min(np.linalg.eigvals(A)) # Minimum eigenvalue
    if lambda_min <= 0:
        A = A - (2 * lambda_min - 1) * np.eye(n) # Add a multiple of the identity matrix

    return A

def power_iteration(matrix, num_iterations=100):
    n = matrix.shape[0]
    x = np.ones(n)
    for _ in range(num_iterations):
        x = matrix.dot(x)
        x /= np.linalg.norm(x)
    eigenvalue = x.dot(matrix.dot(x))
    return eigenvalue

def count_eigenvalues(matrix):
    eigenvalues, _ = np.linalg.eig(matrix)
    eigenvalues_real = np.real(eigenvalues)
    positive_eigenvalues = np.sum(eigenvalues_real > 0)
    negative_eigenvalues = np.sum(eigenvalues_real < 0)
    largest_positive_eigenvalue = np.max(eigenvalues_real[eigenvalues_real > 0]) if positive_eigenvalues > 0 else None
    largest_negative_eigenvalue = np.min(eigenvalues_real[eigenvalues_real < 0]) if negative_eigenvalues > 0 else None
    condition_number = np.linalg.cond(matrix)
    
    print("Number of positive eigenvalues:", positive_eigenvalues)
    print("Number of negative eigenvalues:", negative_eigenvalues)
    print("Largest positive eigenvalue:", largest_positive_eigenvalue)
    print("Largest negative eigenvalue:", largest_negative_eigenvalue)
    print("Condition number:", condition_number)

# Parameters
n = 200
m = 100
max_it = 200
tol = 1e-3
beta = 1  # parameter that controls the coupling
log = False
maxalpha = 10
minalpha = 0.1
maxbeta = 10
minbeta = -10

# Saddle point
# M = randpd(n)
# M = beta * M + randtridiagpd(n)
# Bt = np.random.randn(n, m)
# B = Bt.T
# D = np.zeros((m, m))
# omega = -1e-5

# if np.linalg.matrix_rank(B)!=m:  ## double check!
#     print("Warning!! INF-SUP not satisfied!")
#     print(np.linalg.matrix_rank(B), m)
#     # exit()

# # # Generic coupling
M = randpd(n)
M = randtridiagpd(n)
# D = -np.random.randn(m,m)
D = randtridiagpd(m)/beta
Bt = beta*np.random.randn(n, m)
B = beta*np.random.randn(m, n)
B = -Bt.T
omega = 0


# linear system
print("Assemble and solve block-coupled system A")
A = np.block([[M, Bt], [B, D]])
f = np.random.rand(n + m, 1)
x = np.linalg.solve(A, f)

# Count eigenvalues
count_eigenvalues(A)

# Assemble Schur-based factorisations for Peters' method
print("\n")
print("Assemble and analyse the Schur complements")
Mhm1 = np.linalg.inv(np.diag(np.diag(M)))
S = np.dot(B, np.dot(Mhm1, Bt))    # approximate Schur
Q = np.block([[M, np.zeros((n, m))], [B, -np.eye(m)]])  # Uzawa
R = np.block([[np.eye(n), np.dot(np.linalg.inv(M), Bt)], [np.zeros((m, n)), S-D]])  # Peters paper

# convergence of Peters' modified Uzawa method
Sexact = np.dot(B, np.dot(np.linalg.inv(M), Bt))    # approximate Schur
rho = power_iteration(np.eye(m)-np.dot(np.linalg.inv(S-D),Sexact-D))
print("Spectral radius of Peters' method", rho)
if abs(rho)<1:
    print("Peters' method will converge!")
else:
    print("Peters' method will not converge, check methods with alpha and beta")

print("\n")

# Count eigenvalues
print("Exact Schur complement")
count_eigenvalues(Sexact)

print("\n")

# Count eigenvalues
print("Approximate Schur complement")
count_eigenvalues(S)

# # scipy version
# sA = sparse.csr_matrix(A)
# x = sparse.linalg.spsolve(sA,f)

print("\n")

########################
########################
#  Uzawa
print(" Standard Uzawa with alpha and beta")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
z = np.linalg.solve(Q, r)
z0 = z.copy()
for i in range(max_it):
    alpha = np.dot(z.flatten(),r.flatten())/np.dot(z.flatten(),A.dot(z).flatten())
    if log:
        print('alpha ', alpha)
    alpha = min(max(alpha,minalpha),maxalpha)
    xk = x0 - alpha*z
    r = A.dot(xk) - f
    # print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    z0 = z
    z = np.linalg.solve(Q, r)
    beta = np.dot(z.flatten(),z0.flatten())/np.dot(z0.flatten(),z0.flatten())
    if log:
        print('beta ', beta)
    beta = min(max(beta,minbeta),maxbeta)
    z = z - beta*z0
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")

########################
########################
#  Uzawa
print(" Standard Uzawa ")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
for i in range(max_it):
    xk = x0 - np.linalg.solve(Q, r)
    r = A.dot(xk) - f
    # print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")

########################
########################
# Peters iterative scheme
print("Peters with alpha and beta")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
z = np.linalg.solve(R, np.linalg.solve(Q, r))
z0 = z.copy()
for i in range(max_it):
    alpha = np.dot(z.flatten(),r.flatten())/np.dot(z.flatten(),A.dot(z).flatten())
    if log:
        print('alpha ', alpha)
    alpha = min(max(alpha,minalpha),maxalpha)
    xk = x0 - alpha*z
    r = A.dot(xk) - f
    if log:
        print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    z0 = z
    z = np.linalg.solve(R, np.linalg.solve(Q, r))
    beta = np.dot(z.flatten(),z0.flatten())/np.dot(z0.flatten(),z0.flatten())
    if log:
        print('beta ', beta)
    beta = min(max(beta,minbeta),maxbeta)
    z = z - beta*z0
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")

########################
########################
# Peters iterative scheme
print("Peters with alpha")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
for i in range(max_it):
    z = np.linalg.solve(R, np.linalg.solve(Q, r))
    alpha = np.dot(z.flatten(),r.flatten())/np.dot(z.flatten(),A.dot(z).flatten())
    if log:
        print('alpha ', alpha)
    alpha = min(max(alpha,minalpha),maxalpha)
    xk = x0 - alpha*z
    r = A.dot(xk) - f
    if log:
        print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")

########################
########################
# Peters iterative scheme
print("Peters with alpha and beta in A norm")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
z = np.linalg.solve(R, np.linalg.solve(Q, r))
z0 = z.copy()
for i in range(max_it):
    alpha = np.dot(z.flatten(),A.dot(r).flatten())/np.dot(A.dot(z).flatten(),A.dot(z).flatten())
    if log:
        print('alpha ', alpha)
    alpha = min(max(alpha,minalpha),maxalpha)
    xk = x0 - alpha*z
    r = A.dot(xk) - f
    if log:
        print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    z0 = z
    z = np.linalg.solve(R, np.linalg.solve(Q, r))
    beta = np.dot(z.flatten(),A.dot(z0).flatten())/np.dot(A.dot(z0).flatten(),A.dot(z0).flatten())
    if log:
        print('beta ', beta)
    beta = min(max(beta,minbeta),maxbeta)
    z = z - beta*z0
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")

########################
########################
# Peters iterative scheme
print("Peters with alpha in A norm")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
for i in range(max_it):
    z = np.linalg.solve(R, np.linalg.solve(Q, r))
    alpha = np.dot(z.flatten(),A.dot(r).flatten())/np.dot(A.dot(z).flatten(),A.dot(z).flatten())
    if log:
        print('alpha ', alpha)
    alpha = min(max(alpha,minalpha),maxalpha)
    xk = x0 - alpha*z
    r = A.dot(xk) - f
    if log:
        print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")

########################
########################
# Peters iterative scheme alpha=1
print("Peters original")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
for i in range(max_it):
    z = np.linalg.solve(R, np.linalg.solve(Q, r))
    xk = x0 - z
    r = A.dot(xk) - f
    if log:
        print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")

########################
########################
# preconditioner Q for new idea and preconditioned Uzawa
Q = np.block([[M, np.zeros((n, m))], [B, -S+D]])  # Prec Uzawa
R = np.block([[np.eye(n), np.dot(np.linalg.inv(M), Bt)], [np.zeros((m, n)), np.eye(m)]])
H = np.eye(m + n) - R
Rm1h = H + np.eye(m + n)

########################
########################
# New idea iterative scheme (seems totally equivalent to Peters)
print("New idea with alpha and beta")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
z = np.dot(Rm1h,np.linalg.solve(Q, r))
z0 = z.copy()
for i in range(max_it):
    alpha = np.dot(z.flatten(),r.flatten())/np.dot(z.flatten(),A.dot(z).flatten())
    if log:
        print('alpha ', alpha)
    alpha = min(max(alpha,minalpha),maxalpha)
    xk = x0 - alpha*z
    r = A.dot(xk) - f
    if log:
        print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    z0 = z
    z = np.dot(Rm1h,np.linalg.solve(Q, r))
    beta = np.dot(z.flatten(),z0.flatten())/np.dot(z0.flatten(),z0.flatten())
    if log:
        print('beta ', beta)
    beta = min(max(beta,minbeta),maxbeta)
    z = z - beta*z0
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")

# ########################
# ########################
# # Preconditioned Uzawa
# print("Preconditioned Uzawa BiCG")
# x0 = np.zeros((n + m, 1))
# xx0 = x0.copy()
# r = A.dot(x0) - f
# rr = r.copy()
# p = np.linalg.solve(Q, r)
# pp = p.copy()
# for i in range(max_it):
#     z = np.linalg.solve(Q, r)
#     alpha = np.dot(z.flatten(),rr.flatten())/np.dot(pp.flatten(),A.dot(p).flatten())
#     print('alpha ', alpha)
#     xk = x0 - min(max(0.5,alpha),1)*p
#     xxk = xx0 - min(max(0.5,alpha),1)*pp
#     r0 = r.copy()
#     rr0 = rr.copy()
#     z0 = z.copy()
#     r = A.dot(xk) - f
#     rr = A.T.dot(xxk) - f
#     if np.linalg.norm(r) < tol:
#         print("Converged in ", i, " iterations")
#         break
#     z = np.linalg.solve(Q, r)
#     beta = np.dot(rr.flatten(),z.flatten())/np.dot(rr0.flatten(),z0.flatten())
#     print('beta ', beta)
#     p = z + beta*p
#     pp = np.linalg.solve(Q.T, rr) + beta*pp
#     x0 = xk
#     xx0 = xxk
# if i==max_it-1:
#     print("Not converged")

# print("\n")    

# ########################
# ########################
# # Preconditioned Uzawa
# print("Preconditioned Uzawa CG")
# x0 = np.zeros((n + m, 1))
# r = A.dot(x0) - f
# p = np.linalg.solve(Q, r)
# for i in range(max_it):
#     z = np.linalg.solve(Q, r)
#     alpha = np.dot(z.flatten(),r.flatten())/np.dot(p.flatten(),A.dot(p).flatten())
#     print('alpha ', alpha)
#     xk = x0 - min(max(0.5,alpha),1)*p
#     r0 = r.copy()
#     z0 = z.copy()
#     r = A.dot(xk) - f
#     if np.linalg.norm(r) < tol:
#         print("Converged in ", i, " iterations")
#         break
#     z = np.linalg.solve(Q, r)
#     beta = np.dot(r.flatten(),z.flatten())/np.dot(r0.flatten(),z0.flatten())
#     print('beta ', beta)
#     p = z + beta*p
#     x0 = xk
# if i==max_it-1:
#     print("Not converged")

# print("\n")    

########################
########################
# Preconditioned Uzawa (this is equivalent to our paper with Florin Schur-based partial jacobi)
print("Preconditioned Uzawa with alpha and beta")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
z = np.linalg.solve(Q, r)
z0 = z.copy()
for i in range(max_it):
    alpha = np.dot(z.flatten(),r.flatten())/np.dot(z.flatten(),A.dot(z).flatten())
    if log:
        print('alpha ', alpha)
    alpha = min(max(alpha,minalpha),maxalpha)
    xk = x0 - alpha*z
    r = A.dot(xk) - f
    if log:
        print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    z0 = z
    z = np.linalg.solve(Q, r)
    beta = np.dot(z.flatten(),z0.flatten())/np.dot(z0.flatten(),z0.flatten())
    if log:
        print('beta ', beta)
    beta = min(max(beta,minbeta),maxbeta)
    z = z - beta*z0
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")    

########################
########################
# Preconditioned Uzawa
print("Preconditioned Uzawa")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
for i in range(max_it):
    xk = x0 - np.linalg.solve(Q, r)
    r = A.dot(xk) - f
    if log:
        print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    x0 = xk
if i==max_it-1:
    print("Not converged")

print("\n")    

# ########################
# ########################
# # Gauss-Seidel
# Q = np.block([[M, np.zeros((n, m))], [B, D+omega*np.eye(m)]])  # standard Gauss Seidel + Uzawa/L-relaxation

# # Iterative scheme comparison
# print("Gauss-Seidel/Uzawa with alpha")
# x0 = np.zeros((n + m, 1))
# r = A.dot(x0) - f
# for i in range(max_it):
#     z = np.linalg.solve(Q, r)
#     alpha = np.dot(z.flatten(),r.flatten())/np.dot(z.flatten(),A.dot(z).flatten())
#     # print('alpha ', alpha)
#     xk = x0 - min(max(0.5,alpha),1)*z
#     r = A.dot(xk) - f
#     if np.linalg.norm(r) < tol:
#         print("Converged in ", i, " iterations")
#         break
#     x0 = xk
# if i==max_it-1:
#     print("Not converged")

# print("\n")

# ########################
# ########################
# # Gauss-Seidel
# Q = np.block([[M, np.zeros((n, m))], [B, D+omega*np.eye(m)]])  # standard Gauss Seidel + Uzawa/L-relaxation

# # Iterative scheme comparison
# print("Gauss-Seidel/Uzawa")
# x0 = np.zeros((n + m, 1))
# r = A.dot(x0) - f
# for i in range(max_it):
#     z = np.linalg.solve(Q, r)
#     xk = x0 - z
#     r = A.dot(xk) - f
#     if np.linalg.norm(r) < tol:
#         print("Converged in ", i, " iterations")
#         break
#     x0 = xk
# if i==max_it-1:
#     print("Not converged")

# print("\n")