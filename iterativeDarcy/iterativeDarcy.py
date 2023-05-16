import numpy as np

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

# Parameters
n = 50
m = 10
max_it = 20
tol = 1e-3
beta = 0.1  # parameter that controls the coupling
log = True

# Saddle point
M = randpd(n)
M = beta * M + randtridiagpd(n)
Bt = np.random.randn(n, m)
B = Bt.T
D = np.zeros((m, m))
omega = -1e-5

if np.linalg.matrix_rank(B)!=m:
    print("Warning!! INF-SUP not satisfied!")
    print(np.linalg.matrix_rank(B), m)
    # exit()

# # # Generic coupling
# M = randpd(n)
# M = M + randtridiagpd(n)
# D = -np.random.randn(m,m)
# D = -randtridiagpd(m)/beta
# Bt = beta*np.random.randn(n, m)
# B = beta*np.random.randn(m, n)
# B = Bt.T
# omega = 0

# linear system
A = np.block([[M, Bt], [B, D]])
f = np.random.rand(n + m, 1)
x = np.linalg.solve(A, f)
print("Condition number:", np.linalg.cond(A))


# Iterative method
Mhm1 = np.linalg.inv(np.diag(np.diag(M)))
S = np.dot(B, np.dot(Mhm1, Bt))
Q = np.block([[M, np.zeros((n, m))], [B, -np.eye(m)]])  # Uzawa
R = np.block([[np.eye(n), np.dot(np.linalg.inv(M), Bt)], [np.zeros((m, n)), S-D]])  # Peters paper

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
    alpha = min(max(alpha,0.1),1.0)
    # print('alpha ', alpha)
    xk = x0 - alpha*z
    r = A.dot(xk) - f
    # print(np.linalg.norm(r))
    if np.linalg.norm(r) < tol:
        print("Converged in ", i, " iterations")
        break
    z0 = z
    z = np.linalg.solve(Q, r)
    beta = np.dot(z.flatten(),z0.flatten())/np.dot(z0.flatten(),z0.flatten())
    # print('beta ', beta)
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
    alpha = min(max(alpha,0.1),1.0)
    if log:
        print('alpha ', alpha)
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
    # print('beta ', beta)
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
    alpha = min(max(alpha,0.1),1.0)
    if log:
        print('alpha ', alpha)
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
    alpha = min(max(alpha,0.1),1.0)
    if log:
        print('alpha ', alpha)
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
    # print('beta ', beta)
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
    alpha = min(max(alpha,0.1),1.0)
    if log:
        print('alpha ', alpha)
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
# New idea iterative scheme
print("New idea with alpha and beta")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
z = np.dot(Rm1h,np.linalg.solve(Q, r))
z0 = z.copy()
for i in range(max_it):
    alpha = np.dot(z.flatten(),r.flatten())/np.dot(z.flatten(),A.dot(z).flatten())
    alpha = min(max(alpha,0.1),1.0)
    if log:
        print('alpha ', alpha)
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
    # print('beta ', beta)
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
# Preconditioned Uzawa
print("Preconditioned Uzawa with alpha and beta")
x0 = np.zeros((n + m, 1))
r = A.dot(x0) - f
z = np.linalg.solve(Q, r)
z0 = z.copy()
for i in range(max_it):
    alpha = np.dot(z.flatten(),r.flatten())/np.dot(z.flatten(),A.dot(z).flatten())
    alpha = min(max(alpha,0.1),1.0)
    if log:
        print('alpha ', alpha)
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
    # print('beta ', beta)
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
p = np.linalg.solve(Q, r)
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