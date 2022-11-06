import numpy as np
from scipy import sparse

A = np.zeros((5,5))
A[2,4] = 3.3
A[0,2] = -1.1

x = np.random.rand(5)



sA = sparse.csc_matrix(A)   # Here's the initialization of the sparse matrix.
sB = 4.0*(sA.copy())


print(np.dot(sA.toarray(), x))


A = np.random.rand(2,5)
print(A)
x = np.random.rand(5)
print(A@x)
