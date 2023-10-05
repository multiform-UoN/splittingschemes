fvScalarMatrix Aeqn(fvm::laplacian(-muu, u));
fvScalarMatrix Ceqn(fvm::laplacian(-mvu, u));
fvScalarMatrix Beqn(fvm::laplacian(-muv, v));
fvScalarMatrix Deqn(fvm::laplacian(-mvv, v));

