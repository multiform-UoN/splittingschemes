fvScalarMatrix Aeqn(fvm::Sp( beta, u) - fvm::laplacian(mu, u));
fvScalarMatrix Ceqn(fvm::Sp(-beta, u));
fvScalarMatrix Beqn(fvm::Sp(-beta, v));
fvScalarMatrix Deqn(fvm::Sp( beta, v) - fvm::laplacian(mv, v));
