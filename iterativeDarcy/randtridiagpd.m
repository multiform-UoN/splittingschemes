function A = randtridiagpd(n)
% Generate a random tridiagonal positive definite matrix of size n x n

% Generate random vectors for the main diagonal and the off-diagonals
main_diag = randn(n, 1);
off_diag = randn(n-1, 1);

% Construct the tridiagonal matrix
A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

% Add a multiple of the identity matrix to ensure positive definiteness
lambda_min = min(eig(A)); % Minimum eigenvalue
if lambda_min <= 0
    A = A - (2 * lambda_min - 1) * eye(n); % Add a multiple of the identity matrix
end

end
