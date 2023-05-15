function A = randpd(n)
% Generate a random positive definite matrix of size n x n

% Generate a random matrix
A = randn(n, n);

% Construct a positive definite matrix
A = A * A';

% Add a multiple of the identity matrix to ensure positive definiteness
lambda_min = min(eig(A)); % Minimum eigenvalue
if lambda_min <= 0
    A = A - (2 * lambda_min - 1) * eye(n); % Add a multiple of the identity matrix
end

end