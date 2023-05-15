%% parameters
n = 60;
m = 30;
max_it=1000;
tol = 1e-4;

%% saddle point
M = randpd(n);
M = 0.1*M + randtridiagpd(n);
Bt = rand(n,m);
B=Bt';
A=[M,Bt;B,zeros(m,m)];
f=rand(n+m,1);
x=A\f;
disp(['Condition number ' num2str(cond(A))])


%% iterative method
Mhm1 = inv(diag(diag(M)));
% Mhm1 = inv(diag(sum((M))));
% Mhm1 = eye(n) + (eye(n)-M);

%S = B*inv(M)*Bt;
S = B*Mhm1*Bt;

%Q = A;   % exact
%Q = eye(m+n);  % Richardson
Q = [M,0*Bt;B,-eye(m,m)];  % Uzawa
% Q = [M,0*Bt;B,-S];  % Prec Uzawa

%% Peters paper
R = [eye(n),inv(M)*Bt;0*B,S];

%% iterative scheme - Peters
disp('Peters')
x0=zeros(n+m,1);
r = A*x0 - f;
for i=1:max_it
    xk = x0 - R\(Q\r);
    enorm = norm(x-xk);
%     disp(enorm)
    if isnan(enorm) || enorm>norm(r)/tol
        break
    end
    if (enorm<tol)
        disp ('Converged')
        disp(i)
        break
    end
    x0 = xk;
    r = A*x0 - f;
end

%% New idea
Q = [M,0*Bt;B,-S];  % Prec Uzawa
R = [eye(n),inv(M)*Bt;0*B,eye(m)];
H = eye(m+n)-R;
Rm1h= H+ eye(m+n);

%% iterative scheme - new idea
disp('new idea')
x0=zeros(n+m,1);
r = A*x0 - f;
for i=1:max_it
    xk = x0 - Rm1h*(Q\r);
    enorm = norm(x-xk);
%     disp(enorm)
    if isnan(enorm) || enorm>norm(r)/tol
        break
    end
    if (enorm<tol)
        disp ('Converged')
        disp(i)
        break
    end
    x0 = xk;
    r = A*x0 - f;
end


%% iterative scheme compare
disp("preconditioned uzawa")
x0=zeros(n+m,1);
r = A*x0 - f;
for i=1:max_it
    xk = x0 - (Q\r);
    enorm = norm(x-xk);
%     disp(enorm)
    if isnan(enorm) || enorm>norm(r)/tol
        break
    end
    if (enorm<tol)
        disp ('Converged')
        disp(i)
        break
    end
    x0 = xk;
    r = A*x0 - f;
end