clear variables

% pseudo inverse prodcut rule: https://math.stackexchange.com/questions/2204318/pseudo-inverse-of-product-of-matrices
% pseudio inverse wiki: https://en.wikipedia.org/wiki/Proofs_involving_the_Moore%E2%80%93Penrose_inverse

% J1 = [1 5 1; 1 2 2; 1 8 2];
% J2 = [1 3 1; 1 2 2; 1 2 3; 2 1 3];
% D1 = [1 2 3]';
% D2 = [5 6 7 8]';
% r1 = [0.5 0.3 0.6]';
% r2 = [0.8 0.2 0.1 0.9]';

%J1 = [1 5 4; 7 2 8];
%J2 = [1 3 1; 1 9 2; 8 2 3; 2 1 3; 1 8 3; 5 7 4];
J1 = [1 5 1; 1 2 2; 1 8 3; 1 1 4];
J2 = [1 3 1; 1 2 2; 1 2 3; 2 1 3; 9 2 3; 10 0 2];
D1 = [1 2 5 8; 3 5 1 6]';
D2 = [5 6 7 8 10 3; 11 8 9 2 4 10]';
r1 = [0.5 0.3 0.8 0.2]';
r2 = [0.8 0.2 0.1 0.9 0.6 0.1]';

J = [J1, zeros(size(J1, 1),size(J1, 2)), D1; zeros(size(J2, 1), size(J2, 2)), J2, D2];
r = [r1; r2];

N = J'*J;
nnz(eig(N)) == size(N, 1)

x = J \ r;

% Directly with Matlab's inv
%N1 = J1'*J1;
%N1inv = inv(N1);
%P1 = D1' * J1 * N1inv  * J1';

% Cholesky decomposition
% N1 = J1'*J1;
% L = chol(N1);
% invL = L \ eye(size(J1, 2));
% N1inv = invL*invL';
% P1 = D1'*J1*N1inv*J1';

% QR decomposition
% N1 = J1'*J1; can
% [Q R] = qr(N1);
% N1inv = inv(R)*Q';
% P1 = D1'*J1*N1inv*J1';

% [Q1, Rq1] = qr(J1);
% E  = eye(rank(Rq1));
% n = size(E, 1);
% I1 = zeros(size(Rq1, 1));
% I1(1:n, 1:n) = E;
% P1 = D1' * Q1 * I1 * Q1';

%A1 = D1'*J1*inv(N1)*J1';
%A1 = D1'*J1*pinv(U*S*V')*J1';
R1 = D1'*D1 - P1*D1;
l1 = D1'*r1 - P1*r1;

%%
N2 = J2'*J2;
N2inv = inv(N2);
P2 = D2' * J2 * N2inv  * J2';

%[L2 U2] = lu(J2);
%[Q2 R2] = qr(J2);

% E  = eye(rank(R2));
% n = size(E, 1);
% I2 = zeros(size(R2, 1));
% I2(1:n, 1:n) = E;
% A2 = D2' * Q2 * I2 * Q2';
% %A2 = D2'*(Q2*R2)*inv( (Q2*R2)'*(Q2*R2) )*(Q2*R2)';
% %A2 = D2'*(L2*U2)*inv((L2*U2)'*(L2*U2))*(L2*U2)';
% %A2 = D2'*L2'*L2;
% %A2 = D2'*J2*inv(N2)*J2';
R2 = D2'*D2 - P2*D2;
l2 = D2'*r2 - P2*r2;

%%
R = R1 + R2;
l = l1 + l2;
z = R \ l;

x1 = N1inv*(J1'*r1 - J1'*D1*z);
x2 = N2inv*(J2'*r2 - J2'*D2*z);

xd = [x1; x2; z];
norm(x - xd)

%%


