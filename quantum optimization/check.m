clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

L = 3;
len = 2^L;
M = 100;
dt = 0.2;
beta = zeros(M+1, 1);
phi = zeros(len,M+1);
E = zeros(M+1, 1);
dE = zeros(M, 1);
sigmaz = [1 0; 0 -1];
sigmax = [0 1; 1 0];
sigmay = [0 -1i; 1i, 0];
I2 = eye(2);
Hp = -eye(len) + (kron(kron(sigmaz, sigmaz), I2) + kron(kron(I2, sigmaz), sigmaz))/2;
Hd = kron(kron(sigmax, I2), I2) + kron(kron(I2, sigmax), I2) + kron(kron(I2, I2), sigmax);

C_Hpd = commute(Hd, Hp);
C_Hpd0 = imag(kron(kron(sigmay, sigmaz), I2) + kron(kron(sigmaz, sigmay), I2)...
    + kron(kron(I2, sigmay), sigmaz) + + kron(kron(I2, sigmaz), sigmay));

% [V, D] = eig(Hd);
% phi0 = V(:,1);
[V, D] = eig(sigmax);
V1 = V(:,1);
phi0 = kron(kron(V1, V1), V1);
% phi0 = rand(len,1);
% phi0 = phi0./sqrt(sum(phi0.^2));
phi(:,1) = phi0;
beta0 = imag(phi0'*C_Hpd*phi0);
beta(1) = beta0;
H0 = Hp + beta0*Hd;
E(1) = phi0'*Hp*phi0;

ex = zeros(M+1,1);
ex(1) = beta(1)^2*dt;

phi_it = phi0;
for i = 2:M+1
    H = Hp + beta(i-1)*Hd;
    [V, D] = eig(H);
    e = diag(D);
    trans = V'*phi_it;
    phi_it = V*(exp(-1i*e*dt).*trans);
%     beta(i) = -1i*phi_it'*C_Hpd*phi_it;
%     beta(i) = real(-1i*phi_it'*C_Hpd*phi_it);
    beta(i) = imag(phi_it'*C_Hpd*phi_it);
    ex(i) = beta(i)^2*dt;
    phi(:,i) = phi_it;
    E(i) = real(phi_it'*Hp*phi_it);
    dE(i-1) = E(i) - E(i-1);
end

toc;

function y = commute(A, B)
    y = A*B - B*A;
end