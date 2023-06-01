clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

dt0 = 5;
sigma = 0.01;
nt = 1e3;
L = 8;
K = -4;
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
N_s = 1e4;

dn_final = zeros(N_s,1);
pos_mean_final = zeros(N_s,1);
E_final = zeros(N_s,1);
x = (1:L)';

Tij = gen_H(1,L);

parfor n = 1:N_s

    phi = rand(L,1);
    phi = phi./sqrt(sum(abs(phi).^2));

    nit = abs(phi).^2;
    H = Tij + diag(mu) + K*diag(nit);

    % exact ED %%%%%%%%%%%%%%%%%%%%%%%%

    for i = 2:nt
        dt = max(dt0*(1+randn*sigma),0);
        H = Tij + diag(mu) + K*diag(nit);
        [V,D] = eig(H);
        e = diag(D);
        trans = V'*phi;
        phi = V*(exp(-1i*e*dt).*trans);
        nit = abs(phi).^2;
    end

    dn_final(n) = nit(1)-nit(2);
    pos_mean_final(n) = sum(nit.*x);
    E_final(n) = real(phi'*H*phi);
end

figure
histogram(E_final,100,'Normalization','pdf','DisplayStyle','stairs');

toc;

function Tij = gen_H(s,L)
Tij = zeros(L);
count = 0;
for i = 1:L-1
    Tij(i,i+1) = Tij(i,i+1)-s;
    Tij(i+1,i) = Tij(i+1,i)-conj(s);
    count = count +1;
end
% Tij(L,1) = Tij(L,1)-s;
% Tij(1,L) = Tij(1,L)-conj(s);
count = count +1;
end


function y = wmean(x,phi,dx)
y = sum(x.*phi)*dx;
end