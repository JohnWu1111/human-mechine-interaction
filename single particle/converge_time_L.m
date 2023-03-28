clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

dt = 1;
T = 0:dt:10000*dt;
nt = length(T);
L_all = 50:10:100;
nL = length(L_all);
K = -1;
tol = 1e-6;

step_con = zeros(nL,1);

for n = 1:nL
    L = L_all(n);
    L_it = L/2;


    % Tij = gen_H_2(1,sqrt(L),L);
    Tij = gen_H(1,L);

    phi = zeros(L,1);
    phi(L_it) = 1;

    Et = zeros(1,nt);
    nit = zeros(L,nt);
    phit = zeros(L,nt);
    nit0 = abs(phi).^2;
    nit(:,1) = abs(phi).^2;
    phit(:,1) = phi;

    H = Tij + K*diag(nit(:,1));
    Et(1) = real(phi'*H*phi);

    % exact ED %%%%%%%%%%%%%%%%%%%%%%%%

    for i = 2:nt
        H = Tij + K*diag(nit(:,i-1));
        %     phi = expm(-1i*H*dt)*phi;
        [V,D] = eig(H);
        e = diag(D);
        trans = V'*phi;
        phi = V*(exp(-1i*e*dt).*trans);
        nit(:,i) = abs(phi).^2;
        phit(:,i) = abs(phi);
        Et(i) = real(phi'*H*phi);
        if abs(nit(L_it,i) - nit(L_it,i-1)) < tol
            break;
        end
    end
    step_con(n) = i;
end

figure;
plot(L_all,step_con)

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