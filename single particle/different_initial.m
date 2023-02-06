clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

figure;

dt = 1;
T = 100;
t = 0:dt:T;
nt = length(t);

Sz_all = [0.6 0.4 0.6 0.4];
eta_all = [-sqrt(1.6/0.4) -sqrt(1.4/0.6) sqrt(1.6/0.4) sqrt(1.4/0.6)];
% Sz_all = (eta_all.^2-1)./(1+eta_all.^2);
L = 2;
% initial_all = [2 4 6 12 14 16];
% L = 16;
% L_it = 50;
K= -1;
neta = length(eta_all);
mu = [-0.5 0];
% mu_A = 2;
% mu = mu_A*(2*rand(1,L)-1);
% mu(L_it) = 0;
Tij = gen_H(1,L);

% phi_0 = zeros(L,1);
% phi_0(L_it) = 1;

nit_final = zeros(1,neta);
peak_final = zeros(1,neta);

target = zeros(nt,neta+1);
target(:,1) = t';

% exact ED %%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:neta
    eta = eta_all(n);

    phi = zeros(L,1);
    phi(2) = sqrt(1/(1+abs(eta)^2));
    phi(1) = eta*phi(2);
    test = sum(abs(phi).^2);
%     phi(L_it) = 1;

    nit = zeros(L,nt);
    nit0 = abs(phi).^2;
    nit(:,1) = abs(phi).^2;
    phit = zeros(L,nt);

    pos_mean = zeros(nt,1);
    pos_mean(1) = wmean((1:L)',abs(phi).^2,1);

    Et = zeros(1,nt);
    H = Tij + diag(mu) + K*diag(nit(:,1));
    Et(1) = real(phi'*H*phi);
    target(1,n+1) = Et(1);

    for i = 2:nt
        H = Tij + diag(mu) + K*diag(nit(:,i-1));
        %     phi = expm(-1i*H*dt)*phi;
        [V,D] = eig(H);
        e = diag(D);
        trans = V'*phi;
        phi = V*(exp(-1i*e*dt).*trans);
        nit(:,i) = abs(phi).^2;
        phit(:,i) = phi;
        Et(i) = real(phi'*H*phi);
        pos_mean(i) = wmean((1:L)',nit(:,i),1);
        target(i,n+1) = Et(i);
    end

%     nit_final(n) = mean(nit(L_it,floor(0.9*nt):end));
%     [~,I] = max(nit(:,end));
%     peak_final(n) = I;

    plot(t,Et,'LineWidth',2)
    hold on
end

xlabel('time','FontSize',14)
ylabel('mean position','FontSize',14)

le = cell(1, neta);
for i = 1:neta
    le{i} = strcat('\eta_0 = ', num2str(eta_all(i)));
end
legend(le)

% figure;
% plot(eta_all,peak_final)
% xlabel('eta','FontSize',14)
% ylabel('<n_{it} final>','FontSize',14)

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