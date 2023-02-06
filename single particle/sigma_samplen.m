% clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

figure;
% hold on

dt = 1;
T = 100;
t = 0:dt:T;
nt = length(t);

eta = sqrt(0.55/0.45);
L = 2;
% initial_all = [2 4 6 12 14 16];
% L = 16;
% L_it = 50;
K= -3;
mu = [0 0];
sigma = 0.09;
num = 1e4;
% mu_A = 2;
% mu = mu_A*(2*rand(1,L)-1);
% mu(L_it) = 0;
Tij = gen_H(1,L);

% phi_0 = zeros(L,1);
% phi_0(L_it) = 1;

% nit_final = zeros(1,nsigma);
% peak_final = zeros(1,nsigma);

target = zeros(nt,3);
target(:,1) = t';
Sz_mean = zeros(nt,num);

% exact ED %%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:num
    phi = zeros(L,1);
%     phi(2) = sqrt(1/(1+abs(eta)^2));
%     phi(1) = eta*phi(2);
    phi(2) = 0;
    phi(1) = 1;
    test = sum(abs(phi).^2);
%     phi(L_it) = 1;

    nit = zeros(L,nt);
    nit0 = abs(phi).^2;
    nit(:,1) = abs(phi).^2;
    phit = zeros(L,nt);

    pos_mean = zeros(nt,1);
    pos_mean(1) = wmean((1:L)',abs(phi).^2,1);
    
    Sz_mean(1,n) = (nit(1,1) - nit(2,1))/2;

    Et = zeros(1,nt);
    H = Tij + diag(mu) + K*diag(nit(:,1));
    Et(1) = real(phi'*H*phi);

    for i = 2:nt
        H = Tij + diag(mu) + K*diag(nit(:,i-1)).*(1+randn(L,1)*sigma);
        %     phi = expm(-1i*H*dt)*phi;
        [V,D] = eig(H);
        e = diag(D);
        trans = V'*phi;
        phi = V*(exp(-1i*e*dt).*trans);
        nit(:,i) = abs(phi).^2;
        phit(:,i) = phi;
        Et(i) = real(phi'*H*phi);
        pos_mean(i) = wmean((1:L)',nit(:,i),1);
        Sz_mean(i,n) = (nit(1,i) - nit(2,i))/2;
    end

%     nit_final(n) = mean(nit(L_it,floor(0.9*nt):end));
%     [~,I] = max(nit(:,end));
%     peak_final(n) = I;

    
end

temp = mean(Sz_mean,2);
target(:,3) = temp;

semilogy(t,target(:,3),'LineWidth',2)

xlabel('time','FontSize',14)
ylabel('<S_z>','FontSize',14)

% le = cell(1, nsigma);
% for i = 1:nsigma
%     le{i} = strcat('\sigma = ', num2str(sigma_all(i)));
% end
% legend(le)

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