clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

% figure;
% hold on

dt = 0.1;
T = 1000;
t = 0:dt:T;
nt = length(t);

eta = sqrt(0.55/0.45);
L = 2;
% initial_all = [2 4 6 12 14 16];
% L = 16;
% L_it = 50;
K= -4;
mu = [0 0];
sigma_s = 0.1:0.01:0.15;
ns = length(sigma_s);
num = 1e3;
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

count = zeros(ns,1);
saturate = zeros(ns,1);

H0 = Tij + diag(mu);

% exact ED %%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:ns
    sigma = sigma_s(m);
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
        H = H0 + K*diag(nit(:,1));
        Et(1) = real(phi'*H*phi);

        for i = 2:nt
%             temp = 1+randn*sigma*[1;-1];
%             nit_real = [nit(1,i-1)*temp(1);nit(2,i-1)*temp(2)];
%             if nit(1,i-1)*temp(1) > 1
%                 saturate(m) = saturate(m) + 1;
%                 nit_real = [1;0];
%             elseif nit(1,i-1)*temp(1) < 0
%                 saturate(m) = saturate(m) + 1;
%                 nit_real = [0;1];   
%             end
            temp = 1+randn(2,1)*sigma;
            nit_real = [nit(1,i-1)*temp(1);nit(2,i-1)*temp(2)];
            if nit_real(1) > 1
                saturate(m) = saturate(m) + 1;
                nit(1,i-1) = 1;
            elseif nit_real(1) < 0
                saturate(m) = saturate(m) + 1;
                nit(1,i-1) = 0;
            end
            if nit_real(2) > 1
                saturate(m) = saturate(m) + 1;
                nit(2,i-1) = 1;
            elseif nit_real(2) < 0
                saturate(m) = saturate(m) + 1;
                nit(2,i-1) = 0;
            end
                %             nit_real = nit(:,i-1).*temp;

            H = H0 + K*diag(nit_real);
            %     phi = expm(-1i*H*dt)*phi;
            %         [V,D] = eig(H);
            %         e = diag(D);
            %         trans = V'*phi;
            %         phi = V*(exp(-1i*e*dt).*trans);
            a = H(1,1);
            b = H(2,2);
            genhao = sqrt(a^2+b^2-2*a*b+4);
            exp_p = exp(-0.5*(-a-b+genhao)*1i*dt);
            exp_m = exp(-0.5*(-a-b-genhao)*1i*dt);
            bagenhao = (b-a)/genhao;
            offdiag = (exp_p-exp_m)/genhao;
            h1 = ((bagenhao+1)*exp_p-(bagenhao-1)*exp_m)/2;
            h2 = ((1-bagenhao)*exp_p-(-1-bagenhao)*exp_m)/2;
            %         expH = [h1 offdiag;
            %             offdiag h2];
            expH = zeros(2,2);
            expH(1,1) = h1;
            expH(2,2) = h2;
            expH(2,1) = offdiag;
            expH(1,2) = offdiag;
            phi = expH*phi;
            nit(:,i) = abs(phi).^2;
            phit(:,i) = phi;
            Sz_mean(i,n) = (nit(1,i) - nit(2,i))/2;
            if nit(2,i)-nit(1,i) > 0
                count(m) = count(m) + 1;
                break;
            end
        end

    end
end

% temp = mean(Sz_mean,2);
% target(:,3) = temp;
% 
% semilogy(t,target(:,3),'LineWidth',2)
% 
% xlabel('time','FontSize',14)
% ylabel('<S_z>','FontSize',14)

figure;
loglog(sigma_s,count)

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