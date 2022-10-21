clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed*1e5)

% dt_all = 0.1:0.1:2;
dt_all = [0.1 0.3 0.5 0.8 1 2 3 5 8];
% dt_all = 0.2;
% dt_all = [0.3 0.5 0.8 1 2];
ldt = length(dt_all);
step = 2000;
L = 4;
K = -2;
% mu = [0 0.5];
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
trial_num = 100;
Tij = gen_H(1,L);
H1 = Tij + diag(mu);


nit_mean_final = zeros(trial_num,ldt);

% exact ED %%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:trial_num
    phi = rand(L,1);
    phi = phi./sqrt(sum(abs(phi).^2));
    nit = abs(phi).^2;

    for n = 1:ldt
        dt = dt_all(n);
        t = 0:dt:dt*step;
        nt = step+1;

        nit = zeros(L,nt);
        nit0 = abs(phi).^2;
        nit(:,1) = abs(phi).^2;
        phit = zeros(L,nt);

        H = H1 + K*diag(nit(:,1));

        for i = 2:nt
            H = H1 + K*diag(nit(:,i-1));
            %     phi = expm(-1i*H*dt)*phi;
            [V,D] = eig(H);
            e = diag(D);
            trans = V'*phi;
            phi = V*(exp(-1i*e*dt).*trans);
            nit(:,i) = abs(phi).^2;
            phit(:,i) = phi;
        end

        nit_final = mean(nit(:,floor(0.9*nt):end),2);
        nit_mean_final(j,n) = wmean(1:L,nit_final',1);
    end
end

transit_mean = mean(nit_mean_final);
transit_std = std(nit_mean_final);

figure;
errorbar(dt_all,transit_mean,transit_std,'o')
set(gca,'XScale','log')


toc;

function Tij = gen_H(s,L)
Tij = zeros(L);
count = 0;
for i = 1:L-1    
    Tij(i,i+1) = Tij(i,i+1)-s;
    Tij(i+1,i) = Tij(i+1,i)-conj(s);    
    count = count +1;    
end
Tij(L,1) = Tij(L,1)-s;
Tij(1,L) = Tij(1,L)-conj(s);
count = count +1;
end

function y = wmean(x,phi,dx)
    y = sum(x.*phi)*dx;
end