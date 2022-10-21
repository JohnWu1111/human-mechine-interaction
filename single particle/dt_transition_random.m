clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

% dt_all = 0.1:0.1:2;
% dt_all = [0.2 0.4 0.6 0.8 2 3 5 8];
% dt_all = [0.2 0.4];
dt_all = 0.5;
% dt_all = [0.3 0.5 0.8 1 2];
ldt = length(dt_all);
L = 20;
K = -1;
step_max = 1e6;
% mu = [0 0.5];
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
trial_num = 10;
Tij = gen_H(1,L);
H1 = Tij + diag(mu);
% 
% load('mu_L10_K1_No5.mat')

nit_mean_final = zeros(trial_num,ldt);
nit_std_final = zeros(trial_num,ldt);
step_store = zeros(trial_num,ldt);
osc_store = zeros(trial_num,ldt);
x = (1:L)';

% exact ED %%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:trial_num
    phi_0 = rand(L,1);
    phi_0 = phi_0./sqrt(sum(abs(phi_0).^2));
    nit_0 = abs(phi_0).^2;

    for n = 1:ldt
        dt = dt_all(n);
        phi = phi_0;
        nit = nit_0;

        step = 1;
        while true
            H = H1 + K*diag(nit);
            [V,D] = eig(H);
            e = diag(D);
            trans = V'*phi;
            fact = exp(-1i*e*dt);
            temp = fact.*trans;
            phi = V*temp;
            nit_new = abs(phi).^2;
            step = step + 1;

            if step > step_max
                osc_store(j,n) = 1;
                break;
            end
            
            judge = abs(sum(x.*(nit_new -nit)));
            if judge < 1e-8
                break
            end
            nit = nit_new;
        end

        nit_final = nit(:,end);
        nit_mean_final(j,n) = wmean(1:L,nit_final',1);
        nit_std_final(j,n) = sqrt(wmean(((1:L)-nit_mean_final(j,n)).^2,nit_final',1));
        step_store(j,n) = step;
    end
end

nit_mean_mean = mean(nit_mean_final);
nit_mean_std = std(nit_mean_final);
nit_std_mean = mean(nit_std_final);
nit_std_std = std(nit_std_final);
mean_step = mean(step_store);

figure;
errorbar(dt_all,nit_mean_mean,nit_mean_std,'o')
set(gca,'XScale','log')

figure;
errorbar(dt_all,nit_std_mean,nit_std_std,'o')
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