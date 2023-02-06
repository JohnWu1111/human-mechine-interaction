clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

figure;

% dt_all = [1e-3 1e-2 1];
dt_all = [1e-5 1e-2 1e-1 1];
% dt_all = [0.1 1];
% dt_all = 1:1:10;
ldt = length(dt_all);
T = 100;
L = 2;
L_it = 1;
K = -1;
mu = [-0.5 0.5];
% mu_A = 2;
% mu = mu_A*(2*rand(1,L)-1);
% mu(L_it) = 0;
Tij = gen_H(1,L);

dt = dt_all(2);
t = 0:dt:T;
nt = length(t);
target = zeros(nt,ldt*2);
target = target + nan;

phi_0 = zeros(L,1);
% phi_0(L_it) = 1;
dn = 1;
phi_0(1) = sqrt((1+dn)/2);
phi_0(2) = sqrt((1-dn)/2);

% phi_0 = [-0.512012895578133 + 0.508058501746983i;-0.621197245024349 + 0.306322275289446i];

% phi_0 = 2*rand(L,1)-1;
% phi_0 = phi_0./sqrt(sum(abs(phi_0).^2));

nit_final = zeros(1,ldt);
pos_mean = cell(1,ldt);

for n = 1:ldt
    dt = dt_all(n);
    t = 0:dt:T;
    nt = length(t);

    phi = phi_0;

    nit = zeros(L,nt);
    nit0 = abs(phi).^2;
    nit(:,1) = abs(phi).^2;
    phit = zeros(L,nt);
    phit(:,1) = phi;
    
    pos_mean{n} = zeros(nt,1);
    pos_mean{n}(1) = (1:L)*nit0;

    Et = zeros(1,nt);
    H = Tij + diag(mu) + K*diag(nit(:,1));
    Et(1) = phi'*H*phi;

    target(1,2*n-1) = t(1);
    target(1,2*n) = Et(1);

    % ED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:nt
        H = Tij + diag(mu) + K*diag(nit(:,i-1));
        %     phi = expm(-1i*H*dt)*phi;
        [V,D] = eig(H);
        e = diag(D);
        trans = V'*phi;
        phi = V*(exp(-1i*e*dt).*trans);
        nit(:,i) = abs(phi).^2;
        phit(:,i) = phi;
        pos_mean{n}(i) = (1:L)*nit(:,i);
        Et(i) = real(phi'*H*phi);
        if n == 1
            if mod(i-1,1000) == 0
                target(1+(i-1)/1000,2*n-1) = t(i);
                target(1+(i-1)/1000,2*n) = Et(i);
            end
        else
            target(i,2*n-1) = t(i);
            target(i,2*n) = Et(i);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % TSSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     miu= -pi + (0:L-1)*2*pi/L;
%     miu = miu';
%     pha2 = exp(-1i*dt*2*cos(miu));
%     for i = 2:nt
%         phi1 = exp(-1i*dt*(mu'+K*abs(phi).^2)/2).*phi;
%         phi1f = nufft(phi1,1:L,miu/(2*pi));
%         phi2 = nufft(pha2.*phi1f,-miu/(2*pi),1:L)/L;
%         phi = exp(-1i*dt*(mu'+K*abs(phi).^2)/2).*phi2;
%     
%         nit(:,i) = abs(phi).^2;
%         phit(:,i) = phi;
%         pos_mean{n}(i) = (1:L)*nit(:,i);
%         Et(i) = real(phi'*H*phi);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nit_final(n) = mean(nit(L_it,floor(0.9*nt):end));

    plot(t,Et,'LineWidth',2)
    hold on
end

% plot(t,-3.42*ones(1,nt),'LineWidth',2)

xlabel('time','FontSize',14)
ylabel('energy','FontSize',14)

le = cell(1, ldt);
for i = 1:ldt
    le{i} = strcat('dt = ', num2str(dt_all(i)));
end
legend(le)

% figure;
% plot(dt_all,nit_final)
% xlabel('dt','FontSize',14)
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