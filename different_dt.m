clear;
% close all;
clc;
format long
tic;

% myseed = 1;
% rng(myseed)

figure;

dt_all = [0.001 0.2 0.5 1];
ldt = length(dt_all);
T = 200;
L = 3;
K = -1;
% mu = [-0.5 0];
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
Tij = gen_H(1,L);

% phi_0 = zeros(L,1);
% dn = -0.259;
% phi_0(1) = sqrt((1+dn)/2);
% phi_0(2) = sqrt((1-dn)/2);

phi_0 = 2*rand(L,1)-1;
phi_0 = phi_0./sqrt(sum(abs(phi_0).^2));


% exact ED %%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:ldt
    dt = dt_all(n);
    t = 0:dt:T;
    nt = length(t);

    phi = phi_0;

    nit = zeros(L,nt);
    nit0 = abs(phi).^2;
    nit(:,1) = abs(phi).^2;
    phit = zeros(L,nt);

    Et = zeros(1,nt);
    H = Tij + diag(mu) + K*diag(nit(:,1));
    Et(1) = phi'*H*phi;

    for i = 2:nt
        H = Tij + diag(mu) + K*diag(nit(:,i-1));
        %     phi = expm(-1i*H*dt)*phi;
        [V,D] = eig(H);
        e = diag(D);
        trans = V'*phi;
        phi = V*(exp(-1i*e*dt).*trans);
        nit(:,i) = abs(phi).^2;
        phit(:,i) = phi;
        Et(i) = phi'*H*phi;
    end

    plot(t,nit(1,:),'LineWidth',2)
    hold on
end

xlabel('time','FontSize',14)
ylabel('<n_1>','FontSize',14)

le = cell(1, ldt);
for i = 1:ldt
    le{i} = strcat('dt = ', num2str(dt_all(i)));
end
legend(le)

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