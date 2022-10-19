clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed)

figure;

dt = 0.1;
T = 10000;
t = 0:dt:T;
nt = length(t);

L = 100;
L_it = 50;
K_all = -1:-1:-12;
nK = length(K_all);
% mu = [0 0.5];
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
mu(L_it) = 0;
Tij = gen_H(1,L);

phi_0 = zeros(L,1);
phi_0(L_it) = 1;
% dn = -1;
% phi_0(1) = sqrt((1+dn)/2);
% phi_0(2) = sqrt((1-dn)/2);

% phi_0 = [-0.512012895578133 + 0.508058501746983i;-0.621197245024349 + 0.306322275289446i];

% phi_0 = 2*rand(L,1)-1;
% phi_0 = phi_0./sqrt(sum(abs(phi_0).^2));

nit_final = zeros(1,nK);
peak_final = zeros(1,nK);

% exact ED %%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:nK
    K = K_all(n);

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

    nit_final(n) = mean(nit(L_it,floor(0.9*nt):end));
    [~,I] = max(nit(:,end));
    peak_final(n) = I;

    plot(t,nit(L_it,:),'LineWidth',2)
    hold on
end

xlabel('time','FontSize',14)
ylabel('<n_{it}>','FontSize',14)

le = cell(1, nK);
for i = 1:nK
    le{i} = strcat('K = ', num2str(K_all(i)));
end
legend(le)

figure;
plot(K_all,peak_final)
xlabel('|K|','FontSize',14)
ylabel('<n_{it} final>','FontSize',14)

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