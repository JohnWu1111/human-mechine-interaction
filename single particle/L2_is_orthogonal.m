clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

step = 100;
L = 2;
K = -6;
% mu = 0*(2*rand(L,1)-1);
mu = [0.5 0];
Tij = gen_H(1,L);

L_it = 1;
ni = zeros(L,1);
ni(L_it) = 1;
nit0 = ni./sum(ni);
nit = zeros(L,step);
nit1(:,1) = nit0;
phit1 = zeros(L,step);
phit1(:,1) = sqrt(nit0);
for i = 2:step
    H = Tij + diag(mu) + K*diag(nit1(:,i-1));

    [V,D] = eig(H);

    VV = V.^2;
    [~,peak_pos] = max(VV);
    cand = find(peak_pos == L_it);
    if isempty(cand)
        break
    elseif length(cand) == 1
        it = cand;
    else
        VV_cand = VV(:,cand);
        [~,it] = max(VV(L_it,:));
    end
    
    nit1(:,i) = VV(:,it);
    phit1(:,i) = V(:,it);
end
phi1 = phit1(:,end);

L_it = 2;
ni = zeros(L,1);
ni(L_it) = 1;
nit0 = ni./sum(ni);
nit2 = zeros(L,step);
nit2(:,1) = nit0;
phit2 = zeros(L,step);
phit2(:,1) = sqrt(nit0);
for i = 2:step
    H = Tij + diag(mu) + K*diag(nit2(:,i-1));

    [V,D] = eig(H);

    VV = V.^2;
    [~,peak_pos] = max(VV);
    cand = find(peak_pos == L_it);
    if isempty(cand)
        break
    elseif length(cand) == 1
        it = cand;
    else
        VV_cand = VV(:,cand);
        [~,it] = max(VV(L_it,:));
    end
    
    nit2(:,i) = VV(:,it);
    phit2(:,i) = V(:,it);
end
phi2 = phit2(:,end);

phi1'*phi2

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

function y = myrunge(H,phi,dt)
c1 = H*phi;
c2 = H*(phi+c1.*(dt/2));
c3 = H*(phi+c2.*(dt/2));
c4 = H*(phi+c3.*dt);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = wmean(x,phi,dx)
    y = sum(x.*phi)*dx;
end