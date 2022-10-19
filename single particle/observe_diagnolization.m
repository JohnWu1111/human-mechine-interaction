clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

step = 100;
L = 100;
K = -1;
mu = 2*(2*rand(L,1)-1);
% for i = 1:6
%     ni = rand(L,1);
% end
ni = zeros(L,1);
ni(39) = 1;
nit0 = ni./sum(ni);
nit = zeros(L,step);
nit(:,1) = nit0;
Tij = gen_H(1,L);
E_step = zeros(1,step);

target = 1;
for i = 2:step
    H = Tij + diag(mu) + K*diag(nit(:,i-1));
    [V,D] = eigs(sparse(H),3,'smallestreal');
    nit(:,i) = abs(V(:,target)).^2;
    E_step(i) = V(:,target)'*H*V(:,target);
end

dt = 0.005;
T = 0:dt:100;
nt = length(T);
nit2 = zeros(L,nt);
nit2(:,1) = nit(:,end).*(1 + 0.00.*(2*rand(L,1)-1));
phi = V(:,1);
Et = zeros(1,nt);
Et(1) = phi'*H*phi;
for i = 2:nt
    H = Tij + diag(mu) + K*diag(nit2(:,i-1));
    phi = myrunge(-1i*H,phi,dt);
    nit2(:,i) = abs(phi).^2;
    Et(i) = phi'*H*phi;
end

figure;
% plot(T,nit2(1,:))
plot(T,Et)

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