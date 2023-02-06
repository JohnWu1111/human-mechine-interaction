clear;
% close all;
clc;
format long
tic;

myseed = 5;
rng(myseed)

step = 1000;
L = 50;
K = -1;
mu = 0*(2*rand(L,1)-1);
L_it = 25;
% for i = 1:6
%     ni = rand(L,1);
% end
ni = zeros(L,1);
ni(L_it) = 1;
nit0 = ni./sum(ni);
nit = zeros(L,step);
nit(:,1) = nit0;
Tij = gen_H(1,L);
E_step = zeros(1,step);

target = 1;
for i = 2:step
    H = Tij + diag(mu) + K*diag(nit(:,i-1));
    [V,D] = eigs(sparse(H),1,'smallestreal');
    nit(:,i) = abs(V).^2;
    E_step(i) = D;

%     [V,D] = eig(H);
%     %     [~,it] = max(V(L_it,:).^2);
% 
%     VV = V.^2;
%     [~,peak_pos] = max(VV);
%     cand = find(peak_pos == L_it);
%     if isempty(cand)
%         break
%     elseif length(cand) == 1
%         it = cand;
%     else
%         VV_cand = VV(:,cand);
%         [~,it] = max(VV(L_it,:));
%     end
%     
%     nit(:,i) = abs(V(:,it)).^2;
%     E_step(i) = D(it,it);
end

mean_pos = sum((1:L)'.*nit(:,end));

figure;
% plot(1:step,nit(L_it,:))
plot(1:L,nit(:,end)')

dt = 0.5;
T = 0:dt:1000;
nt = length(T);
nit2 = zeros(L,nt);
nit2(:,1) = nit(:,end).*(1 + 0.0.*(2*rand(L,1)-1));
phi = V(:,it);
Et = zeros(1,nt);
Et(1) = phi'*H*phi;
for i = 2:nt
    H = Tij + diag(mu) + K*diag(nit2(:,i-1));
    [V,D] = eig(H);
    e = diag(D);
    trans = V'*phi;
    phi = V*(exp(-1i*e*dt).*trans);
    nit2(:,i) = abs(phi).^2;
    Et(i) = real(phi'*H*phi);
end

figure;
% plot(T,nit2(L_it,:))
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