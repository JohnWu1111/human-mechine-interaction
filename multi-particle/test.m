clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

dt = 1;
T = 0:dt:1000*dt;
nt = length(T);
L = 2;
J = 5;
L_it = 1;
mu_A = 1;
mu = mu_A*(2*rand(1,L)-1);
K = 3;

nit = zeros(L,nt);
phit = zeros(L,nt);
phi = zeros(L,1);
phi(L_it) = 1;
phit(:,1) = phi;
nit(:,1) = abs(phi).^2;

Tij = gen_H(J,L,K*nit(:,1));
H = Tij + diag(mu);
Et(1) = phi'*H*phi;
pos_mean = zeros(nt,1);
var_x2 = zeros(nt,1);
pos_mean(1) = wmean((1:L)',abs(phi).^2,1);
var_x2(1) = sqrt(wmean(((1:L)'-pos_mean(1)).^2,abs(phi).^2,1));


for i = 2:nt
    Tij = gen_H(J,L,K*nit(:,i-1));
    H = Tij + diag(mu);
%     phi = expm(-1i*H*dt)*phi;
    [V,D] = eig(H);
    e = diag(D);
    trans = V'*phi;
    phi = V*(exp(-1i*e*dt).*trans);
    phit(:,i) = abs(phi);
    nit(:,i) = abs(phi).^2;
    pos_mean(i) = wmean((1:L)',nit(:,i),1);
    var_x2(i) = sqrt(wmean(((1:L)'-pos_mean(i)).^2,nit(:,i),1));
    Et(i) = real(phi'*H*phi);
end

filename = strcat('L = ',num2str(L), ', mu_A = ', num2str(mu_A), ', K = ', num2str(K), ', seed = ', num2str(myseed), ', dt = ', num2str(dt),', L_it = ',num2str(L_it));
figure('Name',filename);
set(gcf, 'position', [250 70 1900 900]);

subplot(2,3,1)
plot(T,var_x2)
xlabel('T')
ylabel('variance')

subplot(2,3,2)
plot(T,nit(L_it,:))
xlabel('T')
ylabel('ni of L_it')

subplot(2,3,3)
plot(T,Et)
xlabel('T')
ylabel('energy')

subplot(2,3,4)
meshc(nit)

subplot(2,3,5)
plot(1:L,mu)
xlabel('position')
ylabel('random field')

dEt = Et(2:end) - Et(1:end-1);
subplot(2,3,6)
% plot(T(floor(nt/2)+1:end),dEt(floor(nt/2):end))
% xlabel('T')
% ylabel('dE')
plot(T,pos_mean)
xlabel('T')
ylabel('pos_mean')

toc;

function Tij = gen_H(J,L,nit)
Tij = zeros(L);
count = 0;
for i = 1:L-1
    s = J - nit(i);
    Tij(i,i+1) = Tij(i,i+1)-s;
    Tij(i+1,i) = Tij(i+1,i)-conj(s);    
    count = count +1;    
end
s = J - nit(L);
% Tij(L,1) = Tij(L,1)-s;
% Tij(1,L) = Tij(1,L)-conj(s);
count = count +1;
end

function y = wmean(x,phi,dx)
    y = sum(x.*phi)*dx;
end