clear;
% close all;
clc;
format long
tic;

rng(1)

mu = [-0.5 0];
delta = mu(1) - mu(2);
K = 1;
dt = 1;
T = 0:dt:100;
nt = length(T);

sigmaz = [1 0; 0 -1];
sigmax = [0 1; 1 0];

% phi = rand(2,1);
phi = [1 0]';
phi = phi./sqrt(sum(abs(phi).^2));

Et = zeros(1,nt);
phit = zeros(2,nt);
phit(:,1) = phi;
phit2(:,1) = abs(phi).^2;
Sz = zeros(1,nt);

Sz(1) = abs(phi(1))^2 - abs(phi(2))^2;
H = (delta-K*Sz(1))*sigmaz/2 - sigmax;
Et(1) = phi'*H*phi;

% exact ED %%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:nt
    H = (delta-K*Sz(i-1))*sigmaz/2 - sigmax;
%     phi = expm(-1i*H*dt)*phi;
    [V,D] = eig(H);
    e = diag(D);
    trans = V'*phi;
    phi = V*(exp(-1i*e*dt).*trans);
    phit(:,i) = phi;
    phit2(:,i) = abs(phi).^2;
    Sz(i) = abs(phi(1))^2 - abs(phi(2))^2;
    Et(i) = phi'*H*phi;
end

figure
set(gcf, 'position', [250 70 1400 900]);

subplot(2,2,1)
plot(T,Sz)
xlabel('T')
ylabel('Sz')

subplot(2,2,2)
plot(T,Et)
xlabel('T')
ylabel('energy')

subplot(2,2,3)
plot(1:2,mu)
xlabel('position')
ylabel('random field')

subplot(2,2,4)
plot(T,phit2(1,:))
xlabel('T')
ylabel('')

toc;