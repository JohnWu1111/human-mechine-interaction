clear;
% close all;
clc;
format long
tic;

omega = 1;
E = 1;
gamma = 0.1;

rho0 = zeros(4,1);
rho0(4) = 1;
Liou = [0 1i*omega -1i*omega gamma;
    1i*omega 1i*E-gamma/2 0 -1i*omega;
    -1i*omega 0 -1i*E-gamma/2 1i*omega;
    0 -1i*omega 1i*omega -gamma];
dt = 0.01;
t = 0:dt:10;
nt = length(t);
phit = zeros(4,nt);
phit(:,1) = rho0;
H = Liou;

%%%%%%%%%% ED %%%%%%%%%%%%
for i = 2:nt
    [V,D] = eig(H);
    e = diag(D);
    phi = phit(:,i-1);
    trans = V\phi;
    phi = V*(exp(e*dt).*trans);
    phit(:,i) = phi;
end

% for i = 2:nt
%     phit(:,i) = myrunge(Liou,phit(:,i-1),dt);
% end

toc;

check = phit(1,:) + phit(4,:);

figure;
plot(t,real(phit(1,:)))
hold on
plot(t,real(phit(4,:)))

function y = myrunge(H,phi,dt)
c1 = H*phi;
c2 = H*(phi+c1.*(dt/2));
c3 = H*(phi+c2.*(dt/2));
c4 = H*(phi+c3.*dt);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end