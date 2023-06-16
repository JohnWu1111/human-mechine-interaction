clear;
% close all;
clc;
format long
tic;

T_max = 10000;
dt = 10;
t = 0:dt:T_max;
nt = length(t);
Jz = 1;
Jx = 0.1;
omega = 0.01;

sigmax = [0 1; 1 0];
sigmay = 1i*[0 -1; 1 0];
sigmaz = [1 0; 0 -1];
I2 = eye(2);

Et = zeros(nt,1);
phit = zeros(2,nt);
Szt = zeros(nt,1);

phi0 = [1;0];
H0 = -Jz*sigmaz + Jx*sigmax;
Et(1) = phi0'*H0*phi0;
Szt(1) = abs(phi0(1))^2-abs(phi0(2))^2;
phi = phi0;

for i = 2:nt
    H = -Jz*cos(omega*2*pi*t(i))*sigmaz + Jx*sigmax;
    [V,D] = eig(H);
    e = diag(D);
    trans = V'*phi;
    phi = V*(exp(-1i*e*dt).*trans);
    Szt(i) = abs(phi(1))^2-abs(phi(2))^2;
    Et(i) = real(phi'*H*phi);
    phit(:,i) = phi;
end

figure
set(gcf, 'position', [100 70 1700 900]);
subplot(1,2,1)
plot(t,Et)
xlabel('t')
ylabel('E')

subplot(1,2,2)
plot(t,Szt)
xlabel('t')
ylabel('Sz')

toc;