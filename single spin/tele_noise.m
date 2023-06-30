clear;
% close all;
clc;
format long
tic;

T_max = 1000;
dt = 1;
t = 0:dt:T_max;
nt = length(t);
Jz = 1;
Jx = 1;
pos = 0.1;

sigmax = [0 1; 1 0];
sigmay = 1i*[0 -1; 1 0];
sigmaz = [1 0; 0 -1];
I2 = eye(2);

Et = zeros(nt,1);
phit = zeros(2,nt);
Szt = zeros(nt,1);
tele = zeros(nt,1);
tele(1) = 1;

phi0 = [1;0];
H0 = -Jz*sigmaz + Jx*sigmax;
Et(1) = phi0'*H0*phi0;
Szt(1) = abs(phi0(1))^2-abs(phi0(2))^2;
phi = phi0;

for i = 2:nt
    H = -Jz*sigmaz + Jx*tele(i-1)*sigmax;
    [V,D] = eig(H);
    e = diag(D);
    trans = V'*phi;
    phi = V*(exp(-1i*e*dt).*trans);
    Szt(i) = abs(phi(1))^2-abs(phi(2))^2;
    Et(i) = real(phi'*H*phi);
    phit(:,i) = phi;
    if rand < pos*dt
        tele(i) = -tele(i-1);
    else
        tele(i) = tele(i-1);
    end
end

tele_f = fft(tele);
tele_f_abs = abs(tele_f);
dw = 2*pi/T_max;
w = 0:dw:2*pi/dt;

figure
set(gcf, 'position', [100 70 1700 900]);
subplot(2,2,1)
plot(t,Szt)
xlabel('t')
ylabel('Sz')

subplot(2,2,2)
plot(t,Et)
xlabel('t')
ylabel('E')

subplot(2,2,3)
plot(t,tele)
xlabel('t')
ylabel('noise')

subplot(2,2,4)
plot(w,tele_f_abs)
xlabel('w')
ylabel('noise_f')

toc;