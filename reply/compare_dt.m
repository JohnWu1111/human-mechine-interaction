clear;
% close all;
clc;
format long
tic;

% Definition of parameters
L = 1; %size
dt = 1e-4;
T = 100;
t = 0:dt:T;
nt = length(t);
J = 1;
delta = 0.5;
K = 1;

S = L/2;

S_z = zeros(L+1,1);
S_p = zeros(L+1);
S_m = zeros(L+1);

% construction of matrice
for m = 1:L+1
    S_z(m) = L/2 - (m-1);
end

for m = 1:L
    S_p(m,m+1) = sqrt(S*(S+1)-S_z(m+1)*(S_z(m+1)+1));
    S_m(m+1,m) = sqrt(S*(S+1)-S_z(m)*(S_z(m)-1));
end

S_x = (S_p + S_m)/2;

phi0 = zeros(L+1,1);
phi0(1) = 1;

Et = zeros(nt,1);
Sz_order = zeros(nt,1);
Sz_order(1) = S_z'*phi0.^2/(L);
H = -2*(K*Sz_order(1)+delta)*diag(S_z) - 2*J*S_x - eye(2)/2;
Et(1) = phi0'*H*phi0;

Et_real = zeros(nt,1);
Et_real(1) = Et(1) + K*L*Sz_order(1)^2;

phi = phi0;
for i = 2:nt
    H = -2*(K*Sz_order(i-1)+delta)*diag(S_z) - 2*J*S_x - eye(2)/2;

    % time revolution
    [V,D] = eig(H);
    e = diag(D);
    temp = V'*phi;
    trans = exp(-1i*e*dt);
    phi1 = temp.*trans;
    phi = V*phi1;

    phit2 = abs(phi).^2;
    Sz_order(i) = S_z'*phit2/(L);
    Et(i) = real(phi'*H*phi);
    Et_real(i) = Et(i) + K*L*Sz_order(i)^2;
end

toc;

figure
% set(gcf, 'position', [100 70 1700 900]);
% subplot(1,2,1)
% % plot(t,Sz_order)
% plot(t,Et_real)
% 
% subplot(1,2,2)
% plot(t,Et)
plot(t,Et_real+(1.2546306-1.0827574),t,Et+(1.2546306-1.0827574))
legend('E_1','E_2')
axis([0 100 -1.4 -1])
xlabel('t')
ylabel('E_n')


