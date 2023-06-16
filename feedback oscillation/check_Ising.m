clear;
% close all;
clc;
format long
tic;

T_max = 100;
dt = 1e-3;
t = 0:dt:T_max;
nt = length(t);
m0 = 0;
J = 1;
beta = 1;
c = 0.1;

phi = zeros(2,nt);
phi(1,1) = 1;
phi(2,1) = 0.1;

for i = 2:nt
    phi(:,i) = myFTCS(phi(:,i-1), dt, J, beta, c, m0);
end

toc;

figure;
plot(t, phi)

function y = myFTCS(x, dt, J, beta, c, m0)
c1 = f(x, J, beta, c, m0);
c2 = f(x+c1*dt/2, J, beta, c, m0);
c3 = f(x+c2*dt/2, J, beta, c, m0);
c4 = f(x+c3*dt, J, beta, c, m0);
y = x + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = f(x, J, beta, c, m0)
y = x;
y(1) = -x(1) + tanh(beta*(J*x(1)+x(2)));
y(2) = -c*(x(1)-m0);
end