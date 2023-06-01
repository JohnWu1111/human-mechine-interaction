clear;
% close all;
clc;
format long
tic;

dt = 0.01;
T_max = 100;
t = 0:dt:T_max;
nt = length(t);
phi = zeros(3,nt);
phi(1,1) = 1;

for i = 2:nt
    phi(:,i) = myrunge(phi(:,i-1),dt);
end

ob = phi(2,:);

figure;
plot(t,ob)

toc;

function y = myrunge(x,dt)
    c1 = f(x);
    c2 = f(x + c1 * dt / 2);
    c3 = f(x + c2 * dt / 2);
    c4 = f(x + c3 * dt);
    y = x + dt * (c1 + 2 * c2 + 2 * c3 + c4) / 6;
end

function y = f(x)
    y = zeros(3,1);
    eta1 = 0.2;
    eta2 = 0.2;
    gamma1 = 0.1;
    gamma2 = 0.1;
    y(1) = -eta1*(x(1)+gamma1);
    y(2) = eta1*(x(1)+gamma1) - eta2*(x(2)+gamma2);
    y(3) = eta2*(x(2)+gamma2);
end

