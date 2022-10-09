% Starting from a bound state, slowly move the potential and study the
% properties of the bound state, using TSSP method to perform the time
% evolution.

clear;
tic;

global sigma
sigma = 0.05;
% tol = 1e-11;
dt = 0.01;
dx = 0.2;
T = 20;
t = 0:dt:T;
L = 10;
kappa = 10;
x = -L:dx:L-dx;
nt = length(t);
nx = length(x);
phi0 = zeros(nt,nx);
var_x2 = zeros(nt,1);

phi_0 = exp(-x.^2)/sqrt(2*pi);
phi0(1,:) = phi_0./sqrt(sum(abs(phi_0).^2));
var_x2(1) = wmean(x.^2,abs(phi0(1,:)).^2,dx);

miu = zeros(1,nx);
pha2 = zeros(1,nx);
for i = 1:nx
    miu(i) = 2*pi*(-nx/2+i-1)/(2*L);
%     miu(i) = 2*pi*(nx/2+i-1)/(2*L);
    pha2(i) = exp(-1i*dt*miu(i)^2/2);
end
% V(1,:) = f(x);

for i = 2:nt
%     V(i,:) = f(x+v*t(i));
    phi1 = exp(-1i*dt*(x.^2/2+kappa*abs(phi0(i-1,:)).^2)/2).*phi0(i-1,:);
%     phi1f = phi1*exp(-1i*(x'+L)*miu);
%     phi2 = pha2.*phi1f*exp(1i*miu'*(x+L))/nx;
    phi1f = nufft(phi1,x+L,miu/(2*pi));
    phi2 = nufft(pha2.*phi1f,-miu/(2*pi),x+L)/nx;    
    
    phi0(i,:) = exp(-1i*dt*(x.^2/2+kappa*abs(phi0(i-1,:)).^2)/2).*phi2;
    var_x2(i) = wmean(x.^2,abs(phi0(i,:)).^2,dx);
%     phi0(i,:) = phi0(i,:)./sqrt(s*dx);
end
phi00 = abs(phi0).^2;
phi_norm = sum(phi00,2);
% std_phi = sqrt(var_phi);
Ek = sum(abs(phi0(:,2:end)-phi0(:,1:end-1)).^2/dx^2,2)/2;
Ev = sum(x.^2/2.*phi00 + kappa*phi00.^2,2);
E = Ek + Ev;
toc;

figure;
% plot(t,var_x2);
plot(t,E);


function y = f(x)
    global sigma;
    y = -exp(-(x/sigma).^2/2)/(sqrt(2*pi)*sigma);
end

function y = var1(phi,x)
    y = 0;
    len = length(x);
    for i = 1:len
        y = y + x(i)^2*phi(i);
    end
    y = y/len;
end

function y = wmean(x,phi,dx)
    y = 0;
    len = length(x);
    for i = 1:len
        y = y + x(i)*phi(i);
    end
    y = y*dx;
end