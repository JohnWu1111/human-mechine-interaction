clear;
% close all;
clc;
format long
tic;

% Definition of parameters
L = 2^6; %size
dt = 1e-3;
T = 100;
t = 0:dt:T;
nt = length(t);
J = 1;
hx1 = 0.3;
hx2 = 0.5;
hz = 1e-5;

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

H0 = diag(-J*(S_z.^2)/(L) + hz*S_z) + hx1*S_x;

[V,~] = eig(H0);
phi0 = V(:,1);

% construction of Hamiltonian
H = diag(-J*(S_z.^2)/(L)) + hx2*S_x;

% time revolution
[V,D] = eig(H);
e = diag(D);
temp = V'*phi0;
trans = exp(-1i*e*t);
phi1 = temp.*trans;
phit = V*phi1;

phit2 = abs(phit).^2;
Sz_order = S_z'* phit2/L;

toc;

figure
plot(t,Sz_order)

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end

function y = kron3(a,b,c)
y = (kron(kron(a,b),c));
end

function y = kron_p(a,b)
la = length(a);
lb = length(b);
y = zeros(la*lb,1);
for i = 1:la
    for j = 1:lb
        y((i-1)*lb+j) = a(i) + b(j);
    end
end
end

function y = kron_p4(a,b,c,d)
y = kron_p(kron_p(kron_p(a,b),c),d);
end