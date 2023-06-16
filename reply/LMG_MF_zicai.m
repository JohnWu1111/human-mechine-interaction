clear;
% close all;
clc;
format long
tic;

% Definition of parameters
L = 2^8; %size
dt = 1e-2;
T = 100;
t = 0:dt:T;
nt = length(t);
J = 1;
hx1 = 0.3;
hx2 = 0.5;
hz = 1e-5;

S = L;

S_z = zeros(2*L+1,1);
S_p = zeros(2*L+1);
S_m = zeros(2*L+1);

% construction of matrice
for m = 1:2*L+1
    S_z(m) = L - (m-1);
end

for m = 1:2*L
    S_p(m,m+1) = sqrt(S*(S+1)-S_z(m+1)*(S_z(m+1)+1));
    S_m(m+1,m) = sqrt(S*(S+1)-S_z(m)*(S_z(m)-1));
end

S_x = (S_p + S_m)/2;

H0 = diag(-J*(S_z.^2)/(2*L) + hz*S_z) + hx1*S_x;

[V,~] = eig(H0);
phi0 = V(:,1);

Et = zeros(nt,1);
Et(1) = phi0'*H0*phi0;
Sz_order = zeros(nt,1);
Sz_order(1) = S_z'*phi0.^2/(2*L);

phi = phi0;
for i = 2:nt
    H = -J*Sz_order(i-1)*diag(S_z) + hx2*S_x;

    % time revolution
    [V,D] = eig(H);
    e = diag(D);
    temp = V'*phi;
    trans = exp(-1i*e*dt);
    phi1 = temp.*trans;
    phi = V*phi1;

    phit2 = abs(phi).^2;
    Sz_order(i) = S_z'*phit2/(2*L);
    Et(i) = phi'*H0*phi;
end

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