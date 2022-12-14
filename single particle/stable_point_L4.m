clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed*1e5)

dn = 0.2;
dt = 0.2;
T = 0:dt:10000*dt;
nt = length(T);
L = 4;
% L_it = floor(L/2);
L_it = 1;
K = -1;
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
% mu = [0 0.5 0.7 1];
Tij = gen_H(1,L);
sign_s = [1 1 1 1 1 1 1 1;1 1 1 1 -1 -1 -1 -1;1 1 -1 -1 1 1 -1 -1;1 -1 1 -1 1 -1 1 -1];

ln = 1/dn+1;
ln_total = ln*(ln+1)*(ln+2)/6;
result = cell(ln_total,8);
result_t = cell(ln_total,8);
mean_pos = zeros(ln_total,8);
std_pos = zeros(ln_total,8);
n0_distri = zeros(ln_total,L);

it = 1;
for m = 1:ln
    for n= 1:ln-m+1
        for p = 1:ln-m-n+2
            n1 = (m-1)*dn;
            n2 = (n-1)*dn;
            n3 = (p-1)*dn;
            n4 = 1-n1-n2-n3;

            nn = [n1,n2,n3,n4]';
            n0_distri(it,:) = nn';
            for k = 1:1
                phi0 = sign_s(:,k).*sqrt(nn);
                phi = phi0;
                nit = zeros(L,nt);
                nit0 = abs(phi).^2;
                nit(:,1) = abs(phi).^2;

                for i = 2:nt
                    H = Tij + diag(mu) + K*diag(nit(:,i-1));
                    %     phi = expm(-1i*H*dt)*phi;
                    [V,D] = eig(H);
                    e = diag(D);
                    trans = V'*phi;
                    phi = V*(exp(-1i*e*dt).*trans);
                    nit(:,i) = abs(phi).^2;
                end
                result{it,k} = nit(:,end);
                result_t{it,k} = nit;
                mean_pos(it,k) = wmean(1:L,nit(:,end)',1);
                std_pos(it,k) = sqrt(wmean(((1:L)-mean_pos(it,k)).^2,nit(:,end)',1));
            end
            it = it + 1;
        end
    end
end

toc;

function Tij = gen_H(s,L)
Tij = zeros(L);
count = 0;
for i = 1:L-1
    Tij(i,i+1) = Tij(i,i+1)-s;
    Tij(i+1,i) = Tij(i+1,i)-conj(s);
    count = count +1;
end
Tij(L,1) = Tij(L,1)-s;
Tij(1,L) = Tij(1,L)-conj(s);
count = count +1;
end

function Tij = gen_H_2(s,L,len)
Tij = zeros(len);
count = 0;
for i = 1:L-1
    for j = 1:L-1
        pos = (i-1)*L+j;
        Tij(pos,pos+1) = Tij(pos,pos+1)-s;
        Tij(pos+1,pos) = Tij(pos+1,pos)-conj(s);
        Tij(pos,pos+L) = Tij(pos,pos+L)-1;
        Tij(pos+L,pos) = Tij(pos+L,pos)-1;
        count = count +1;
    end
    pos = i*L;
    Tij(pos,pos-L+1) = Tij(pos,pos-L+1)-s;
    Tij(pos-L+1,pos) = Tij(pos-L+1,pos)-conj(s);
    Tij(pos,pos+L) = Tij(pos,pos+L)-1;
    Tij(pos+L,pos) = Tij(pos+L,pos)-1;
    count = count +1;
end
for j = 1:L-1
    pos = (L-1)*L+j;
    Tij(pos,pos+1) = Tij(pos,pos+1)-s;
    Tij(pos+1,pos) = Tij(pos+1,pos)-conj(s);
    Tij(pos,pos+L-len) = Tij(pos,pos+L-len)-1;
    Tij(pos+L-len,pos) = Tij(pos+L-len,pos)-1;
    count = count +1;
end
Tij(len,len-L+1) = Tij(len,len-L+1)-s;
Tij(len-L+1,len) = Tij(len-L+1,len)-conj(s);
Tij(len,L) = Tij(len,L)-1;
Tij(L,len) = Tij(L,len)-1;
count = count +1;
end

function y = myrunge(H,phi,dt)
c1 = H*phi;
c2 = H*(phi+c1.*(dt/2));
c3 = H*(phi+c2.*(dt/2));
c4 = H*(phi+c3.*dt);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = wmean(x,phi,dx)
y = sum(x.*phi)*dx;
end