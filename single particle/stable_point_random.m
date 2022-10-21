clear;
% close all;
clc;
format long
tic;

myseed = 3;
rng(myseed)

dt = 0.5;
L = 20;
K = -1;
step_max = 1e6;
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
trial_num = 10;
Tij = gen_H(1,L);
H1 = Tij + diag(mu);

result = zeros(trial_num,L);
result_t = cell(trial_num,1);
osc_store = zeros(trial_num,1);
step_store = zeros(trial_num,1);
nit_mean_final = zeros(trial_num,1);
x = (1:L)';

it = 1;
for j = 1:trial_num
    phi_0 = rand(L,1);
%     phi_0 = rand(L,1);
    phi_0 = phi_0./sqrt(sum(abs(phi_0).^2));
    nit_0 = abs(phi_0).^2;

    phi = phi_0;
    nit = nit_0;

    step = 1;
    while true
        H = H1 + K*diag(nit);
        [V,D] = eig(H);
        e = diag(D);
        trans = V'*phi;
        fact = exp(-1i*e*dt);
        temp = fact.*trans;
        phi = V*temp;
        nit_new = abs(phi).^2;
        step = step + 1;

        if step > step_max
            osc_store(j,1) = 1;
            break;
        end

        judge = abs(sum(x.*(nit_new -nit)));
        if judge < 1e-8
            break
        end
        nit = nit_new;
    end

    result_t{j} = nit(:,end);
    result(j,:) = nit(:,end)';

    nit_mean_final(j) = wmean(x,nit(:,end),1);
    %         nit_std_final(j,n) = sqrt(wmean(((1:L)-nit_mean_final(j,n)).^2,nit_final',1));
    step_store(j) = step;
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