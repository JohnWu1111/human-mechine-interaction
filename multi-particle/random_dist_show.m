clear;
% close all;
clc;
format long
tic;

myseed = 4;
rng(myseed)

dt = 1;
T = 0:dt:2000*dt;
nt = length(T);
L = 100;
% L_it = floor(L/2);
K = -6;
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
% mu = [0.5 0];
Tij = gen_H(1,L);
H1 = Tij + diag(mu);
num = 3;

AF = zeros(L,1);
order = zeros(num,nt);

for n = 1:num

    temp = randperm(L);

    nit0 = zeros(L,1);
    for i = 1:L/2
        nit0(temp(i)) = 1;
    end

    G = diag(nit0);
    nit_it = nit0;
    % order(1) = nit(:,1)'*AF/L;
    order(n,1) = nit_it'*((-1).^(nit0+1))/L;

    % H = H1 + K*diag(nit(:,1));
    % Et(1) = cal_energy(G,mu,K);

    % expH %%%%%%%%%%%%%%%%%%%%%%%%

    for i = 2:nt
        H = H1 + K*diag(nit_it);
        %     expH = expm(-1i*H*dt);
        %     G = expH'*G*expH;
        [V,D] = eig(H);
        e = diag(D);
        expH = exp(-1i*e*dt);
        V_trans = V';
        expHV = expH.*V_trans;
        G = V_trans*G*V;
        G = expHV'*G*expHV;
        nit_it = real(diag(G));
        %     Et(i) = cal_energy(G,mu,K);
        %     order(i) = nit(:,i)'*AF/L;
        order(n,i) = nit_it'*((-1).^(nit0+1))/L;
    end
end

% repeating ED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diag_step = 1000;
% nit_diag = zeros(L,diag_step);
% nit_diag(:,1) = nit0;
% for i = 2:diag_step
%     H = H1 + K*diag(nit_diag(:,i-1));
%     [V,D] = eig(H);
%     e = diag(D);
%     nit_diag(:,i) = sum(abs(V(:,1:L/2)).^2,2);
% end
% E_diag = sum(e(1:L/2));

filename = strcat('L = ',num2str(L), ', mu_A = ', num2str(mu_A), ', K = ', num2str(K), ', seed = ', num2str(myseed), ', dt = ', num2str(dt));
figure('Name',filename);
% set(gcf, 'position', [250 70 1500 900]);
% histogram(order_mean,50,'BinLimits',[0 0.5],'Normalization','pdf','DisplayStyle','stairs')
f = plot(T,order);
xlabel('t')
ylabel('I_n')
legend('run1','run2','run3')

% subplot(2,2,1)
% plot(T,Et)
% xlabel('T')
% ylabel('energy')
%
% subplot(2,2,2)
% plot(T,order)
% xlabel('T')
% ylabel('AF order')
%
% subplot(2,2,3)
% plot(1:L,mu)
% xlabel('position')
% ylabel('random field')
%
% subplot(2,2,4)
% plot(1:L,nit(:,end),1:L,nit_diag(:,end))
% xlabel('position')
% ylabel('nit')
% legend('final nit','GS nit')

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

function y = cal_energy(G,mu,K)
y = 0;
L = length(G);
for i = 1:L-1
    y = y - G(i,i+1) - G(i+1,i);
end
y = y - G(L,1) - G(1,L);
y = y + (mu+K)*diag(G);
y = real(y);
end