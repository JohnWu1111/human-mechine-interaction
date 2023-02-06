clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

dt = 0.01;
T = 0:dt:10000*dt;
nt = length(T);
L = 2;
% L_it = floor(L/2);
L_it = 1;
K = -0;
mu_A = 2;
% mu = mu_A*(2*rand(1,L)-1);
mu = [0 0];
Tij = gen_H(1,L);
gamma = 0.1;

phi = zeros(L,1);
phi(L_it) = 1;
% dn = 1;
% phi(1) = sqrt((1+dn)/2);
% phi(2) = sqrt(1-(1+dn)/2);
% phi = rand(L,1);
% phi = phi./sqrt(sum(abs(phi).^2));

rho = phi*phi';
rho = rho';
rho = rho(:);

Et = zeros(1,nt);
nit = zeros(L,nt);
rhot = zeros(L^2,nt);
nit0 = abs(phi).^2;
nit(:,1) = abs(phi).^2;
rhot(:,1) = rho;
etat = zeros(nt,1);
etat(1) = phi(1)/phi(2);

H = Tij + diag(mu) + K*diag(nit(:,1));
Et(1) = phi'*H*phi;
pos_mean = zeros(nt,1);
var_x2 = zeros(nt,1);
pos_mean(1) = wmean((1:L)',abs(phi).^2,1);
var_x2(1) = sqrt(wmean(((1:L)'-pos_mean(1)).^2,abs(phi).^2,1));

S = zeros(1,nt);
% rho_M = col2mat(rho);
S(1) = 0;


% exact ED %%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:nt
    H = Tij + diag(mu) + K*diag(nit(:,i-1));
%     phi = expm(-1i*H*dt)*phi;
    Liou = gen_liou(H,gamma);
    [V,D] = eig(Liou);
    e = diag(D);
    trans = V\rho;
%     rho = V*(exp(e*dt).*trans);
    rho = expm(Liou*dt)*rho;
    temp = (1:L).^2;
    nit(:,i) = real(rho(temp'));
    rhot(:,i) = rho;
    etat(i) = phi(1)/phi(2);
    pos_mean(i) = wmean((1:L)',nit(:,i),1);
    var_x2(i) = sqrt(wmean(((1:L)'-pos_mean(i)).^2,nit(:,i),1));
    H = H';
    Et(i) = real(sum(H(:).*rho));
    rho_M = col2mat(rho);
    S(i) = -trace(rho_M*logm(rho_M));
end

% runge-kuta %%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 2:nt
%     H = Tij + diag(mu) + K*diag(nit(:,i-1));
% %     phi = myrunge(-1i*H,phi,dt);
%     phi = phi-1i*H*dt*phi;
%     nit(:,i) = abs(phi).^2;
%     pos_mean(i) = wmean((1:L)',nit(:,i),1);
%     var_x2(i) = sqrt(wmean(((1:L)'-pos_mean(i)).^2,nit(:,i),1));
%     Et(i) = real(phi'*H*phi);
% end

filename = strcat('L = ',num2str(L), ', mu_A = ', num2str(mu_A), ', K = ', num2str(K), ', seed = ', num2str(myseed), ', dt = ', num2str(dt),', L_it = ',num2str(L_it));
figure('Name',filename);
set(gcf, 'position', [250 70 1900 900]);

subplot(2,3,1)
plot(T,var_x2)
xlabel('T')
ylabel('variance')

subplot(2,3,2)
plot(T,nit(L_it,:))
xlabel('T')
ylabel('ni of L_it')

subplot(2,3,3)
plot(T,Et)
xlabel('T')
ylabel('energy')

subplot(2,3,4)
meshc(nit)

subplot(2,3,5)
plot(1:L,mu)
xlabel('position')
ylabel('random field')

dEt = Et(2:end) - Et(1:end-1);
subplot(2,3,6)
% plot(T(floor(nt/2)+1:end),dEt(floor(nt/2):end))
% xlabel('T')
% ylabel('dE')
plot(T,pos_mean)
xlabel('T')
ylabel('pos_mean')

% saveas(gcf,strcat('figures\',filename,'.fig'))

mean_dE = mean(dEt(floor(nt*0.9):end));

toc;

function Tij = gen_H(s,L)
Tij = zeros(L);
count = 0;
for i = 1:L-1    
    Tij(i,i+1) = Tij(i,i+1)-s;
    Tij(i+1,i) = Tij(i+1,i)-conj(s);    
    count = count +1;    
end
% Tij(L,1) = Tij(L,1)-s;
% Tij(1,L) = Tij(1,L)-conj(s);
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

function y = gen_liou(H,gamma)
len = length(H);
y1 = zeros(len^2);

for a = 1:len % # of row in left
    for b = 1:len % # of col in right
        for j = 1:len % # of multiplier
            Liou_row = (a-1)*len+b;
            Liou_col1 = (j-1)*len+b;
            Liou_col2 = (a-1)*len+j;
            y1(Liou_row,Liou_col1) = y1(Liou_row,Liou_col1) + H(a,j);
            y1(Liou_row,Liou_col2) = y1(Liou_row,Liou_col2) - H(j,b);
        end
    end
end

y2 = zeros(len^2);
for i = 1:len
    for j = i+1:len
        D_row1 = (i-1)*len+j;
        D_row2 = (j-1)*len+i;
        y2(D_row1,D_row1) = 1;
        y2(D_row2,D_row2) = 1;
    end
end
y = -1i*y1 - y2*gamma;
end

function y = col2mat(A)
    len = length(A);
    L = sqrt(len);
    if mod(L,1) ~= 0
        error("error!")
    else
        y = zeros(L);
        for i = 1:L
            y(:,i) = A(1+(i-1)*L:i*L);
        end
    end
end