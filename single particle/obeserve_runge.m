clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

dt = 1;
T = 0:dt:10000*dt;
nt = length(T);
L = 2;
% L_it = floor(L/2);
% L_it = floor(L/2)+1;
L_it = 1;
K = -4;
mu_A = 0;
mu = mu_A*(2*rand(1,L)-1);
% mu = zeros(L,1);
% Tij = gen_H_2(1,sqrt(L),L);
Tij = gen_H(1,L);

phi = zeros(L,1);
% ss = 8;
% phi(L_it-ss) = 1;
% phi(L_it+ss) = 1;
phi(L_it) = 1;

% dn = 0.1;
% phi(1) = sqrt((1+dn)/2);
% phi(2) = sqrt(1-(1+dn)/2);
% phi = rand(L,1);
% phi = rand(L,1) + rand(L,1)*1i;
phi = phi./sqrt(sum(abs(phi).^2));

Et = zeros(1,nt);
nit = zeros(L,nt);
phit = zeros(L,nt);
nit0 = abs(phi).^2;
nit(:,1) = abs(phi).^2;
phit(:,1) = phi;
% etat = zeros(nt,1);
% etat(1) = phi(1)/phi(2);

H = Tij + diag(mu) + K*diag(nit(:,1));
Et(1) = real(phi'*H*phi);
pos_mean = zeros(nt,1);
var_x2 = zeros(nt,1);
pos_mean(1) = wmean((1:L)',abs(phi).^2,1);
var_x2(1) = sqrt(wmean(((1:L)'-pos_mean(1)).^2,abs(phi).^2,1));


% exact ED %%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:nt
    H = Tij + diag(mu) + K*diag(nit(:,i-1));
%     phi = expm(-1i*H*dt)*phi;
    [V,D] = eig(H);
    e = diag(D);
    trans = V'*phi;
    phi = V*(exp(-1i*e*dt).*trans);
    nit(:,i) = abs(phi).^2;
    phit(:,i) = abs(phi);
%     etat(i) = phi(1)/phi(2);
    pos_mean(i) = wmean((1:L)',nit(:,i),1);
    var_x2(i) = sqrt(wmean(((1:L)'-pos_mean(i)).^2,nit(:,i),1));
    Et(i) = real(phi'*H*phi);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TSSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% miu= -pi + (0:L-1)*2*pi/L;
% miu = miu';
% pha2 = exp(-1i*dt*2*cos(miu));
% for i = 2:nt
%     phi1 = exp(-1i*dt*(mu'+K*abs(phi).^2)/2).*phi;
%     phi1f = nufft(phi1,1:L,miu/(2*pi));
%     phi2 = nufft(pha2.*phi1f,-miu/(2*pi),1:L)/L;    
%     
%     phi = exp(-1i*dt*(mu'+K*abs(phi).^2)/2).*phi2;
%     nit(:,i) = abs(phi).^2;
%     pos_mean(i) = wmean((1:L)',nit(:,i),1);
%     var_x2(i) = sqrt(wmean(((1:L)'-pos_mean(i)).^2,nit(:,i),1));
%     Et(i) = real(phi'*H*phi);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = strcat('L = ',num2str(L), ', mu_A = ', num2str(mu_A), ', K = ', num2str(K), ', seed = ', num2str(myseed), ', dt = ', num2str(dt),', L_it = ',num2str(L_it));
figure('Name',filename);
set(gcf, 'position', [100 70 1900 900]);

subplot(2,3,1)
plot(T,var_x2)
xlabel('T')
ylabel('variance')

subplot(2,3,2)
plot(1:L,nit(:,end))
xlabel('N')
ylabel('final ni')

subplot(2,3,3)
plot(T,Et)
xlabel('T')
ylabel('energy')

subplot(2,3,4)
meshc(nit)

subplot(2,3,5)
[~,it] = max(nit(:,end));
plot(T,nit(it,:))
xlabel('T')
ylabel('ni of peak')

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
% Tij(len,len-L+1) = Tij(len,len-L+1)-s;
% Tij(len-L+1,len) = Tij(len-L+1,len)-conj(s);
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