clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

J = 1;
dt = 1;
M = 1;
T = 0:dt:10000*dt;
dt = dt/M;
nt = length(T);
L = 32;
% L_it = floor(L/2);
L_it = 9;
Tij = gen_H(1,L);
K1 = -1;
K2 = 1;
K = zeros(L,1);
K(1:floor(L/4)) = K1;
K(floor(L/4)+1) = (K1 + K2)/2;
K(floor(L/4)+2:floor(3*L/4)) = K2;
K(floor(3*L/4)+1) = (K1 + K2)/2;
K(floor(3*L/4)+2:end) = K1;
x = 1:L;

phi = zeros(L,1);
phi(L_it) = 1;
% dn = 0.1;
% phi(1) = sqrt((1+dn)/2);
% phi(2) = sqrt(1-(1+dn)/2);
% phi = rand(L,1);
% phi = phi./sqrt(sum(abs(phi).^2));

Et = zeros(1,nt);
nit = zeros(L,nt);
phit = zeros(L,nt);
nit0 = abs(phi).^2;
nit_it = nit0;
nit(:,1) = nit_it;
phit(:,1) = phi;
% etat = zeros(nt,1);
% etat(1) = phi(1)/phi(2);

H = Tij + K.*diag(nit(:,1));
Et(1) = phi'*H*phi;
pos_mean = zeros(nt,1);
var_x2 = zeros(nt,1);
pos_mean(1) = wmean(x',abs(phi).^2,1);
var_x2(1) = sqrt(wmean((x'-pos_mean(1)).^2,abs(phi).^2,1));


% TSSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

miu= -pi + (0:L-1)*2*pi/L;
miu = miu';
pha2 = exp(-1i*dt*2*J*cos(miu));
miu_f = miu/(2*pi);
for i = 2:nt
    for j = 1:M
        phi1 = exp(-1i*dt*(K.*nit_it)/2).*phi;
        phi1f = nufft(phi1,x,miu_f);
        pha2_phi1f = pha2.*phi1f;
        phi2 = nufft(pha2_phi1f,-miu_f,x);
        phi2 = phi2/L;

        phi = exp(-1i*dt*(K.*nit_it)/2).*phi2;
        
    end
    nit_it = abs(phi).^2;
    nit(:,i) = nit_it;
    phit(:,i) = phi;
    pos_mean(i) = wmean(x',nit_it,1);
    %     var_x2(i) = sqrt(wmean((x'-pos_mean(i)).^2,nit(:,i),1));
    Et(i) = real(phi'*H*phi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = strcat('L = ',num2str(L), ', dt = ', num2str(dt),', M = ',num2str(M));
figure('Name',filename);
set(gcf, 'position', [250 70 1900 900]);

subplot(2,3,1)
plot(1:L,nit(:,end))
xlabel('T')
ylabel('final nit')

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
plot(1:L,K)
xlabel('position')
ylabel('nonlinear factor')

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