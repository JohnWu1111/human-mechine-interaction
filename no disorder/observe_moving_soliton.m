clear;
% close all;
clc;
format long
tic;

myseed = 14;
rng(myseed)

dt = 0.01;
T = 0:dt:1000000*dt;
nt = length(T);
L = 101;
% L_it = floor(L/2);
L_it = floor(L/2)+1;
K = -1;
mu = zeros(L,1);
% Tij = gen_H_2(1,sqrt(L),L);
Tij = gen_H(1,L);
x = (1:L)';

%%%%%%%%%%%% ED %%%%%%%%%%%
step = 1000;
ni = zeros(L,1);
ni(L_it) = 1;
nit0 = ni./sum(ni);
nit = zeros(L,step);
nit(:,1) = nit0;
E_step = zeros(1,step);

target = 1;
for i = 2:step
    H = Tij + diag(mu) + K*diag(nit(:,i-1));
    [V,D] = eigs(sparse(H),1,'smallestreal');
    nit(:,i) = abs(V).^2;
    E_step(i) = D;
end

wp = nit(:,end);
k = 0.1;
ss = 0;
% phi = sqrt(circshift(wp,ss)).*exp(-1i*k*x) + sqrt(circshift(wp,-ss)).*exp(1i*k*x);
% phi = sqrt(circshift(wp,ss)) + sqrt(circshift(wp,-ss));
% phi = sqrt(circshift(wp,ss)+circshift(wp,-ss));
phi = sqrt(wp).*exp(1i*k*x);

% phi = zeros(L,1);
% ss = 2;
% left = L_it-ss;
% right = L_it+ss;
% phi(left) = 1;
% phi(right) = 1;
% phi(L_it) = 1;

Et = zeros(1,nt);
nit = zeros(L,nt);
phit = zeros(L,nt);
nit0 = abs(phi).^2;
nit(:,1) = abs(phi).^2;
phit(:,1) = phi;
nit_sum = zeros(1,nt);
nit_sum(1) = 2-sum(nit0);
% etat = zeros(nt,1);
% etat(1) = phi(1)/phi(2);

H = Tij + K*diag(nit(:,1));
Et(1) = real(phi'*H*phi);
flux = zeros(L,nt);
dflux = zeros(L,nt);
flux(:,1) = imag(conj(phi).*(phi-phi([2:end,1])));
dflux(:,1) = flux(:,1) - flux([2:end,1],1);
% pos_mean = zeros(nt,1);
% var_x2 = zeros(nt,1);
% pos_mean(1) = wmean((1:L)',abs(phi).^2,1);
% var_x2(1) = sqrt(wmean(((1:L)'-pos_mean(1)).^2,abs(phi).^2,1));


% exact ED %%%%%%%%%%%%%%%%%%%%%%%%

phase = exp(1i*2*pi*x/L);
% Tij = gen_H(exp(1i*pi/2),L);
for i = 2:nt
    H = Tij + K*diag(nit(:,i-1));
%     phi = expm(-1i*H*dt)*phi;
    [V,D] = eig(H);
    e = diag(D);
    trans = V'*phi;
    phi = V*(exp(-1i*e*dt).*trans);
    nit(:,i) = abs(phi).^2;
    phit(:,i) = abs(phi);
    flux(:,i) = imag(conj(phi).*(phi-phi([2:end,1])));
    dflux(:,i) = flux(:,i) - flux([2:end,1],1);
%     etat(i) = phi(1)/phi(2);
%     pos_mean(i) = wmean((1:L)',nit(:,i),1);
%     var_x2(i) = sqrt(wmean(((1:L)'-pos_mean(i)).^2,nit(:,i),1));
    Et(i) = real(phi'*H*phi);
    nit_sum(i) = 2-sum(nit(:,i));
end

filename = strcat('L = ',num2str(L), ', K = ', num2str(K), ', seed = ', num2str(myseed), ', dt = ', num2str(dt),', L_it = ',num2str(L_it));
figure('Name',filename);
set(gcf, 'position', [250 70 1500 900]);

subplot(2,2,1)
plot(1:L,nit(:,end))
xlabel('N')
ylabel('final ni')

subplot(2,2,2)
plot(T,Et)
xlabel('T')
ylabel('energy')

subplot(2,2,3)
meshc(nit)

subplot(2,2,4)
% [~,it] = max(nit(:,end));
% plot(T,nit(it,:))
% xlabel('T')
% ylabel('ni of peak')
% plot(x,flux(:,end))
meshc(flux)
xlabel('T')
ylabel('flux')

% saveas(gcf,strcat('figures\',filename,'.fig'))

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

function Tij = gen_HJ(L,phase)
Tij = zeros(L);
count = 0;
for i = 1:L-1    
    Tij(i,i+1) = Tij(i,i+1)-phase(i);
    Tij(i+1,i) = Tij(i+1,i)-conj(phase(i));    
    count = count +1;    
end
Tij(L,1) = Tij(L,1)-phase(L);
Tij(1,L) = Tij(1,L)-conj(phase(L));
count = count +1;
end

function y = wmean(x,phi,dx)
    y = sum(x.*phi)*dx;
end