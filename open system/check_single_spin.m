clear;
% close all;
clc;
format long
tic;

myseed = 14;
% rng(myseed)

dt = 0.01;
M = 1;
T = 0:dt:1000*dt;
nt = length(T);
nt_t = (nt-1)*M+1;
dt = dt/M;
L = 2;
% L_it = floor(L/2);
L_it = 1;
K = -0;
mu_A = 2;
% mu = mu_A*(2*rand(1,L)-1);
mu = [0 0];
gamma = 0.1;
Tij = gen_H(1,L);
num = 1e2;

Et_final = zeros(1,nt);
nit_final = zeros(1,nt);
pos_mean_final = zeros(1,nt);
var_x2_final = zeros(1,nt);

for n = 1:num

    phi = zeros(L,1);
    phi(L_it) = 1;
    % dn = 0.0;
    % phi(1) = sqrt((1+dn)/2);
    % phi(2) = sqrt(1-(1+dn)/2);
    % phi = rand(L,1);
    % phi = phi./sqrt(sum(abs(phi).^2));

    Et = zeros(1,nt);
    nit = zeros(L,nt);
    phit = zeros(L,nt);
    nit0 = abs(phi).^2;
    nit(:,1) = abs(phi).^2;
    nit_now = nit(:,1);
    phit(:,1) = phi;
    etat = zeros(nt,1);
    etat(1) = phi(1)/phi(2);

    rhot = zeros(L^2,nt);
    rho = phi*phi';
    rhot(:,1) = rho(:);

    H00 = Tij + diag(mu);
    H0 = H00 + K*diag(nit(:,1));
    Et(1) = phi'*H0*phi;
    pos_mean = zeros(nt,1);
    var_x2 = zeros(nt,1);
    pos_mean(1) = wmean((1:L)',abs(phi).^2,1);
    var_x2(1) = sqrt(wmean(((1:L)'-pos_mean(1)).^2,abs(phi).^2,1));

    % exact ED %%%%%%%%%%%%%%%%%%%%%%%%

    count = 1;
    for i = 2:nt_t
        envir = gamma*(randn(1))*[1;-1];
        envir = diag(envir);

%         if mod(i-1,M) == 0

            H1 = H00 + K*diag(nit_now) + envir;
            %     phi = expm(-1i*H*dt)*phi;
            [V1,D1] = eig(H1);
            e1 = diag(D1);
            trans = V1'*phi;
            phi1 = V1*(exp(-1i*e1*dt).*trans);
            nit1 = abs(phi1).^2;

            H2 = H00 + K*diag(nit1) + envir;
            [V2,D2] = eig(H2);
            e2 = diag(D2);
            trans = V2'*phi;
            phi2 = V2*(exp(-1i*e2*dt).*trans);
            phi = (phi1+phi2)/2;
        if mod(i-1,M) == 0
            count = count + 1;
            nit(:,count) = abs(phi).^2;
            nit_now = nit(:,count);
            phit(:,count) = abs(phi);
            etat(count) = phi(1)/phi(2);
            pos_mean(count) = wmean((1:L)',nit(:,count),1);
            var_x2(count) = sqrt(wmean(((1:L)'-pos_mean(count)).^2,nit(:,count),1));
            Et(count) = real(phi'*H1*phi);

            rho = phi*phi';
            rhot(:,count) = rho(:);
            H0 = Tij + diag(mu) + K*diag(nit_now);
%         else
%             H = H0 + envir;
%             %     phi = expm(-1i*H*dt)*phi;
%             [V,D] = eig(H);
%             e = diag(D);
%             trans = V'*phi;
%             phi1 = V*(exp(-1i*e*dt).*trans);
%             phi2 = V*(exp(-1i*e*2*dt).*trans);
%             phi = (phi1+phi2)/2;
        end

    end
    Et_final = Et_final + Et;
    nit_final = nit_final + nit;
    pos_mean_final = pos_mean_final + pos_mean;
    var_x2_final = var_x2_final + var_x2;
end
Et_final = Et_final/num;
nit_final = nit_final/num;
pos_mean_final = pos_mean_final/num;
var_x2_final = var_x2_final/num;


filename = strcat('L = ',num2str(L), ', mu_A = ', num2str(mu_A), ', K = ', num2str(K), ', seed = ', num2str(myseed), ', dt = ', num2str(dt),', L_it = ',num2str(L_it));
figure('Name',filename);
set(gcf, 'position', [250 70 1900 900]);

subplot(2,3,1)
plot(T,var_x2_final)
xlabel('T')
ylabel('variance')

subplot(2,3,2)
plot(T,nit_final(L_it,:))
xlabel('T')
ylabel('ni of L_it')

subplot(2,3,3)
plot(T,Et_final)
xlabel('T')
ylabel('energy')

subplot(2,3,4)
meshc(nit_final)

subplot(2,3,5)
plot(1:L,mu)
xlabel('position')
ylabel('random field')

dEt = Et(2:end) - Et(1:end-1);
subplot(2,3,6)
% plot(T(floor(nt/2)+1:end),dEt(floor(nt/2):end))
% xlabel('T')
% ylabel('dE')
plot(T,pos_mean_final)
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