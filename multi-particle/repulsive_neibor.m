clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed)

dt = 1;
T = 0:dt:1000*dt;
nt = length(T);
L = 8;
len = 2^L;
% L_it = floor(L/2);
K = -8;
mu_A = 2;
mu = mu_A*(2*rand(1,L)-1);
lambda = -20;
% mu = [0.5 0];

% generate H1 matrix
len_s = factorial(L)/(factorial(L/2)^2);
H1 = zeros(len_s,len_s);
phi_store = zeros(L,len_s);
pos_ori = zeros(len,1);
count = 1;
for i = 1:len
    temp = tentotwo(i,L);
    temp_num = sum(temp);
    if temp_num == L/2
        phi_store(:,count) = temp;
        pos_ori(i) = count;
        count = count +1;
    end  
end

for k = 1:len_s
    phi_k = phi_store(:,k);
    for i = 1:L-1
        j = i+1;
        if phi_k(i) == 1
            H1(k,k) = H1(k,k) + mu(i);
            if phi_k(j) == 0
                phi_b = phi_k;
                phi_b(j) = phi_k(i);
                phi_b(i) = phi_k(j);
                b = twototen(phi_b);
                b = pos_ori(b);
                H1(k,b) = H1(k,b) - 1;
                H1(b,k) = H1(b,k) - 1;
            else
                H1(k,k) = H1(k,k) + lambda;
            end
        end
    end
    i = L;
    j = 1;
    if phi_k(i) == 1
        H1(k,k) = H1(k,k) + mu(i);
        if phi_k(j) == 0
            phi_b = phi_k;
            phi_b(j) = phi_k(i);
            phi_b(i) = phi_k(j);
            b = twototen(phi_b);
            b = pos_ori(b);
            H1(k,b) = H1(k,b) + (-1)^(L/2);
            H1(b,k) = H1(b,k) + (-1)^(L/2);
        else
            H1(k,k) = H1(k,k) + lambda;
        end
    end
end

AF = zeros(L,1);

Et = zeros(1,nt);
order = zeros(1,nt);

nit0 = zeros(L,1);
for i = 1:L
    nit0(i) = ((-1)^(i)+1)/2;
    AF(i) = (-1)^i;
end
% for i = L/2+1:L
%     nit0(i) = 1;
% end

nit = zeros(L,nt);
nit(:,1) = nit0;
order(1) = nit(:,1)'*AF/L;
phi = zeros(len_s,1);
phi(pos_ori(twototen(nit0))) = 1;

% generate H2 matrix
H2 = zeros(len_s,len_s);
for k = 1:len_s
    phi_k = phi_store(:,k);
    for i = 1:L
        if phi_k(i) == 1
            H2(k,k) = H2(k,k) + K*nit(i,1);
        end
    end
end

H = H1 + H2;
Et(1) = real(phi'*H*phi);


% exact ED %%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:nt
    % generate H2 matrix
    H2 = zeros(len_s,len_s);
    for k = 1:len_s
        phi_k = phi_store(:,k);
        for j = 1:L
            if phi_k(j) == 1
                H2(k,k) = H2(k,k) + K*nit(j,i-1);
            end
        end
    end
    H = H1 + H2;
    [V,D] = eig(H);
    e = diag(D);
    trans = V'*phi;
    phi = V*(exp(-1i*e*dt).*trans);
    phi2 = abs(phi).^2;
    nit(:,i) = sum(phi2.*phi_store')';
    Et(i) = real(phi'*H*phi);
    order(i) = nit(:,i)'*AF/L;
end

nit_GS = sum(abs(V(:,1)).^2.*phi_store')';
order_GS = nit_GS'*AF/L;

% repeating ED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diag_step = 1000;
nit_diag = zeros(L,diag_step);
nit_diag(:,1) = nit0;
for i = 2:diag_step
    H2 = zeros(len_s,len_s);
    for k = 1:len_s
        phi_k = phi_store(:,k);
        for j = 1:L
            if phi_k(j) == 1
                H2(k,k) = H2(k,k) + K*nit_diag(j,i-1);
            end
        end
    end
    H = H1 + H2;
    [V,D] = eigs(sparse(H),1,'smallestreal');
    nit_diag(:,i) = sum(abs(V).^2.*phi_store')';
end
E_diag = D;


filename = strcat('L = ',num2str(L), ', mu_A = ', num2str(mu_A), ', K = ', num2str(K), ', seed = ', num2str(myseed), ', dt = ', num2str(dt), ', lambda = ', num2str(lambda));
figure('Name',filename);
set(gcf, 'position', [250 70 1500 900]);

subplot(2,2,1)
plot(T,Et)
xlabel('T')
ylabel('energy')

subplot(2,2,2)
plot(T,order)
xlabel('T')
ylabel('AF order')

subplot(2,2,3)
plot(1:L,mu)
xlabel('position')
ylabel('random field')

subplot(2,2,4)
plot(1:L,nit(:,end),1:L,nit_diag(:,end))
xlabel('position')
ylabel('nit')
legend('final nit','GS nit')

toc;

function y = wmean(x,phi,dx)
    y = sum(x.*phi)*dx;
end

function y = cal_energy(G,mu,K)
y = 0;
L = length(G);
for i = 1:L-1
    y = y + G(i,i+1) + G(i+1,i);
end
y = y + G(L,1) + G(1,L);
y = y + (mu+K)*diag(G);
y = real(y);
end

function y = twototen(phi)
    len = length(phi);
    y = 0;
    for i = 1:len
        y = 2*y+phi(i);
    end
    y = y+1;
end

function y = tentotwo(n,L)
    y = zeros(L,1);
    m = n-1;
    for i = L:-1:1
        y(i) = mod(m,2);
        m = (m - y(i))/2;
    end
end