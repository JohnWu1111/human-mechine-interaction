clear;
tic;

L = 2000;
dt = 1;
sigma = 0.01;
M = 1;
T = 0:M*dt:1000;
nt = length(T);

k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
E_k = -2*cos(k');
% E_k = -2*rand(length(k),1);

Vf = 3;
Vi = 1;

phi = ones(L,1);
phi = phi*sqrt(L/(2*sum(abs(phi).^2)));
m = zeros(nt,1);
phi_sum_store = zeros(nt,1);
phi1 = phi(1:L/2);
phi2 = phi(L/2+1:L);
step = 100;
m0 = zeros(step,1);
m0(1) = 2*phi1'*phi2/L;

for i = 2:step
%     E0 = 0;
    for j = 1:L/2
        H = [E_k(j) -2*m0(i-1)*Vi;-2*m0(i-1)*Vi -E_k(j)];
        [V,D] = eig(H);
        phi1(j) = V(1,1);
        phi2(j) = V(2,1);
%         E0 = E0 + V(:,1)'*H*V(:,1);
    end
    m0(i) = 2*phi1'*phi2/L;
    E0 = -2*m0(i).^2*Vi*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
end
phi = [phi1;phi2];

m_iGS = m0(end);
m(1) = m_iGS;
m_it = m_iGS;

Et = zeros(nt,1);
Et(1) = E0;

count = 2;
t_it = 0;
T = zeros(1,nt);
T(1) = 0;
for i = 2:nt*M
    t_it = t_it + dt*rand;
    T(i) = t_it;
    
    b = 2*m_it*Vf;

    fact = sqrt(E_k.^2+b^2);
    ft = fact*dt;
    ss = sin(ft);
    ss = ss./fact;
    cc = cos(ft);
    Es = E_k.*ss;
    bs = b*ss;
    phi1n = (cc-1i*Es).*phi1 +1i*bs.*phi2;
    phi2 = (cc+1i*Es).*phi2 +1i*bs.*phi1;
    phi1 = phi1n;
    m_it = (phi1'*phi2 + phi2'*phi1)/L;
    
    if mod(i-1,M) == 0
        m(count) = real(m_it);
        Et(count) = Et(count) -2*m(count).^2*Vf*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
%         phi_sum_store(count) = phi_sum -1;
        count = count + 1;
    end
    
end

phi = ones(L,1);
phi = phi*sqrt(L/(2*sum(abs(phi).^2)));
phi1 = phi(1:L/2);
phi2 = phi(L/2+1:L);
mf = zeros(step,1);
mf(1) = 2*phi1'*phi2/L;
for i = 2:step
%     E0 = 0;
    for j = 1:L/2
        H = [E_k(j) -2*mf(i-1)*Vf;-2*mf(i-1)*Vf -E_k(j)];
        [V,D] = eig(H);
        phi1(j) = V(1,1);
        phi2(j) = V(2,1);
%         E0 = E0 + V(:,1)'*H*V(:,1);
    end
    mf(i) = 2*phi1'*phi2/L;
    Ef = -2*mf(i).^2*Vf*L + sum(E_k.*(abs(phi1).^2 - abs(phi2).^2));
end
mf_end = mf(end);

toc;

filename = strcat('L = ',num2str(L), ', V_f = ', num2str(Vf), ', V0 = ', num2str(Vi), ', sigma = ', num2str(sigma));
figure('Name',filename);
set(gcf, 'position', [250 70 1500 900]);

subplot(1,2,1)
plot(T,m)
xlabel('t')
ylabel('m')

subplot(1,2,2)
plot(T,Et)
xlabel('T')
ylabel('energy')
