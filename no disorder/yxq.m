clear;
tic;

L = 1000;
dt = 1e-3;
M = 1;
T = 0:M*dt:100;
nt = length(T);

k = -pi/2 + 2*pi/L:2*pi/L:pi/2;
E_k = -2*cos(k');

Vf = 3;
delta = 0.8;
Vi = 10;
w = 1;

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
    for j = 1:L/2
        H = [E_k(j) -2*m0(i-1)*Vi;-2*m0(i-1)*Vi -E_k(j)];
        [V,D] = eig(H);
        phi1(j) = V(1,1);
        phi2(j) = V(2,1);
    end
    m0(i) = 2*phi1'*phi2/L;
end
phi = [phi1;phi2];

m_it = m0(end);
m(1) = m_it;
% m(1) = phi1'*phi2/L;
% S1 = [1:L/2,L/2+1:L,1:L/2,L/2+1:L];
% S2 = [1:L/2,L/2+1:L,L/2+1:L,1:L/2];
% expH = spalloc(L,L,2*L);
% expH = gpuArray(expH);
% S3 = zeros(1,2*L);

% expH =zeros(2);
% temp = zeros(2,1);

G = cell(L/2,1);
for i = 1:L/2
    G{i} = zeros(2);
    G{i}(1,1) = phi1(i)^2;
    G{i}(2,2) = phi2(i)^2;
    G{i}(1,2) = phi1(i)*phi2(i);
    G{i}(2,1) = phi1(i)*phi2(i);
end

count = 2;
t_it = 0;
for i = 2:nt*M
    V_it = Vf+delta*cos(2*pi*w*t_it);
    t_it = t_it + dt;
%         H = sparse(Tij + Vij*m(i-1)*V_it);
%         phi = myrunge(-1i*H,phi,dt);
%         phi1 = phi(1:L);
%         phi2 = phi(L+1:2*L);
%     phi1 = phi1 - 1i*(E_k.*phi1+m(i-1)*V_it*phi2)*dt;
%     phi2 = phi2 - 1i*(-E_k.*phi2+m(i-1)*V_it*phi1)*dt;

%     c11 = - 1i*(E_k.*phi1-2*m(i-1)*V_it*phi2);
%     c12 = - 1i*(-E_k.*phi2-2*m(i-1)*V_it*phi1);
%     c21 = - 1i*(E_k.*phi1-2*m(i-1)*V_it*phi2) - 1i*(E_k.*c11-2*m(i-1)*V_it*c12)*dt/2;
%     c22 = - 1i*(-E_k.*phi2-2*m(i-1)*V_it*phi1) - 1i*(-E_k.*c12-2*m(i-1)*V_it*c11)*dt/2;
%     c31 = - 1i*(E_k.*phi1-2*m(i-1)*V_it*phi2) - 1i*(E_k.*c21-2*m(i-1)*V_it*c22)*dt/2;
%     c32 = - 1i*(-E_k.*phi2-2*m(i-1)*V_it*phi1) - 1i*(-E_k.*c22-2*m(i-1)*V_it*c21)*dt/2;
%     c41 = - 1i*(E_k.*phi1-2*m(i-1)*V_it*phi2) - 1i*(E_k.*c31-2*m(i-1)*V_it*c32)*dt;
%     c42 = - 1i*(-E_k.*phi2-2*m(i-1)*V_it*phi1) - 1i*(-E_k.*c32-2*m(i-1)*V_it*c31)*dt;
%     phi1 = phi1 + dt*(c11+2*c21+2*c31+c41)/6;
%     phi2 = phi2 + dt*(c12+2*c22+2*c32+c42)/6;
    
    b = 2*m_it*V_it;
%     j = 1;
%     while j <= L/2
%         a = E_k(j);
%         fact = sqrt(a^2+b^2);
%         ss = sin(fact*dt)/fact;
%         cc = cos(fact*dt);
%         iEs = 1i*a*ss;
%         ibs = 1i*b*ss;
% %         expH = [cos(fact*dt)-1i*a*sin(fact*dt)/fact -1i*b*sin(fact*dt)/fact;
% %             -1i*b*sin(fact*dt)/fact cos(fact*dt)+1i*a*sin(fact*dt)/fact];
%         expH(1,1) = cc-iEs;
%         expH(2,2) = cc+iEs;
%         expH(1,2) = ibs;
%         expH(2,1) = ibs;
%         temp(1) = phi1(j);
%         temp(2) = phi2(j);
%         temp = expH*temp;
%         phi1(j) = temp(1);
%         phi2(j) = temp(2);
%         j =j +1;
%     end

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

%     m_it = 0;
%     for j = 1:L/2
%         a = E_k(j);
%         fact = sqrt(a^2+b^2);
%         ss = sin(fact*dt)/fact;
%         cc = cos(fact*dt);
%         iEs = 1i*a*ss;
%         ibs = 1i*b*ss;
%         expH(1,1) = cc-iEs;
%         expH(2,2) = cc+iEs;
%         expH(1,2) = ibs;
%         expH(2,1) = ibs;
%         G{j} = G{j}*expH;
%         G{j} = expH\G{j};
%         m_it  = m_it + (G{j}(1,2)+G{j}(2,1))/L;
%     end


%     fact = sqrt(E_k.^2+b^2);
%     ss = sin(fact*dt)./fact;
%     cc = cos(fact*dt);
%     iEs = 1i*E_k.*ss;
%     ibs = 1i*b*ss;
%     S3(1:L/2) = cc-iEs;
%     S3(L/2+1:L) = cc+iEs;
%     S3(L+1:3*L/2) = ibs;
%     S3(3*L/2+1:2*L) = ibs;
%     expH = sparse(S1,S2,S3);
%     phi = expH*phi;
%     phi1 = phi(1:L/2);
%     phi2 = phi(L/2+1:L);

%     phi_sum = sum(abs(phi1).^2+abs(phi2).^2);
%     phi_sum = sqrt(phi_sum*2/L);
%     phi1 = phi1/phi_sum;
%     phi2 = phi2/phi_sum;
    m_it = (phi1'*phi2 + phi2'*phi1)/L;
    if mod(i-1,M) == 0
        m(count) = real(m_it);
%         phi_sum_store(count) = phi_sum -1;
        count = count + 1;
    end
end

toc;

figure
% plot(T(floor(95*nt/100):end),m(floor(95*nt/100):end))
plot(T,m)

function y = myrunge(H,phi,dt)
c1 = H*phi;
c2 = H*(phi+c1.*(dt/2));
c3 = H*(phi+c2.*(dt/2));
c4 = H*(phi+c3.*dt);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end