clear;
% close all;
clc;
format long
tic;

L = 10;
No = 1;
K_all = -1:-1:-6;
nK = length(K_all);

rng(1e5*No)

step_max = 1000;
mu = 2*(2*rand(L,1)-1);
Tij = gen_H(1,L);
H1 = Tij + diag(mu);

x = (1:L)';
stable_store = zeros(L,2);
step_store = zeros(L,1);
result = zeros(L,L);
pos_mean = zeros(L,1);

dt = 0.5;
T = 0:dt:dt*10000;
nt = length(T);

num_stable = zeros(nK,1);
num_substable = zeros(nK,1);

for k = 1:nK
    K = K_all(k);

    for n = 1:L
        target = n;

        phi = zeros(L,1);
        phi(target) = 1;
        nit = zeros(L,1);
        nit(target) = 1;
        step = 0;

        while true
            H = H1 + K*diag(nit);
            %     [V,D] = eigs(sparse(H),3,'smallestreal');
            [V,D] = eig(H);
            VV = V.^2;
            [~,peak_pos] = max(VV.^2);
            cand = find(peak_pos == n);
            if isempty(cand)
                break
            elseif length(cand) == 1
                it = cand;
            else
                VV_cand = VV(:,cand);
                [~,it] = max(VV(n,:));
            end

            nit_new = abs(V(:,it)).^2;
            step = step + 1;

            if step > step_max
                break;
            end

            judge = abs(sum(x.*(nit_new -nit)));
            if judge < 1e-7
                nit_it = nit(:,end);
                pos_mean_it = sum(x.*nit_it);
                [~,peak] = max(nit_it);
                if (n == 1 || abs(pos_mean_it - pos_mean(n-1)) > 0.1) && peak == n

                    stable_store(n,1) = 1;
                    result(n,:) = nit_it';
                    step_store(n) = step;
                    pos_mean(n) = pos_mean_it;

                    nit2 = nit_it.*(1 + 0.2.*(2*rand(L,1)-1));
                    phi = V(:,it);
                    Et0 = phi'*H*phi;
                    for i = 2:nt
                        H = H1 + K*diag(nit2);
                        [V,D] = eig(H);
                        e = diag(D);
                        trans = V'*phi;
                        phi = V*(exp(-1i*e*dt).*trans);
                        nit2 = abs(phi).^2;
                    end
                    Et_end = real(phi'*H*phi);

                    if abs(Et_end - Et0) < 1e-2
                        stable_store(n,2) = 1;
                    end
                end
                break
            end
            nit = nit_new;
        end

    end

    num_stable(k) = sum(stable_store(:,2));
    num_substable(k) = sum(stable_store(:,1));

end

clear i n k

save(strcat('main_L', num2str(L), '_K', num2str(K), '_No', num2str(No), '.mat'), '-v7.3');

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