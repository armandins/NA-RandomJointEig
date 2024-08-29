% Plots Figure 3 (left) in the paper Haoze He, Daniel Kressner, 
% Bor Plestenjak, "Randomized methods for computing joint eigenvalues, with 
% applications to multiparameter eigenvalue problems and root finding" 

% Bor Plestenjak 2024

da = [1 1 1 2 2 2 3]';  
db = [1 2 3 1 2 3 3]';

rng(1)
N = 1000000;
n = length(da);
p = 2;
class_t = 'double';
dif = [];
z = [];
w = [];
d = zeros(N,1);

for i = 1:n
    dif(i,1) = norm([da(i)-da(1)  db(i)-db(1)]);
end

mu = randn(N,2,class_t) + 1i*randn(N,2,class_t);

for j = 1:N
    % mu = randn(p,1,class_t) + 1i*randn(p,1,class_t);
    mu(j,:) = mu(j,:)/norm(mu(j,:));
    z = mu(j,1)*da + mu(j,2)*db;
    w = abs(z(2:n)-z(1))./dif(2:n);
    d(j) = min(w);
end

C1 = min(d)*5;
C2 = max(d);

PlotSettings

M = 500;
sp = exp(linspace(log(C1),log(C2),M));
pr = zeros(1,M);
for j = 1:M
   pr(1,j) = length(find(d>sp(j)))/N;
end
figure
loglog(1./sp,1-pr)
hold on
loglog(1./sp,(n-1)*sp.^2,'--')
hold off
ylabel('Probability')
xlabel('R')
legend('empirical probability','upper bound')