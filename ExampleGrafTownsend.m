% Plots Figure 7 in the paper Haoze He, Daniel Kressner, Bor Plestenjak, 
% "Randomized methods for computing joint eigenvalues, with applications to 
% multiparameter eigenvalue problems and root finding" 

% This example is related to example from E. Graf and A. Townsend, 
% "Numerical instability of algebraic rootfinding methods", 
% arXiv:2408.02805, 2024.

% Bor Plestenjak 2024

% We take a random 2x2 orthogonal matrix Q and sigma>0 and form polynomials
%
% p1(x,y) = (x-1/3)^2 + sigma*Q(1,1)*(x-1/3) + sigma*Q(1,2)*(y-1/3)
% p2(x,y) = (y-1/3)^2 + sigma*Q(2,1)*(x-1/3) + sigma*Q(2,2)*(y-1/3)
%
% we use biroots to find the roots

if ~isempty(strfind(path, ['OldMultiParEig', pathsep]))
    rmpath OldMultiParEig
end
addpath NewMultiParEig
addpath BiRoots

N = 31;
delta_eksp = linspace(-8,-0.5,N);
runs = 1000;
opts = [];
opts.refine = 0;
opts.twosidedRQ = 1;
method = 0; % this uses the default bipoly_detrep_quintic method for the linearization

mediana = [];

for j = 1:N
    fprintf('Run new %d/%d\n',j,N)
    sigma = 10^delta_eksp(j);
    napaka = [];
    for k=1:runs
        rng(k);
        A = randn(2);
        [Q,~] = qr(A);
        P1 = [1/9-sigma*Q(1,1)/3-sigma*Q(1,2)/3  sigma*Q(1,2)  0; -2/3+sigma*Q(1,1) 0 0; 1 0 0];
        P2 = [1/9-sigma*Q(2,1)/3-sigma*Q(2,2)/3  -2/3+sigma*Q(2,2)  1; sigma*Q(2,1) 0 0; 0 0 0];
        [x,y] = biroots(P1,P2,method,opts);
        razlika = [x y]-[1/3 1/3];
        napaka(k,1) = min(sqrt(diag(razlika*razlika')));
    end
    napaka = sort(napaka);
    mediana(j,1) = napaka(runs/2);
end

% revert to old multipareig
rmpath NewMultiParEig
addpath OldMultiParEig

mediana2 = [];

for j = 1:N
    fprintf('Run old %d/%d\n',j,N)
    sigma = 10^delta_eksp(j);
    napaka = [];
    for k=1:runs
        rng(k);
        A = randn(2);
        [Q,~] = qr(A);
        P1 = [1/9-sigma*Q(1,1)/3-sigma*Q(1,2)/3  sigma*Q(1,2)  0; -2/3+sigma*Q(1,1) 0 0; 1 0 0];
        P2 = [1/9-sigma*Q(2,1)/3-sigma*Q(2,2)/3  -2/3+sigma*Q(2,2)  1; sigma*Q(2,1) 0 0; 0 0 0];
        [x,y] = biroots(P1,P2,method,opts);
        razlika = [x y]-[1/3 1/3];
        napaka(k,1) = min(sqrt(diag(razlika*razlika')));
    end
    napaka = sort(napaka);
    mediana2(j,1) = napaka(runs/2);
end

rmpath OldMultiParEig
rmpath BiRoots

PlotSettings
f1 = figure;
plot(delta_eksp,log10(1./mediana),'*b',delta_eksp,log10(1./mediana2),'*r',delta_eksp,16+2*delta_eksp,'r--',delta_eksp,16+delta_eksp,'b--')
ylabel('Digits of Accuracy')
xlabel('log(\sigma)')
legend('New MultiParEig','Old MultiParEig','Theoretical bound','Stable performance','Location','southeast')

