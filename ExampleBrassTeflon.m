% Plots Figure 10 in the paper Haoze He, Daniel Kressner, Bor Plestenjak, 
% "Randomized methods for computing joint eigenvalues, with applications to 
% multiparameter eigenvalue problems and root finding" 

% We use the same data and solvers as in Section 6.2 in
% H. Gravenkamp, B. Plestenjak, D. A. Kiefer, and J. Elias: Computation 
% of leaky waves in layered structures coupled to unbounded media by 
% exploiting multiparameter eigenvalue problems. arXiv:2404.15277, 2024.

% Bor Plestenjak, Hauke Gravenkamp, Daniel Kiefer 2024

clear
close all
clc

%% load matrices 
load example_brassTeflon_Data.mat
R1 = R{1};
R2 = R{2};

w = 2*pi*linspace(0.01, 3, 150).'; % frequency

% options for the new solver for multiparameter eigenvalue problems
opts1 = [];
opts1.refine = 0;
opts1.rand_orth = 1;
opts1.twosidedRQ = 1;
opts1.solver = 'eig';

% options for the old solver for multiparameter eigenvalue problems
opts2 = [];
opts2.refine = 0;

rng(1) % to make results repeatable 

n = size(E0,1);
NN = nan(length(w), 2^3*n);

%% solution 1 - using the new MultiParEig package
if ~isempty(strfind(path, ['OldMultiParEig', pathsep]))
    rmpath OldMultiParEig
end
addpath NewMultiParEig

k1 = NN + 1i*NN;
kyL1 = k1;
kyS1 = k1;
for i = 1:length(w)
    kappa = w(i)./c;
    kappal = kappa(1);
    kappat = kappa(2);
    [lambda1, tmp_lambda1] = eig_LeakySolid(E0,E1,-E2,M,R1,R2,kappal,kappat,w(i),opts1);
    k1(i,1:numel(lambda1))  = lambda1;
    kyL1(i,1:numel(lambda1))  = tmp_lambda1(:,2);
    kyS1(i,1:numel(lambda1))  = tmp_lambda1(:,3);
end

%% solution 2 - using the old MultiParEig package

rmpath NewMultiParEig
addpath OldMultiParEig

k2 = NN + 1i*NN;
kyL2 = k2;
kyS2 = k2;
for i = 1:length(w)
    kappa = w(i)./c;
    kappal = kappa(1);
    kappat = kappa(2);
    [lambda2, tmp_lambda2] = eig_LeakySolid(E0,E1,-E2,M,R1,R2,kappal,kappat,w(i),opts2);
    k2(i,1:numel(lambda2))  = lambda2;
    kyL2(i,1:numel(lambda2))  = tmp_lambda2(:,2);
    kyS2(i,1:numel(lambda2))  = tmp_lambda2(:,3);
end

fh1 = w/2/pi.*ones(size(k1));
fh2 = w/2/pi.*ones(size(k2));

%% filter modes
indRemove1 = (real(kyL1)>-1e-2)|(real(kyS1)>-1e-2);
k1(indRemove1) = nan + 1i*nan;
indRemove2 = (real(kyL2)>-1e-2)|(real(kyS2)>-1e-2);
k2(indRemove2) = nan + 1i*nan;

%% plot wave numbers:
PlotSettings
figure 
hold on
plot(fh1(:), real(k1(:)), 'k.', 'MarkerSize',18,'DisplayName','new solver')
plot(fh2(:), real(k2(:)), 'ro', 'MarkerSize',6,'DisplayName','old solver');
ylabel('wavenumber in rad/mm')
xlabel('frequency in Mhz')
legend('new solver','old solver','Location','southwest')
ylim([0.8, 2.9]) 
xlim([0.5, 3])

rmpath OldMultiParEig