% Several numerical examples from the paper Haoze He, Daniel Kressner, 
% Bor Plestenjak, "Randomized methods for computing joint eigenvalues, 
% with applications to multiparameter eigenvalue problems and root finding"
%
% To start an example, set variable example to:
%   - example=1.1 : Left Figure 2 from Example 5.1 (default example)
%   - example=1.2 : Right Figure 2 from Example 5.1 
%   - example=1.3 : Figure 3 from Example 5.1
%   - example=1.4 : Right Figure 4 from Example 5.2 (empirical failure of 1-S RQ)
%   - example=2.1 : Left Figure 5 from Example 5.3
%   - example=2.2 : Right Figure 5 from Example 5.3
%   - example=3   : Example 5.4, set also deltadif
%   - example=4.1 : Left Figure 8 from Example 5.6
%   - example=4.2 : Right Figure 8 from Example 5.6
%
% For example=3, set also the variable deltadif to 1e-2, 1e-4, ..., 1e-16
%
% To change the sample size (default is 10000), set variable Nsamples
%
% 

% Bor Plestenjak 2024
 
if exist("example")~= 1
    example = 1.1 % this is the default example 
end

if exist("deltadif")~=1
    deltadif = 1e-6 % change this parameter (applies only to example=3) to 1e-2, 1e-4, ..., 1e-14, 1e-16
end

if exist("Nsamples")~=1
    Nsamples = 10000 % number of experiments for each noise level
end

fprintf('Running example %3.1f with deltadif = %5.1e and Nsamples=%d  \n',example,deltadif,Nsamples)

%% ======================================  
% Try not to change this part

% Default parameters of the experiments 
rng(1) % to make results repetable
index = [1]; % default cluster A (can be a row of indices)
noise_level = [-Inf -10 -8 -6];
no_backward = 0; % set 1 or 2 to compute eigenvectors of A(mu) in high precision to simulate no backward and roundoff error (for double precision)
noise_for_Amu =  0; %1e-12; % how much noise we add to A(mu) matrix (relative to the norm of A1) in an attempt to regularize a case with Jordan blocks
plot_figures = 1:9; % selection of figures that we want to plot
use_mp = 0; % use multiple precision in Advanpix Multiple Precision Toolbox
uniform_conditioned = 1; % all parts of eigenvector metrix X should have similar condition

if use_mp 
    class_t = 'mp';
else
    class_t = 'double';
end

if (example > 1) && (example < 2)
    % nice example, all eigenvalues are well separated and cond(X) is small
    % Example 5.1 (Figure 2 left) in the manuscript
    da = [1 1 1 2 2 2 3]';  
    db = [1 2 3 1 2 3 3]';
    DA = diag(da);
    DB = diag(db);
    noise_level = [-inf -14 -12 -10];
    kappa = 1e2;
    string_eig = '(1,1)';
    if example == 1.1
        plot_figures = [1];
    elseif example == 1.2
        string_eig = '(1,1)';
        plot_figures = [1];
        no_backward = 2;  % we use MP for eigs and Rayleigh quotient
    elseif example == 1.3
        string_eig = '(1,1)';
        noise_level = [-14 -10];
        plot_figures = [2 4];    
    elseif example == 1.4
        string_eig = '(1,1)';
        noise_level = [-12 -11 -10 -9];
        if exist("Nsamples")~=1
            Nsamples = 1000000;
        end
        plot_figures = [10];    
    end
elseif (example > 2) && (example < 3)
    % eigenvalues are well separated but cond(X) is very large. 
    % eigenvalue A is much less ill-conditioned than eigenvalue B and
    % the eigenvalue error for 2RQ is smaller
    % the residuals are small for both 1RQ and 2RQ
    da = [1 1 1 2 2 2 3]';  
    db = [1 2 3 1 2 3 3]';
    DA = diag(da);
    DB = diag(db);
    noise_level = [-inf -12 -10 -8];
    kappa = 1e4;  % condition number of eigenvalue matrix X (try 1e2, 1e4, 1e6, 1e8)
    uniform_conditioned = 0;
    plot_figures = [1];    
    if example == 2.1
        index = 1; string_eig = '(1,1)';
    elseif example == 2.2
        index = 4; string_eig = '(2,1)';
    end
elseif example == 3
    % cluster of close eigenvalues 
    %da = [1+deltadif 1-deltadif 1 2 2 2 3]';  
    %db = [1-deltadif 1+deltadif 1 1 2 3 3]';
    da = [1 1+deltadif 1-deltadif 2 2 2 3]';  
    db = [1 1-deltadif 1+deltadif 1 2 3 3]';
    DA = diag(da);
    DB = diag(db);
    kappa = 1e4  % condition number of eigenvalue matrix X (try 1e2, 1e4, 1e6, 1e8)
    noise_level = [-14:2:-6];
    index = 1; 
    string_eig = '(1,1)';
    plot_figures = [1];
elseif (example > 4) && (example < 5)
    % Jordan blocks of size 2 and 3, for multiple eigenvalues it makes sense to look
    % clusters of indices
    da = [1 1 1 2 3 4]'; 
    db = [1 1 1 4 3 2]';
    kappa = 1e1;  % condition number of eigenvalue matrix X (try 1e2, 1e4, 1e6, 1e8)
    DA = diag(da)+ diag([1 1 0 1 0],1);
    DB = diag(db)+ diag([1 1 0 1 0],1);    
    noise_level = [-Inf -14 -12 -10];
    plot_figures = [1];     
    if example == 4.1
        index = 1; string_eig = '(1,1)';
    elseif example == 4.2
        index = 4; string_eig = '(2,4)';
    end
else
    error('This example is not implemented')
end

if use_mp
    da = mp(da);
    db = mp(db);
    DA = mp(DA);
    DB = mp(DB);
    kappa = mp(kappa);
    noise = mp(noise);
end

% Set exact eigenvalues with da and db, matrices A and B are constructed
% as A = X1*diag(da)*inv(X1) and B = X1*diag(db)*inv(X1), where X1 is
% a random matrix with prescribed condition number kappa. If
% uniform_conditioned=1, then all parts of X1 should have similar condition
% number, if uniform_conditioned=0, then first two columns are well
% conditioned and the second part is ill-conditioned (depending on kappa).
% Set:
%   - indexA : index of the first eigenvalue than we are interested in
%   - indexB : index of the second eigenvalue than we are interested in
%   - kappa : condition number of eigenvector matrix
%   - da, db: eigenvalues of A and B
%   - noise_level: vector of exponents used for the noise that is added to A(mu)
%   - N: number of experiments for each noise level

%% ======================================  
% Try not to change this part

exact_eig = [da db];
n = length(da);

center_eig = sum(exact_eig(index,:),1)/length(index);
if length(index)==1
    delta = 0;
else
    delta = norm(exact_eig(index(1),:)-center_eig);
    for j = 2:length(index)
        delta = max(delta, norm(exact_eig(index(j),:)-center_eig));
    end
end
delta = delta*sqrt(length(index));

rng(1)
if kappa == 1
    X1 = eye(n,class_t);
else
    if uniform_conditioned == 1
        X1 = condmatX(n,kappa,class_t); 
    else
        % first part with two vectors is well conditioned
        X2 = condmatX(n-2,kappa,class_t); 
        X1 = zeros(n,class_t); X1(1:2,1:2) = eye(2,class_t); X1(3:n,3:n) = X2; % ill conditioned X but small kappa1
        X1 = randn(n,class_t)*X1;
    end
end
A1 = X1*DA*inv(X1);
B1 = X1*DB*inv(X1);
Aset = {A1,B1};
Y1 = inv(X1');
d = length(Aset);
realfamily = 1;
for j = 1:d
    realfamily = and(realfamily,isreal(Aset{j}));
end

gapA = Inf;
for j = 1:n
    if j~=index
        gapA = min(norm(exact_eig(j,:)-exact_eig(index,:)),gapA);
    end
end

materr2SA = [];
materr1SA = [];
materr2X = [];
materr1X = [];
cndhist = [];
cndevA = [];
D2A = [];
S2A = [];
D2Aord = [];
S2Aord = [];
T2Aord = [];
T2A = [];
X2normA = [];
resA = [];
res1A = [];
nX1A = [];
nX2A = [];
nY1A = [];
nY2A = [];
nY1Aord = [];
nY2Aord = [];
nDA = [];
nSA = [];
nSB = [];
mY1A = [];
mY2A = [];
mX1A = [];
mX2A = [];
bound_simple1 = [];
bound_simple2_eps = [];
bound_simple2 = [];
bound_schur1 = [];
bound_schur2 = [];
ratio = [];
ratio5 = [];
if use_mp
    noise = 10.^mp(noise_level);
else
    noise = 10.^noise_level;
end
legend_text = cell(1,length(noise_level));
for j = 1:length(noise_level)
    legend_text{j} = sprintf('e%d',noise_level(j));
end
if isinf(noise_level(1))
    legend_text{1} = 'no noise';
end

min_noise = norm([norm(A1) norm(B1)])*eps(class_t)/2;

% random noise matrices we add to A to simmulate near commutativity
PA = cell(d,1);
for j = 1:d 
    if realfamily
        PA{j} = randn(n,class_t);
    else
        PA{j} = randn(n,class_t) + 1i*randn(n,class_t);
    end
    PA{j} = PA{j}/norm(PA{j},'fro');
end

for ind = 1:length(noise_level)
    fprintf('computing for noise %5.2e \n',noise(ind))
    A = Aset;
    for j = 1:d 
        A{j} = A{j} + noise(ind)/sqrt(d)*PA{j};
    end
    use_noise = max(noise(ind),min_noise);
    [rq2, rq1, res, err, cnd, d2, normsX, SA] = testRayleigh(A,Nsamples,exact_eig,index,1,noise_for_Amu,no_backward,Aset,0);
    [~,ord2] = sort(rownorms(err(:,1)));
    [~,ord1] = sort(rownorms(err(:,2)));
    materr2SA = [materr2SA sort(rownorms(err(:,1)))];
    materr1SA = [materr1SA sort(rownorms(err(:,2)))];
    materr2X = [materr2X sort(rownorms(err(:,3)))];
    materr1X = [materr1X sort(rownorms(err(:,4)))];
    ratio = [ratio length(find(rownorms(err(:,1))<rownorms(err(:,2))))/Nsamples];
    ratio5 = [ratio5 length(find(rownorms(err(:,1))<5*rownorms(err(:,2))))/Nsamples];
    cndhist = [cndhist sort(cnd)];
    boundRQ1 = normsX(:,1).^2.*(delta + use_noise*(1 + sqrt(2)*normsX(:,2).*normsX(:,4).*d2(:,1)));
    boundRQ2_eps = use_noise*normsX(:,1).*normsX(:,3);
    boundRQ2 = boundRQ2_eps.*(1 + use_noise*normsX(:,2).*normsX(:,4).*(sqrt(2)*d2(:,3)));
    boundRQ1Schur = use_noise*(1 + sqrt(2)*SA(:,2).*SA(:,4).*SA(:,5));
    D2A =   [D2A sort(rownorms(d2(:,1)))];
    S2A =   [S2A sort(rownorms(d2(:,2)))];
    D2Aord =   [D2Aord rownorms(d2(ord2,1))];
    S2Aord =   [S2Aord rownorms(d2(ord2,2))];
    T2Aord =   [T2Aord rownorms(d2(ord2,3))];
    T2A =   [T2A sort(rownorms(d2(:,3)))];
    X2normA =   [X2normA sort(normsX(:,2))];
    resA  = [resA  sort(res(:,1))];
    res1A = [res1A sort(res(:,2))];
    nX1A = [nX1A sort(normsX(:,1))]; 
    nX2A = [nX2A sort(normsX(:,2))];
    nY1A = [nY1A sort(normsX(:,3))]; 
    nY2A = [nY2A sort(normsX(:,4))];
    nY1Aord = [nY1Aord normsX(ord2,3)]; 
    nY2Aord = [nY2Aord normsX(ord2,4)]; 
    nDA = [nDA sort(normsX(:,5))];
    nSA = [nSA sort(SA(:,5))];
    mX1A = [mX1A sort(SA(:,1))];
    mX2A = [mX2A sort(SA(:,2))];
    mY1A = [mY1A sort(SA(:,3))];
    mY2A = [mY2A sort(SA(:,4))];
    bound_simple1 = [bound_simple1 boundRQ1(ord1)];
    bound_simple2_eps = [bound_simple2_eps boundRQ2_eps(ord2)];
    bound_simple2 = [bound_simple2 boundRQ2(ord2)];
    bound_schur1 = [bound_schur1 sort(boundRQ1Schur)];
end
  
PlotSettings

% -------------- figure 1
% Error distribution of 1-S and 2-S RQ
if find(plot_figures==1)
    f1 = figure;
    Ax1(1) = axes(f1);
    semilogy(materr2X(:,1))
    hold on
    for ind = 2:length(noise_level)
        semilogy(materr2X(:,ind))
    end
    set(gca,'ColorOrderIndex',1)
    for ind = 1:length(noise_level)
        semilogy(materr1X(:,ind),'--')
    end
    lgd1 = legend(legend_text,'Location','southeast');
    title(strcat('error distribution for \lambda=  ',string_eig))
    Ax1(2) = copyobj(Ax1(1),gcf);
    delete(get(Ax1(2),'children'))
    
    % plot helper data, but invisible
    hold on
    H1 = plot(nan, nan, '-', 'Color', [0 0 0], 'Parent', Ax1(2),'Visible', 'on');
    H2 = plot(nan, nan, '--', 'Color', [0 0 0], 'Parent', Ax1(2),'Visible', 'on');
    hold off
    % make second axes invisible
    set(Ax1(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd2 = legend([H1 H2], '2-sided RQ', '1-sided RQ', 'Location', 'south');
    set(lgd2,'color','none')
    set(lgd2, 'TextColor', 'black');
    
    % to export figure use
    % exportgraphics(f1,'figures\Ex2_B.pdf')
end

% -------------- figure 2
% Error distribution of 2-S vs. the bound
if find(plot_figures==2)
    f2 = figure;
    Ax2(1) = axes(f2);
    semilogy(materr2SA(:,1))
    hold on
    for ind = 2:length(noise_level)
        semilogy(materr2SA(:,ind))
    end
    set(gca,'ColorOrderIndex',1)
    for ind = 1:length(noise_level)
        semilogy(bound_simple2(:,ind),'--')
    end
    lgd21 = legend(legend_text,'Location','southeast');
    title(strcat('2-RQ error and bound distribution for \lambda=  ',string_eig))
    Ax2(2) = copyobj(Ax2(1),gcf);
    delete(get(Ax2(2),'children'))
    
    % plot helper data, but invisible
    hold on
    H1 = plot(nan, nan, '-', 'LineWidth', 2, 'Color', [0 0 0], 'Parent', Ax2(2),'Visible', 'on');
    H2 = plot(nan, nan, '--', 'LineWidth', 2, 'Color', [0 0 0], 'Parent', Ax2(2),'Visible', 'on');
    hold off
    % make second axes invisible
    set(Ax2(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd22 = legend([H1 H2], 'error', 'bound', 'Location', 'south');
    set(lgd22,'color','none')
    set(lgd22, 'TextColor', 'black');
end

% -------------- figure 3
% Error distribution of 2-S vs. the bound (dashed)
if find(plot_figures==3)
    f2b = figure;
    Ax2b(1) = axes(f2b);
    semilogy(materr2SA(:,1))
    hold on
    for ind = 2:length(noise_level)
        semilogy(materr2SA(:,ind))
    end
    set(gca,'ColorOrderIndex',1)
    for ind = 1:length(noise_level)
        semilogy(bound_simple2_eps(:,ind),'.')
    end
    lgd21b = legend(legend_text,'Location','southeast');
    title(strcat('2-RQ error and bound distribution for \lambda=  ',string_eig))
    Ax2b(2) = copyobj(Ax2b(1),gcf);
    delete(get(Ax2b(2),'children'))
    
    % plot helper data, but invisible
    hold on
    H1 = plot(nan, nan, '-', 'LineWidth', 2, 'Color', [0 0 0], 'Parent', Ax2b(2),'Visible', 'on');
    H2 = plot(nan, nan, '--', 'LineWidth', 2, 'Color', [0 0 0], 'Parent', Ax2b(2),'Visible', 'on');
    hold off
    % make second axes invisible
    set(Ax2b(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd22b = legend([H1 H2], 'error', 'bound', 'Location', 'south');
    set(lgd22b,'color','none')
    set(lgd22b, 'TextColor', 'black');
end

% -------------- figure 4
% Error distribution of 2-S vs. the bound (points like for 1-S)
if find(plot_figures==4)
    f3 = figure;
    Ax3(1) = axes(f3);
    semilogy(materr1SA(:,1))
    hold on
    for ind = 2:length(noise_level)
        semilogy(materr1SA(:,ind))
    end
    set(gca,'ColorOrderIndex',1)
    for ind = 1:length(noise_level)
        semilogy(bound_simple1(:,ind),'.','MarkerSize',8)
    end
    lgd31 = legend(legend_text,'Location','southeast');
    title(strcat('1-RQ error and bound distribution for \lambda=  ',string_eig))
    Ax3(2) = copyobj(Ax3(1),gcf);
    delete(get(Ax3(2),'children'))
    
    % plot helper data, but invisible
    hold on
    H1 = plot(nan, nan, '-', 'LineWidth', 2, 'Color', [0 0 0], 'Parent', Ax3(2),'Visible', 'on');
    H2 = plot(nan, nan, '.', 'LineWidth', 2, 'Color', [0 0 0], 'Parent', Ax3(2),'Visible', 'on');
    hold off
    % make second axes invisible
    set(Ax3(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd32 = legend([H1 H2], 'error', 'bound', 'Location', 'south');
    set(lgd32,'color','none')
    set(lgd32, 'TextColor', 'black');
end

% -------------- figure 5
% Eigenvector matrix condition number
if find(plot_figures==5)
    figure
    plot(log10(cndhist(:,1)))
    hold on
    for ind = 2:length(noise_level)
        plot(log10(cndhist(:,ind)))
    end
    legend(legend_text)
    title('Eigenvector matrix cond')
end

% -------------- figure 6
% norm of matrix D_2 from exact eigs and computed
if find(plot_figures==6)
    figure
    plot(log10(D2A))
    hold on
    set(gca,'ColorOrderIndex',1)
    plot(log10(nDA),'--')
    legend(legend_text)
    title('norm of matrix D_2 from exact eigs and computed')
end

% -------------- figure 7
% condition number of X1
if find(plot_figures==7)
    figure
    semilogy(nY1A.*nX1A)
    hold on
    set(gca,'ColorOrderIndex',1)
    legend(legend_text)
    title('kappa(X1)')
end

% -------------- figure 8
% norm of matrix X1
if find(plot_figures==8)
    figure
    semilogy(nX1A)
    hold on
    set(gca,'ColorOrderIndex',1)
    legend(legend_text)
    title('norm of matrix X1')
end

% -------------- figure 9
% norm of residuals
if find(plot_figures==9)
    figure
    plot(log10(resA(:,1)))
    hold on
    for ind = 2:length(noise_level)
        plot(log10(resA(:,ind)))
    end
    set(gca,'ColorOrderIndex',1)
    for ind = 1:length(noise_level)
        plot(log10(res1A(:,ind)),'--')
    end
    legend(legend_text)
    title(sprintf('Residual(A), s=%5.1e, nrm=%5.1e',...
        gapA,norm(A1)))
end

Med = round(Nsamples/2);
median = [noise; materr2SA(Med,:); materr1SA(Med,:); cndhist(Med,:); resA(Med,:); res1A(Med,:); ...
          nX1A(Med,:); nX2A(Med,:); nY1A(Med,:); nY2A(Med,:); D2A(Med,:); bound_simple2(Med,:); bound_simple1(Med,:); bound_simple2_eps(Med,:); S2A(Med,:); T2A(Med,:); ...
          materr2X(Med,:); materr1X(Med,:)]';

normXY = [noise; mX1A(Med,:); mX2A(Med,:); mY1A(Med,:); mY2A(Med,:); nSA(Med,:); bound_schur1(Med,:);]';

fprintf('||A1||=%7.2e, ||B1||=%7.2e, delta=%7.2e, min_noise = %7.2e \n',norm(A1),norm(B1),delta,min_noise)
index

fprintf('Median error and condition numbers for selected eigenvalue (cluster)   \n')
fprintf('eps      | 1S RQ err | 2S RQ err | ratio     | ratio 5   | 1S RQ cls | 2S RQ cls | cond(X)   | 1S RQ res | 2S RQ res | 2S RQ eps | 2S RQ bnd | 1S RQ bnd | \n')
fprintf('---------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n')
for k=1:size(median,1)
    fprintf('%5.1e  |  %7.1e  |  %7.1e  |  %6.4f   |  %6.4f   |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e  \n',...
        median(k,1),median(k,18),median(k,17),ratio(k),ratio5(k),median(k,3),median(k,2),median(k,4),median(k,6),median(k,5),median(k,14),median(k,12),median(k,13))
end

fprintf('Norms of matrices X1, X2, Y1, Y2, D2 \n')
fprintf('eps      | ||X1||    | ||X2||    | ||Y1||    | ||Y2||    | ||D2||    | ||S2||    | ||T2||    |\n')
fprintf('---------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n')
for k=1:size(median,1)
    fprintf('%5.1e  |  %7.1e  |  %7.1e  |  %7.2e |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e \n',...
        median(k,1),median(k,7),median(k,8),median(k,9),median(k,10),median(k,11),median(k,15),median(k,16))
end

fprintf('Norms of Schur matrices X1, X2, Y1, Y2, and D2  \n')
fprintf('eps      | ||X1||(S) | ||X2||(S) | ||Y1||(S) |||Y2||(S)  | ||D2||(A) | 1S RQ bnd |\n')
fprintf('---------|-----------|-----------|-----------|-----------|-----------|-----------|\n')
for k=1:size(normXY,1)
    fprintf('%5.1e  |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e  |  %7.1e  \n',...
        normXY(k,1),normXY(k,2),normXY(k,3),normXY(k,4),normXY(k,5),normXY(k,6),normXY(k,7))
end

% -------------- figure 10
% empirical probability of failure for 1-sided RQ
if find(plot_figures==10)
    M = 5000;
    maxC1 = max(max(materr1X)./noise); % /kappa;
    C1 = exp(linspace(0,log(maxC1),M));
    NL = length(noise);
    pr1 = zeros(M,NL);
    for ind = 1:NL
        for j = 1:M
           pr1(j,ind) = length(find(materr1X(:,ind)>(1+C1(j))*noise(ind)))/Nsamples;
        end
    end
    kap = median(k,8)*median(k,10); % ||X_2||*||Y_2||

    figure
    loglog(C1,pr1(:,1))
    hold on
    for ind = 2:NL
        loglog(C1,pr1(:,ind))
    end
    loglog(C1,kap^2*2*(n-1)./(C1.^2),'--')
    legend_text2 = legend_text;
    legend_text2{length(noise)+1}='upper bound';
    legend(legend_text2)
    ylabel('Probability')
    xlabel('R')
end

