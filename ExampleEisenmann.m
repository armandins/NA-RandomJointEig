% Plots Figure 7 in the paper Haoze He, Daniel Kressner, Bor Plestenjak, 
% "Randomized methods for computing joint eigenvalues, with applications to 
% multiparameter eigenvalue problems and root finding" 

% This numerical example is related to example in H. Eisenmann, "A Newton 
% method for solving locally definite multiparameter eigenvalue problems 
% by multiindex", arXiv:2408.04194, 2024.

% Bor Plestenjak 2024

% matrices of the three-parameter eigenvalue problem
M = 3;
A = cell(M,M+1); 
nspan = 4:16;  % span of matrix sizes
runs = 10;
good_bound = 1e-6; % if residual is smaller, we tag the eigenvalue as correct
opts = [];
opts.refine = 0;
opts.twosidedRQ = 1;

maxerr = [];
time1 = [];
maxres = zeros(runs,1);
errors = zeros(runs,1);

% First batch - we use new MultiParEig and twosided Raqyleigh quotients
if ~isempty(strfind(path, ['OldMultiParEig', pathsep]))
    rmpath OldMultiParEig
end
addpath NewMultiParEig
for nn = 1:length(nspan)
    n = nspan(nn);
    fprintf("\n =================\n n=%d\n=================\n",n)

    tic
    for step = 1:runs
        rng(step)
        for i = 1:M
            T = randn(n);
            A{i,1} = triu(T) + triu(T,1)';
            for j = 1:M
                P = randn(n);
                [Q,~] = qr(P);
                D = diag(linspace(-1/(2*M),1/(2*M),n));
                A{i,j+1} = Q*D*Q';
            end
            A{i,i+1} = A{i,i+1} + eye(n);
        end
        [lambda,X] = multipareig(A,opts);

        % compute maximum residual
        m = length(lambda);
        res = zeros(m,M);
        for j = 1:m
            for i = 1:M
                S = A{i,1};
                for p = 1:M
                    S = S - lambda(j,p)*A{i,p+1};
                end
                res(j,i) = norm(S*X{j,i});
            end
        end
        maxres(step) = max(max(res));
        errors(step) = length(find(max(res')>good_bound));
        fprintf('Step %d, maximal residual : %8.3e,  errors: %3d\n', step, maxres(step), errors(step))
    end
    fprintf('Maximal overall residual : %8.3e\n', max(maxres))
    fprintf('Missed eigenvalues : %d out of %d in %d cases \n', sum(errors), n^3*runs,length(find(maxres>good_bound)))
    t1 = toc;
    fprintf('Time for one run: %10.4g\n', t1/runs)
    time1 = [time1 t1/runs];
    maxerr = [maxerr maxres];
end

% Second batch - we use new MultiParEig and onesided Raqyleigh quotients
opts.twosidedRQ = 0;
maxerr3 = [];
time3 = [];
maxres = zeros(runs,1);
errors = zeros(runs,1);
for nn = 1:length(nspan)
    n = nspan(nn);
    fprintf("\n =================\n n=%d (one-sided) \n=================\n",n)

    tic
    for step = 1:runs
        rng(step)
        for i = 1:M
            T = randn(n);
            A{i,1} = triu(T) + triu(T,1)';
            for j = 1:M
                P = randn(n);
                [Q,~] = qr(P);
                D = diag(linspace(-1/(2*M),1/(2*M),n));
                A{i,j+1} = Q*D*Q';
            end
            A{i,i+1} = A{i,i+1} + eye(n);
        end
        [lambda,X] = multipareig(A,opts);

        % compute maximum residual
        m = length(lambda);
        res = zeros(m,M);
        for j = 1:m
            for i = 1:M
                S = A{i,1};
                for p = 1:M
                    S = S - lambda(j,p)*A{i,p+1};
                end
                res(j,i) = norm(S*X{j,i});
            end
        end
        maxres(step) = max(max(res));
        errors(step) = length(find(max(res')>good_bound));
        fprintf('Step %d, maximal residual : %8.3e,  errors: %3d\n', step, maxres(step), errors(step))
    end
    fprintf('Maximal overall residual : %8.3e\n', max(maxres))
    fprintf('Missed eigenvalues : %d out of %d in %d cases \n', sum(errors), n^3*runs,length(find(maxres>good_bound)))
    t1 = toc;
    fprintf('Time for one run: %10.4g\n', t1/runs)
    time3 = [time3 t1/runs];
    maxerr3 = [maxerr3 maxres];
end

% Third batch - we use old MultiParEig
rmpath NewMultiParEig
addpath OldMultiParEig
opts = [];
opts.refine = 0;
maxerr2 = [];
maxres = zeros(runs,1);
errors = zeros(runs,1);
time2 = [];
for nn = 1:length(nspan)
    n = nspan(nn);
    fprintf("\n =================\n n=%d old method \n=================\n",n)
        
    tic
    for step = 1:runs
        rng(step)
        for i = 1:M
            T = randn(n);
            A{i,1} = triu(T) + triu(T,1)';
            for j = 1:M
                P = randn(n);
                [Q,~] = qr(P);
                D = diag(linspace(-1/(2*M),1/(2*M),n));
                A{i,j+1} = Q*D*Q';
            end
            A{i,i+1} = A{i,i+1} + eye(n);
        end
        [lambda,X] = multipareig(A,opts);
        
        % compute maximum residual
        m = length(lambda);
        res = zeros(m,M);
        for j = 1:m
            for i = 1:M
                S = A{i,1};
                for p = 1:M
                    S = S - lambda(j,p)*A{i,p+1};
                end
                res(j,i) = norm(S*X{j,i});
            end
        end
        maxres(step) = max(max(res));
        errors(step) = length(find(max(res')>good_bound));
        fprintf('Step %d, maximal residual : %8.3e,  errors: %3d\n', step, maxres(step), errors(step))
    end
    fprintf('Maximal overall residual : %8.3e\n', max(maxres))
    fprintf('Missed eigenvalues : %d out of %d in %d cases \n', sum(errors), n^3*runs,length(find(maxres>good_bound)))
    t1 = toc;
    fprintf('Time for one run: %10.4g\n', t1/runs)
    time2 = [time2 t1/runs];
    maxerr2 = [maxerr2 maxres];
end

PlotSettings

figure
bigX = [];
bigY = [];
for k = 1:length(nspan)
    bigX = [bigX; nspan(k)*ones(runs,1)];
    bigY = [bigY; maxerr(:,k)];
end 
semilogy(bigX-0.1,bigY,'^b')
bigY2 = [];
for k = 1:length(nspan)
    bigY2 = [bigY2; maxerr2(:,k)];
end 
bigY3 = [];
for k = 1:length(nspan)
    bigY3 = [bigY3; maxerr3(:,k)];
end 
hold on
semilogy(bigX+0.1,bigY3,'og')
semilogy(bigX,bigY2,'sr')
axis([nspan(1)-1 nspan(end)+1 1e-16 1e2])
ylabel('error')
xlabel('size of matrices n')
legend('New MultiParEig 2S','New MultiParEig 1S','Old MultiParEig','Location','southeast')
hold off

figure
semilogy(nspan,time1,'^b-')
hold on
semilogy(nspan,time3,'og--')
semilogy(nspan,time2,'sr:')
legend('New MultiParEig 2S','New MultiParEig 1S','Old MultiParEig','Location','southeast')
ylabel('time in seconds')
xlabel('size of matrices n')
hold off

rmpath OldMultiParEig

