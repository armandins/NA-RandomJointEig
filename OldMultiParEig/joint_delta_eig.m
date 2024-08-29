function [lambda,report] = joint_delta_eig(Delta,opts)

%JOINT_DELTA_EIG  Solves a joined system of generalized eigenvalue problems
%
% [lambda,report] = joint_delta_eig(Delta,opts) returns eigenvalues
% of a joint system of generalized eigenvalue problems
%
% Delta{2} z = lambda(1) Delta{1} z 
% ...
% Delta{k+1} z = lambda(k) Delta{1} z
%
% where Delta{1},...,Delta{k+1} are operator determinant matrices related
% to a k-parameter eigenvalue problem. If Delta{1} is nonsingular, then
% matrices inv(Delta{1})*Delta{2},...,inv(Delta{1})*Delta{k+1} commute
%
% Input:
%   - Delta : cell array of size k+1 of matrices Delta{i}, all matrice
%             have to of the same size
%   - opts : options
%
% Options in opts:
%   - singular (0): set to 1 for a singular problem, i.e., det(Delta0)=0
%   - epscluster (1e-6): relative distance between eigenvalues in a cluster
%   - fast (1): use fast algorithm (can fail for multiple eigenvalues) 
%     or slow algorithm (0) with clustering
%   - rrqr (0): for singular problems only, set to 1 to use rank revealing qr
%   - rng_seed (0) : if different from zero, rng(rng_seed) is used at start
%   - all options of auxiliary functions
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - test_shape (1): shows warning if simultaneous block triangularization
%     is not obtained, set to 0 for a slightly faster evaulation but no test
%
% Output:
%   - lambda : matrix of size m x k, each row is an eigenvalue
%   - report: details of compression method in case of a singular problem 
%       rows are [mode m n r s(1) s(r) s(r+1) choice)], 
%       where mode is 1: CR step 1, 2: CR step 2, 3: RC step 1, 4. RC step 2
%       m,n: size of Delta matrices, r: rank
%       s(1), s(r), s(r+1): singular values
%       choice: how was rank determined (see NUMRANK for details)

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 06.11.2022 : extracted from multipareig 

% Validate number of input parameters
narginchk(1, 2);
if nargin<2, opts=[]; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(Delta{:});
end

if isfield(opts,'singular'),   singular   = opts.singular;   else, singular   = 0;                         end
if isfield(opts,'epscluster'), epscluster = opts.epscluster; else, epscluster = numeric_t('1e-6',class_t); end
if isfield(opts,'fast'),       fast       = opts.fast;       else, fast       = 0;                         end
if isfield(opts,'rng_seed'),   rng_seed   = opts.rng_seed;   else, rng_seed   = 0;                         end
if isfield(opts,'test_shape'), test_shape = opts.test_shape; else, test_shape = 1;                         end
refineeps = numeric_t('eps',class_t);

% Default outputs
lambda = numeric_t([],class_t); 
report = numeric_t([],class_t);

k = length(Delta); % number of parameters + 1

[m,n] = size(Delta{1});
for j = 2:k
   [m1,n1] = size(Delta{j});
   if [m1 n1] ~= [m n]
      error('Matrices must be of same size')
   end
end

% Make sure all inputs are of the same numeric type.
for j = 1:numel(Delta)
    if ~isa(Delta{j}, class_t)
         Delta{j} = numeric_t(Delta{j},class_t);
    end
end

if singular
    [Delta, report] = extract_regular_part_np(Delta, opts);
end

if (size(Delta{1},1) == 0) || (size(Delta{1},1)~=size(Delta{1},2))   % no regular part was found
    return
end

% orthogonal transformation of eigenvalues, we assume that the orthogonal 
% transformation that puts inv(Delta0)*Delta1 into triangular form, does 
% this to other matrices inv(Delta0)*Deltaj as well
alpha = orth(rand(k-1,class_t));
N = size(Delta{1},1);
DeltaH = cell(1,k-1);
for i = 1:k-1
   DeltaH{i} = zeros(N,class_t);
   for j = 1:k-1
      DeltaH{i} = DeltaH{i} + alpha(i,j)*Delta{j+1};
   end
end

out_of_shape = 0;

lambdaH = zeros(N,k-1,class_t);
if fast
    % Factorize Delta{1} only once and use factorization to solve against
    % different righ-hand-sides later on. This boosts the speed slightly.
    [L,U,p] = lu(Delta{1},'vector');
    Gamma1 = U\(L\DeltaH{1}(p,:));
    [Q,T] = schur(Gamma1,'complex');
    lambdaH(:,1) = diag(T);
    for r = 2:k-1
        % Avoid O(n^3) operations in the loop
        Gammar = (U\(L\DeltaH{r}(p,:)))*Q;
        if test_shape
            QG = Q'*Gammar; 
            out_of_shape(r) = norm(tril(QG,-1),'fro')/norm(QG,'fro');
            lambdaH(:,r) = diag(QG);
        else
            for i = 1:N
                lambdaH(i,r) = Q(:,i)'*Gammar(:,i);
            end
        end
    end
 else
    [S1,S0,Q,Z,order,start,csize,tlambda] = clustered_qz(DeltaH{1},Delta{1},epscluster); %#ok<*ASGLU>
    lambdaH(:,1) = tlambda;
    if max(csize)==1
        % there are no clusters
        for r = 2:k-1
            DZ = DeltaH{r}*Z;
            if test_shape
                QDZ = Q*DZ;
                out_of_shape(r) = norm(tril(QDZ,-1),'fro')/norm(QDZ,'fro');
                mu = diag(QDZ);
            else
                mu = zeros(N,1,class_t);
                 for i = 1:N
                     mu(i) = Q(i,:)*DZ(:,i);
                 end
            end
            lambdaH(:,r) = mu./diag(S0);
        end
    else
        % we extract eigenvalues block by block, for each block we take an
        % average eigenvalue as a representative
        for r = 2:k-1
            mu = numeric_t([],class_t); 
            S2 = Q*DeltaH{r}*Z;
            if test_shape
                out_of_shape(r) = out_of_shape_norm(S2,start,csize)/norm(S2,'fro');
            end
            for p = 1:length(start)
                if csize(p)==1
                    mu = [mu; S2(start(p),start(p))/S0(start(p),start(p))];
                else
                    partS0 = S0(start(p):(start(p)+csize(p)-1),start(p):(start(p)+csize(p)-1));
                    partS2 = S2(start(p):(start(p)+csize(p)-1),start(p):(start(p)+csize(p)-1));
                    tmpeig = eig(partS2,partS0);
                    avgmu = sum(tmpeig)/csize(p)*ones(csize(p),1);
                    mu = [mu; avgmu];
                end
            end
            lambdaH(:,r) = mu;
        end
    end
end

if max(out_of_shape)>1e-6
    warning(sprintf('Simultaneous triangularization relative error: %0.2g',max(out_of_shape)))
end

% Inverse unitary transformation to obtain eigenvalues
lambda = zeros(N,k-1,class_t);
for r=1:size(lambdaH,1)
    tmp = alpha\(lambdaH(r,:).');
    lambda(r,:) = tmp.';
end

end % multipareig

% auxiliary function returns Frobenious norm of the part outside the block
% triangular shape
function off = out_of_shape_norm(A,start,csize)
    off = 0;
    for k=1:length(start)
        off = off + norm(A(start(k)+csize(k):end,start(k):start(k)+csize(k)-1),'fro')^2;
    end
    off = sqrt(off);
end