function [lambda,mu,XR,YR,XL,YL,report] = twopareig(A1,B1,C1,A2,B2,C2,opts)

%TWOPAREIG   Solve a two-parameter eigenvalue problem
%
% [lambda,mu,XR,YR,XL,YL] = TWOPAREIG(A1,B1,C1,A2,B2,C2,opts) returns
% eigenvalues and eigenvectors of the two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y
%
% Input:
%   - A1, B1, C1, A2, B2, C2: matrices
%   - opts: options (see below)
%
% Output: 
%   - lambda, mu: eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - XR, YR: components of decomposable right eigenvectors
%     (eigenvector is kron(XR(:,j),YR(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)*XR(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)*YR(:,j)=0
%   - XL, YL: components of decomposable left eigenvectors 
%     (eigenvector is kron(XL(:,j),YL(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)'*XL(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)'*YL(:,j)=0
%   - report: details of compression method in case of a singular problem 
%       rows are [mode m n r s(1) s(r) s(r+1) choice)], 
%       where mode is 1: CR step 1, 2: CR step 2, 3: RC step 1, 4. RC step 2
%       m,n: size of Delta matrices, r: rank
%       s(1), s(r), s(r+1): singular values
%       choice: how was rank determined (see NUMRANK for details)
%
% Operator determinants Delta0, Delta1, and Delta2 are used, where
% Delta0 = kron(C2, B1) - kron(B2, C1)
% Delta1 = kron(C2, A1) - kron(A2, C1)
% Delta2 = kron(A2, B1) - kron(B2, A1)
%
% Options in opts:
%   - singular (0): set to 1 for a singular problem, i.e., det(Delta0)=0
%   - rrqr (0): for singular problems, set to 1 to use rank revealing qr
%   - epscluster (1e-6): relative distance between eigenvalues in a cluster
%   - fast (1): use fast algorithm (can fail for multiple eigenvalues) 
%     or slow algorithm (0) with clustering
%   - inviter (1): use inverse iteration for eigenvectors or slow svd (0)
%   - all options of auxiliary functions
%   - novectors (0): set to 1 when report is important and vectors are not
%   - maxgensize (0): maximal size of the regular part (0 means no bound)
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - superior type of input data
%   - refine (1): Newton refinement steps to improve the accuracy of simple
%     eigenvalues of a regular 2EP - requires eigenvectors in output 
%   - rng_seed (0) : if different from zero, rng(rng_seed) is used at start
%   - all options of auxiliary functions
%   - test_shape (1): shows warning if simultaneous block triangularization
%     is not obtained, set to 0 for a slightly faster evaluation without test
%
% See also: TWOPAREIGS, TWOPAREIGS_JD, TWOPAREIGS_SI, THREEPAREIG, MULTIPAREIG

% Reference: M. E. Hochstenbach, T. Kosir, B. Plestenjak: A Jacobi-Davidson 
% type method for the two-parameter eigenvalue problem, SIAM J. Matrix Anal. 
% Appl. 26 (2005) 477-497

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 06.09.2015 : use extract_regular_part_np
% BP 24.08.2016 : double clustering in QZ (lambda first and then mu in blocks)
% BP 03.11.2016 : small speedup in computation od mu (inspired by Pavel Holoborodko's changes in multipareig)
% PH 22.11.2016 : modified to be precision-independent, i.e. be able to work with 
%                 any numeric type (whether it is built-in 'double'/'single' or custom like 'mp').
% BP 26.11.2016 : option fp_type
% PH 26.11.2016 : code simplifications, small fixes and clean-ups.
% BP 22.10.2018 : Newton refinement
% BP 03.04.2022 : random combination for joint triangularization
% Last revision: 03.04.2022

narginchk(6, 7);

% Analyse user supplied options, if any
if nargin < 7, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,A2,B2,C2);
end

if isfield(opts,'epscluster'), epscluster = opts.epscluster;   else, epscluster = numeric_t('1e-6',class_t);   end
if isfield(opts,'fast'),       fast       = opts.fast;         else, fast       = 1;                           end
if isfield(opts,'inviter'),    inviter    = opts.inviter;      else, inviter    = 1;                           end
if isfield(opts,'singular'),   singular   = opts.singular;     else, singular   = 0;                           end
if isfield(opts,'novectors'),  novectors  = opts.novectors;    else, novectors  = 0;                           end
if isfield(opts,'maxgensize'), maxgensize = opts.maxgensize;   else, maxgensize = 0;                           end
if isfield(opts,'refine'),     refine     = opts.refine;       else, refine     = 1;                           end
if isfield(opts,'rng_seed'),   rng_seed   = opts.rng_seed;     else, rng_seed   = 0;                           end
if isfield(opts,'test_shape'), test_shape = opts.test_shape;   else, test_shape = 1;                           end
refineeps = numeric_t('eps',class_t);

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end

% Default outputs
lambda = numeric_t([],class_t); 
mu     = numeric_t([],class_t); 
XR     = numeric_t([],class_t);  
YR     = numeric_t([],class_t);  
XL     = numeric_t([],class_t);  
YL     = numeric_t([],class_t); 
report = numeric_t([],class_t);

% if required, set the random generator
if rng_seed
    rng(rng_seed)
end

% Compute delta matrices 
[Delta0, Delta1, Delta2] = twopar_delta(A1,B1,C1,A2,B2,C2);

if singular
    [DeltaCell, report] = extract_regular_part_np({Delta0,Delta1,Delta2}, opts);
    Delta0 = DeltaCell{1}; Delta1 = DeltaCell{2}; Delta2 = DeltaCell{3}; 
end

% Quick return in degenerate cases
if ((size(Delta0,1) == 0) || (size(Delta0,1) ~= size(Delta0,2))) || ...  % no regular part was found
   ((maxgensize      > 0) && (size(Delta0,1)  > maxgensize))             % regular part is too large
     return
end

% Random rotation of (lambda,mu) to enable simultaneous triangularization
n = size(Delta0,1);
fi = rand(1,class_t)*numeric_t(pi,class_t);
c = cos(fi); s = sin(fi);
Delta1r =  c*Delta1 + s*Delta2;
Delta2r = -s*Delta1 + c*Delta2;

out_of_shape = 0;

if fast
    tmp = Delta0\[Delta1r Delta2r];
    Gamma1 = tmp(:,1:n);
    Gamma2 = tmp(:,n+1:end); 
    [Q1,D1] = schur(Gamma1,'complex');
    lambda = diag(D1);
    GQ = Gamma2*Q1;
    if test_shape
        QGQ = Q1'*GQ; 
        out_of_shape = norm(tril(QGQ,-1),'fro')/norm(QGQ,'fro');
        mu = diag(QGQ);
    else
        mu = zeros(n,1,class_t);
        for i = 1:n
            mu(i) = Q1(:,i)'*GQ(:,i);
        end
    end
else
    [S1,S0,Q,Z,order,start,csize,lambda] = clustered_qz(Delta1r,Delta0,epscluster); %#ok<*ASGLU>
    if max(csize)==1
        DZ = Delta2r*Z;
        if test_shape
            QDZ = Q*DZ;
            out_of_shape = norm(tril(QDZ,-1),'fro')/norm(QDZ,'fro');
            mu = diag(QDZ);
        else
            mu = zeros(n,1,class_t);
            for i = 1:n
                mu(i) = Q(i,:)*DZ(:,i);
            end
        end
        mu = mu./diag(S0);
    else
        S2 = Q*Delta2r*Z;
        if test_shape
            out_of_shape = out_of_shape_norm(S2,start,csize)/norm(S2,'fro');
        end
        for k = 1:length(start)
            if csize(k)==1
                mu = [mu; S2(start(k),start(k))/S0(start(k),start(k))];
            else
                partS0 = S0(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
                partS2 = S2(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
                tmpeig = eig(partS2,partS0);
                avgmu = sum(tmpeig)/csize(k)*ones(csize(k),1);
                mu = [mu; avgmu];
            end
        end
    end
end

% if out_of_shape>1e-06
%     warning(sprintf('Simultaneous triangularization relative error: %0.2g',out_of_shape))
% end

% Inverse rotation to obtain (lambda,mu)
lambdar = lambda; mur = mu;
lambda = c*lambdar - s*mur;
mu = s*lambdar + c*mur;

if (~novectors) && (nargout > 2) && all(isfinite(lambda))
    % extraction of eigenvectors (individually using inverse iteration or SVD)    
    
    % Generate initial vectors (in case of inverse iteration) only
    % once. This gives us a bit of speed-up, especially in extended precision
    % case, where multiple calls to randn might be noticeable.
    if inviter
        x10 = randn(size(A1,1),1,class_t);
        x20 = randn(size(A2,1),1,class_t);    
    else
        x10 = numeric_t([],class_t);
        x20 = numeric_t([],class_t);
    end
   
    for k = 1:length(lambda)
        [xr,xl] = min_sing_vec(A1-lambda(k)*B1-mu(k)*C1,inviter,x10,x10);
        XR(:,k) = xr; XL(:,k) = xl;
        [yr,yl] = min_sing_vec(A2-lambda(k)*B2-mu(k)*C2,inviter,x20,x20);
        YR(:,k) = yr; YL(:,k) = yl;
    end   
    if refine>0
        for k=1:length(lambda)
            [newlambda,newmu,newXR,newYR,res,flag] = twopareig_refine(A1,B1,C1,A2,B2,C2,lambda(k),mu(k),XR(:,k),YR(:,k),refine,refineeps);
            if flag
                lambda(k) = newlambda;
                mu(k) = newmu;
                XR(:,k) = newXR;
                YR(:,k) = newYR;
            end
        end
    end
end
end % twopareig

%------------------------------------------------------------------------

% auxiliary function return Frobenious norm of the part outside the block
% triangular shape
function off = out_of_shape_norm(A,start,csize)
    off = 0;
    for k=1:length(start)
        off = off + norm(A(start(k)+csize(k):end,start(k):start(k)+csize(k)-1),'fro')^2;
    end
    off = sqrt(off);
end