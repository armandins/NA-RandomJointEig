function [x,y,iter,lin_iter] = biroots(P1,P2,method,opts)

%biroots Roots of a system of two bivariate polynomials
%
% [x,y,iter,lin_iter] = biroots(P1,P2,method) returns solutions of a system 
% of bivariate polynomials p1(x,y) = p2(x,y) = 0. Polynomials are given 
% with coefficients in square matrices P1 and P2 so that Pi(j,k) is a 
% coefficient of pi(x,y) at x^(j-1)*y^(k-1). The elements below the 
% antidiagonal are not important. 
%
% Example: To solve the system
%
% 1 + 2x + 3y + 4x^2 + 5xy + 6y^2 + 7x^3 + 8x^2y + 9xy^2 + 10y^3 = 0
% 10 + 9x + 8y + 7x^2 + 6xy + 5y^2 + 4x^3 + 3x^2y + 2xy^2 +  y^3 = 0
%
% use
%
% P1 = [ 1 3 6 10; 2 5 9 0; 4 8 0 0; 7 0 0 0];
% P2 = [10 8 5  1; 9 6 2 0; 7 3 0 0; 4 0 0 0];
% [x,y] = biroots(P1,P2)
%
% method: 
%   - 1: linearization 1 (matrices of size n^2/4) - no computation [1]
%   - 2: linearization 2 (matrices of size n^2/6) - computation [1]
%   - 3: linearization 2 with special construction for degree 3
%        and 4 (matrices of size n^2/6) - computation [1]
%   - 4: bipoly_detrep_quantic - size n for degree 5 or less [2]
%   - 5: bipoly_detret_nxn - size n for a square-free polynomial [3]
%   - 6: bipoly_detrep_unif - matrices of size 2n-1 - no computation [4]  
%   - 0: (default) not specified, default method is 4 for degree 5 or 
%        less, 5 for degree 9 or less, and 6 otherwise
%
% options in opts:
%    - normalize: normalize polynomials by dividing them with the largest
%         coefficient (1 - default), 0: no normalization
%    - refine: number of Newton refinement steps in the end (5 - default),
%         set to 0 in case of multiple roots
%    - fast: use fast solver for the singular MEP (1 - default), set to 0
%         if there are some roots with the same x component
%    - generic: (1 - default) the number of roots is equal to the product
%         of degrees of p1 and p2. In this case the method restarts with
%         different settings if the number of computed roots is different.
%         Set to 0 if this is not a generic pair and the number of roots is
%         smaller than the product of the degrees.
%    - all options of twopareig from MultiParEig toolbox
%
% Output
%     - x,y: the roots
%     - iter: number of additional restarts (-1 when the method fails)
%     - lin_iter : number of restarts in the linearization
%
% See also: bipoly_detrep

% References: 
%   [1] B. Plestenjak, M. E. Hochstenbach: Roots of bivariate 
%       polynomial systems via determinantal representations, SIAM J. Sci. 
%       Comput. 38 (2016) A765-A788
%   [2] A. Buckley, B. Plestenjak: Simple determinantal representations of 
%       up to quintic bivariate polynomials, arXiv.1609.00498
%   [3] B. Plestenjak: Minimal determinantal representations of 
%       bivariate polynomials, arXiv.1607.03969
%   [4] A. Boralevi, J. van Doornmalen, J. Draisma, M. E. Hochstenbach, 
%       B. Plestenjak: Uniform determinantal representations, 
%       arXiv.1607.04873

% This package requires package MultiParEig, available at
% http://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 16.10.2016 Bor Plestenjak

narginchk(2,4)

if nargin<3, method = 0; end
if nargin<4, opts = [];  end

if isfield(opts,'normalize'),  normalize = opts.normalize; else normalize = 1; end
if isfield(opts,'refine'),     refine = opts.refine;       else refine = 5;    end
if ~isfield(opts,'fast'),      opts.fast = 0;      end
if ~isfield(opts,'generic'),   opts.generic = 1;   end
if isfield(opts,'mingap'),     mingap = opts.mingap;       else mingap = 1;    end
if isfield(opts,'gapmargin'),  gapmargin = opts.gapmargin; else gapmargin = 2;    end
if isfield(opts,'maxiter'),    maxiter = opts.maxiter;     else maxiter = 6;   end
defrankeps = (~isfield(opts,'rankeps'));

if normalize
    P10 = P1/max(max(abs(P1)));
    P20 = P2/max(max(abs(P2)));
else
    P10 = P1;
    P20 = P2;
end

bestx = [];
besty = [];
bestgap = 0;

% If the method fails to find solutions (or the exact number of solutions
% when opts.generic=1), it attempts again up to three times with shifted
% polynomials and interchanged roles of x and y. In each additional attempt 
% the tolerance for the staircase algorithm in twopareig from MultiParEig 
% is reduced. Please note that even this heuristics does not always lead 
% to solutions for some challenging problems.

n1 = size(P1,2) - 1;
n2 = size(P2,2) - 1;

success = 0; 
run = 0; 
tran = 0;
alpha = [0 0 rand(1,max(1,maxiter-2))];      % shifts for additional attempts
beta =  [0 0 rand(1,max(1,maxiter-2))];      % shifts for additional attempts
defshifts = [0.5 0.5 rand(1,max(1,maxiter-2))];    % default shifts for the 2nd linearization algorithm
rankepsilon = [1e-11 1e-11 1e-10 1e-10 1e-9*ones(1,max(1,maxiter-4))]; % tolerance for twopareig
opts.novectors = 1;

if (method==0) && (max(n1,n2)>9)
    method = 6;
end

while (run<maxiter) && ~success
    run = run + 1;
    tran = mod(run+1,2); % for each second attempt we switch variables x and y 
    opts.defshift = defshifts(run);
    if defrankeps
        opts.rankeps = rankepsilon(run);
    end
    if tran 
        P1S = P10.';
        P2S = P20.';
    else
        P1S = P10;
        P2S = P20;
    end
    if run>2
        % for additional attempts we shift the polynomials
        P1S = bipolyshift(P1S,alpha(run),1);
        P2S = bipolyshift(P2S,alpha(run),1);
        P1S = bipolyshift(P1S,beta(run),3);
        P2S = bipolyshift(P2S,beta(run),3);
    end
    [A1, B1, C1, met1, flag1] = bipoly_detrep(P1S,method,opts);
    % [A1, C1, B1, met1, flag1] = bipoly_detrep(P1S',method,opts); % alternative linearization of P1 with swapped x and y 
    [A2, B2, C2, met2, flag2] = bipoly_detrep(P2S,method,opts);
    % [A2, C2, B2, met2, flag2] = bipoly_detrep(P2S',method,opts); % alternative linearization of P2 with swapped x and y 
    lin_iter = flag1 + flag2;
    if (flag1>=0) && (flag2>=0)
        m = size(A1,1)*size(A2,1);
        if m==n1*n2
            Delta0 = kron(B1,C2) - kron(C1,B2);
            opts.singular = rank(Delta0)<m;
        else
            opts.singular = 1;
        end
        % The roots are eigenvalues of the corresponding two-parameter eigenvalue problem
        if (met1==6) && (met2==6) && opts.generic && (max(size(P10,1),size(P20,1))>2)
            % In case of uniform representations we have heuristic formula
            % for the ranks in the staircase algorithm (for the generic
            % pair of bivariate polynomials)
            RS = generic_sizes_unif(size(P10,1)-1,size(P20,1)-1);
            opts.ranksequence = [ones(size(RS,1),1) RS(:,3:4)];
        else
            opts.ranksequence = [];
            opts = rmfield(opts,'ranksequence');
        end
        opts.maxgensize = n1*n2;
        [x,y,tilda1,tilda2,tilda3,tilda4,hist] = twopareig(A1,-B1,-C1,A2,-B2,-C2,opts);   
        if opts.singular 
            if (length(x)>0) && (length(x)<=n1*n2)
                ming = log10(min(hist(:,6)./hist(:,7)));
                if ming>bestgap
                    % we memorize the solution with the minimal rank gap
                    tmp = [x y]*[1 alpha(run);beta(run) 1+alpha(run)*beta(run)];
                    bestx = tmp(:,1);
                    besty = tmp(:,2);
                    bestgap = ming;
                end
                success = ming>gapmargin;
            else
                success = 0;
            end
        else
            success = 1;
            tmp = [x y]*[1 alpha(run);beta(run) 1+alpha(run)*beta(run)];
            bestx = tmp(:,1);
            besty = tmp(:,2);
        end
    end
    if ~success && (method==0)
        % if the method is not fixed, we switch to minunif in next runs 
        method = 6;
    end
    
end

if tran
    y = bestx;
    x = besty;
else
    x = bestx;
    y = besty;
end

if ~isempty(x)
    [x,y] = biroots_newton(P1,P2,x,y,refine,1e-12);
    iter = run-1; % number of additional attempts
else
    iter = -1; 
end

end

%------------------------------------------------------------------------
%  Ranks in staircase algorithm when uniform determinantal representation 
%  is applied to a pair of generic bivariate polynomials
%------------------------------------------------------------------------

function M = generic_sizes_unif(n1,n2)

% There is no proof for it, but, by observing ranks in the staircase
% algorithm applied to a singular two-parameter eigenvalue problems from
% uniform determinantal representations of generic bivariate polynomials,
% one can find a pattern. This can be useful for more sensitive systems,
% where rank might be determined wrong otherwise.

n = min(n1,n2);
k = max(n1,n2) - n;

inisize = (2*n-1)*(2*(n+k)-1);
M = [];
r1 = 2*((n+k)*n-1);
r2 = 2*(n+k-1)*(n-1);
M = [inisize inisize r1 r2];

if n<=(n+k)/2
    st = 2*n-2;
    for j = 1:n-1
        M = [M;0 0 0 st];
        st = st - 1;
    end
    M = [M;zeros(k-n+1,3) n*ones(k-n+1,1)];
    st = n+1;
    for j = 1:n-1
        M = [M; 0 0 0 st];
        st = st - 1;
    end
else
    st = 2*n-2;
    for j = 1:k
        M = [M;0 0 0 st];
        st = st - 1;
    end
    st = st + 1;
    for j = 1:n-k-1
        M = [M; 0 0 0 st];
        st = st - 2;
    end
    st = k +2;
    for j = 1:k
        M = [M; 0 0 0 st];
        st = st - 1;
    end
end

dif = zeros(n+k-1,1);
dif(1) = 1;
if n<=(n+k)/2
    for j=1:n-1
        dif(j) = dif(j) - 1;
        dif(n+k-j) = 1;
    end
else
    for j=1:k
        dif(j) = dif(j) - 1;
        dif(n+k-j) = 1;
    end
end

for j = 2:n+k
    M(j,1) = M(j-1,1) - M(j-1,4);
    M(j,2) = M(j-1,3);
    M(j,3) = M(j,2) - M(j,4) + dif(j-1);
end

if n == 1
    M = M(2:end,:);
    M(1,2)  = M(1,1);
end

end

% ---------------------------------------------------------------------
% Final refinement with the Newton method
%------------------------------------------------------------------------

function [xn,yn] = biroots_newton(P1,P2,x,y,maxsteps,epsilon)

xf = [];
yf = [];

for step = 1:maxsteps
    v1 = bipolyval(P1,x,y);
    v2 = bipolyval(P2,x,y);
    [dx1, dy1] = bipolyder(P1,x,y);
    [dx2, dy2] = bipolyder(P2,x,y);
    
    xn = x;  
    yn = y;
    xc = []; 
    yc = [];
    
    for j=1:length(x)
        M = [dx1(j) dy1(j); dx2(j) dy2(j)];
        b = [v1(j); v2(j)];
        delta = M\b;
        xn(j) = xn(j) - delta(1);
        yn(j) = yn(j) - delta(2);
        if norm(delta)>=epsilon*norm([xn(j) yn(j)])
            % we refine this point in the next round again
            xc = [xc; xn(j)];
            yc = [yc; yn(j)];
        else
            % this root has converged
            xf = [xf; xn(j)];
            yf = [yf; yn(j)];
        end
    end
    % points for the next round
    x = xc;
    y = yc;
end
% all nonconverged points are added
xn = [xf;x];
yn = [yf;y];    
    
end