function [A,B,C,iter] = bipoly_detrep_nxn(P,opts)

%bipoly_detrep_nxn  Determinantal representation of a bivariate polynomial
%
% [A,B,C] = bipoly_detrep_nxn(P,opts)
%
% returns matrices A,B,C such that det(A + x*B + y*C) = p(x,y), where
% polynomial p is given with coefficients in a square matrix P so that
% P(j,k) is a coefficient of p(x,y) at x^(j-1)*y^(k-1). The elements 
% below the antidiagonal are not important. 
%
% The algorithm should return an n x n representation for a square-free 
% polynomials of degree n. However, due to numerical computation, it might 
% fail to return a representation good enough in the prescribed number of 
% attempts. Use it only for polynomials of degree <= 10.
%
% Example: For a determinantal representation of the polynomial
%
% p(x,y) = 1 + 2x + 3y + 4x^2 + 5xy + 6y^2 + 7x^3 + 8x^2y + 9xy^2 + 10y^3
%
% use
%
% P = [ 1 3 6 10; 2 5 9 0; 4 8 0 0; 7 0 0 0 ];
% [A,B,C] = bipoly_detrep_nxn(P)
%
% options in opts:
%    - max_tries: (3) maximum number of attempts with a random orthogonal
%                     transformation to compute a stable representation 
%    - max_linerr: (1e-10) maximum linearization error (using Monte Carlo
%                     type estimation)
%    - best_try: (1)  return attempt with minimal error (if all attempts
%                     have error larger than max_error)
% Output:
%    - A, B, C: matrices of the representation
%    - iter: -max_tries: when the method fails, 
%             0,1,...,max_tries-1: number of additional attemtps required, 
%             max_tries: the solution with the smallest error was returned
% 
% See also: bipoly_detrep

% Reference: B. Plestenjak: Minimal determinantal representations of 
% bivariate polynomials, arXiv.1607.03969

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 07.10.2016 Bor Plestenjak

narginchk(1,2)

if nargin<2, opts = [];  end
if isfield(opts,'max_tries'),  max_tries = opts.max_tries;   else max_tries = 5;     end
if isfield(opts,'max_linerr'), max_linerr = opts.max_linerr; else max_linerr = 1e-10; end
if isfield(opts,'best_try'),   best_try = opts.best_try;     else best_try = 1;      end

n = size(P,1)-1; % n is degree, size is n+1

minerr = Inf; 
A = []; B = []; C = [];

% if any of the bounds that should preserve numerical stability is
% exceeded, a random orthogonal transformation is applied up to max_tries 
% times before the method fails
flag = 1;
iter = 0;
while (flag~=0) && (iter<max_tries)
    if iter>0
        QS = orth(rand(3)); % random orthogonal transformation of (x,y,z)
        P1 = bipolysubstitution(P,QS');
    else
        P1 = P;
    end
    [A1,B1,C1,flag] = detrep_nxn(P1,opts);
    if flag==0
        linerror = test_linearization(P1,A1,B1,C1);
        if linerror<max_linerr
            if iter == 0
                A = A1; B = B1; C = C1;
            else
                A = QS(1,3)*B1 + QS(2,3)*C1 + QS(3,3)*A1;
                B = QS(1,1)*B1 + QS(2,1)*C1 + QS(3,1)*A1;
                C = QS(1,2)*B1 + QS(2,2)*C1 + QS(3,2)*A1;
            end
        else
            flag = 1;
            if best_try && (linerror < minerr)
                if iter == 0
                    A = A1; B = B1; C = C1;
                else
                    A = QS(1,3)*B1 + QS(2,3)*C1 + QS(3,3)*A1;
                    B = QS(1,1)*B1 + QS(2,1)*C1 + QS(3,1)*A1;
                    C = QS(1,2)*B1 + QS(2,2)*C1 + QS(3,2)*A1;
                end
                minerr = linerror;
            end
        end
    end
    iter = iter + 1;
end

if flag~=0
    if best_try && (minerr<Inf)
        iter = iter;
    else
        iter = -iter;
    end
else
    iter = iter - 1;
end

end

%-----------------------------------------------------------------------
%  Determinantal n x n representation of a generic bivariate polyomial
%-----------------------------------------------------------------------

function [A,B,C,flag] = detrep_nxn(P,opts)

flag = 0;

if nargin<2
    opts=[];
end

if isfield(opts,'defshift'),  defshift = opts.defshift;  else defshift = 0.5;   end

n = size(P,2);
A = zeros(n-1); B = zeros(n-1); C = zeros(n-1);

if n==1
    A = P; B = 0; C = 0;
    return
end

c = diag(rot90(P,-1)); % Coefficients of terms with max degree on antidiagonal
shift = 0;
% if the coefficient at x^n is small, we perform the shift
if abs(c(1)) < 1e-3*norm(c(2:end))
    shift = defshift;
    if shift ~= 0
        P = bipolyshift(P, shift, 1);
        c = diag(rot90(P,-1));
    end
end

if abs(P(1,n))>1e-14*norm(P)
    % we put coefficients at y^n and y^(n-1) to zero
    r = roots(c);                        % Root(s) of cubic terms
    derv = polyval(polyder(c),r);        % values of derivative in roots
    % Smallest real root
    [tmp, index] = min(abs(r));
    tau = r(index);
    P = bipolyshift(P, tau, 3);          % x -> x + ry to eliminate y^n
    if abs(derv(index)) < 1e-4
        % multiple zero, shift does not exist
        flag = 1;
        return
    end
else
    tau = 0;
end

sigma = -P(1,n-1) / P(2,n-1);
P = bipolyshift(P, sigma, 4);  % x -> x + c  to eliminate y^(n-1)
c = diag(rot90(P,-1));

zeta = roots(c(1:end-1));  % roots of x-y polynomial (without zero)
derv = polyval(polyder(c(1:end-1)),zeta);  % values of derivative in roots
if min(abs(derv)) < 1e-4
     % multiple zero, shift does not exist
     flag = 1;
     return
end
[tmp,ord] = sort(-abs(zeta));
zeta = zeta(ord);
lead = c(1);             % leading coefficient (at x^n)         
Fa = [0; ones(n-2,1)];
Fb = [0; -zeta];
gamma = zeros(1,n-1);

% initial polynomials q_0,...,q_(n-1) (indices go from 1 to n) 
Q = zeros(n-1,n-1,n-1);

Q(1,1,1) = 1;
for j = 2:n-1
    Q(1:j,1:j,j) = Fa(j,1)*[zeros(1,j); Q(1:j-1,1:j-1,j-1) zeros(j-1,1)]...
        + Fb(j,1)*[zeros(j-1,1) Q(1:j-1,1:j-1,j-1); zeros(1,j)];
end

% first residual (degree reduced to n-1)
R = P - lead*[zeros(1,n); Q(:,:,n-1) zeros(n-1,1)];

% main loop
for k = 2:n-2 % n-1
    % we compute alpha and beta 
    c = diag(rot90(R,-1),k-1);  % coefficient of degree n-k+1
    h = c(1:end-1)/lead;      % coefficients of polynomial h
    Fa(k+1:n-2,k) = ones(n-k-2,1);
    Fa(n-1,k) = h(1) - n + k + 2;
    
    W = zeros(n-1-k);
    for m = 1:n-1-k
        for ell = 1:m
            z = 1;
            for i = 1:ell-1
                z = z*(zeta(m)-zeta(i));
            end
            for j = ell+k:n-2
                z = z*(zeta(m)-zeta(j));
            end
            W(m,ell) = z;
        end
    end
    b = zeros(n-1-k,1);
    tmp = polyval(h,zeta(1:n-1-k));
    for m = 1:n-1-k
        b(m) = tmp(m);
        for ell = 1:m
            if ell < n-1-k
                b(m) = b(m) - zeta(m)*W(m,ell);
            else
                b(m) = b(m) - Fa(n-1,k)*zeta(m)*W(m,ell);
            end
        end
    end
    Fb(k+1:n-1,k) = W\b;

    % new computation of polynomials Q 
    for j = k+1:n-1
        Q(1:j,1:j,j) = zeros(j);
        Q(2:j,1:j-1,j) = Fa(j,1)* Q(1:j-1,1:j-1,j-1);
        Q(1:j-1,2:j,j) = Q(1:j-1,2:j,j) + Fb(j,1)*Q(1:j-1,1:j-1,j-1);
        for ell = 2:min(k,j-1)
            Q(2:j-ell+1,1:j-ell,j) = Q(2:j-ell+1,1:j-ell,j) + Fa(j,ell)*Q(1:j-ell,1:j-ell,j-ell);
            Q(1:j-ell,2:j-ell+1,j) = Q(1:j-ell,2:j-ell+1,j) + Fb(j,ell)*Q(1:j-ell,1:j-ell,j-ell);
        end
    end
    
    % new residual r  (degree reduced to n-k)
    R = P - lead*[zeros(1,n); Q(:,:,n-1) zeros(n-1,1)];
    for ell = 2:k-1
        R(1:n-ell,1:n-ell) = R(1:n-ell,1:n-ell) - gamma(n-1-ell)*Q(1:n-ell,1:n-ell,n-ell);
    end
    
    % gamma to set y^(n-k) to zero
    gamma(n-1-k)  = (-1)^(n-1-k)*R(1,n-k);
    for j = 1:n-1-k
        gamma(n-1-k) = gamma(n-1-k)/zeta(j);
    end
      
    % updated residual
    R(1:n-k,1:n-k) = R(1:n-k,1:n-k) - gamma(n-1-k)*Q(1:n-k,1:n-k,n-k);

end

% Now we write the matrices
A = diag([R(1,1) ones(1,n-2)]);
A(1,2:n-2) = gamma(1:n-3);

B(1,1) = R(2,1);
B(1,n-1) = lead;

for j = 2:n-1
    for ell = 1:j-1
        B(j,ell) = -Fa(j,j-ell);
        C(j,ell) = -Fb(j,j-ell);
    end
end
C(n-1,n-1) = 0;

A = A-sigma*B;
C = C-tau*B;

if shift~=0
  B = B - shift*C;
end

end

%-----------------------------------------------------------------------
%  Test of the linearization
%-----------------------------------------------------------------------

function err = test_linearization(P,A,B,C)

m = 200;

x = 0.5-rand(m,1)+1i*(0.5-rand(m,1));
y = 0.5-rand(m,1)+1i*(0.5-rand(m,1));
d1 = zeros(m,1);
d2 = zeros(m,1);
d2 = bipolyval(P,x,y);
for k = 1:m
    d1(k) = det(A+x(k)*B+y(k)*C);
end
err = max(abs(d1-d2)./(1e-4+abs(d2)));

end