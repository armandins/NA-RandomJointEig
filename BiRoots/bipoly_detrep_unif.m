function [A,B,C] = bipoly_detrep_unif(P)

%bipoly_detrep_unif  Uniform representation of a bivariate polynomial
%
% [A,B,C] = bipoly_detrep_unif(P)
%
% returns matrices A,B,C such that det(A + x*B + y*C) = p(x,y), where
% polynomial p is given with coefficients in a square matrix P so that
% P(j,k) is a coefficient of p(x,y) at x^(j-1)*y^(k-1). The elements 
% below the antidiagonal are not important.
%
% The algorithm works for any polynomial of degree n and returns
% a representation with (2n-1) x (2n-1) matrices. No computation is 
% involved in the construction.
%
% Example: For a determinantal representation of the polynomial
%
% p(x,y) = 1 + 2x + 3y + 4x^2 + 5xy + 6y^2 + 7x^3 + 8x^2y + 9xy^2 + 10y^3
%
% use
%
% P = [ 1 3 6 10; 2 5 9 0; 4 8 0 0; 7 0 0 0 ];
% [A,B,C] = bipoly_detrep_unif(P)
%
% See also: bipoly_detrep

% Reference: A. Boralevi, J. van Doornmalen, J. Draisma, M. E. Hochstenbach, 
% B. Plestenjak: Uniform determinantal representations, arXiv.1607.04873

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 07.10.2016 Bor Plestenjak

narginchk(1,1)

n = size(P,2)-1;

if n==0
    A = P;
    B = 0;
    C = 0;
    return
end

M = rot90(triu(rot90(P(1:n,1:n),-1)));

A  = zeros(2*n-1);
B  = zeros(2*n-1);
C  = zeros(2*n-1);

A(1:n,1:n) = rot90(M,2);
d = diag(rot90(P));

B(1:n,1:n) = fliplr(diag(d(end:-1:2)));
C(n,1) = d(1);
A(n+1:end,1:n-1) = eye(n-1);
C(n+1:end,2:n) = -eye(n-1);
A(1:n-1,n+1:end) = eye(n-1);
B(2:n,n+1:end) = -eye(n-1);

if mod(n,2) == 0
    A = -A;
    B = -B;
    C = -C;
end