function [dx,dy] = bipolyder(P,x,y)

%bipolyder  Derivatives of a bivariate polynomial
%
% [dx,dy] = bipolyder(P,x,y) returns derivatives dp/dx and dp/dy of a 
% bivariate polynomial p(x,y) at x,y
%
% Polynomial is given with coefficients in a square matrix P so that 
% P(j,k) is a coefficient of p(x,y) at x^(j-1)*y^(k-1). The elements 
% below the antidiagonal are not important.

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 13.11.2015 Bor Plestenjak

n = size(P,1) - 1;

dx = bipolyval(diag(1:n)*P(2:n+1,1:n),x,y);
dy = bipolyval(P(1:n,2:n+1)*diag(1:n),x,y);
