function v = bipolyval(P,x,y)

%bipolyval  Value of a bivariate polynomial
% 
% v = bipolyval(P,x,y) returns values of a bivariate polynomial P at x,y
%
% Polynomial is given with coefficients in a square matrix P so that 
% P(j,k) is a coefficient of p(x,y) at x^(j-1)*y^(k-1). The elements 
% below the antidiagonal are not important.

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 13.11.2015 Bor Plestenjak

n = size(P,1);

[m1,m2] = size(x);

z = ones(m1,m2);
v = zeros(m1,m2);

for i = 1:n
    v = v + z.*polyval(P(i,n+1-i:-1:1),y);
    z = z.*x;
end


