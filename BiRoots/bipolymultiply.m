function Q = bipolymultiply(P,a,b,c)

%bipolymultiply  Multiply homogeneous bivariate polynomial by ax+by+cz
%
% Q = bipolymultiply(P,a,b,c) returns coefficients of a homogeneous 
% bivariate polynomial q(x,y,z), that we obtain by multiplying a 
% homogeneous bivariate polynomial p(x,y,z) by the linear factor 
% (ax+by+cz). Polynomial is given with coefficients in an (n+1)x(n+1) square 
% matrix P so that P(j,k) is a coefficient of p(x,y,z) at 
% x^(j-1)*y^(k-1)*z^(n-k-j+2). The elements below the antidiagonal are not 
% important.
%
% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 04.10.2016 Bor Plestenjak

n = size(P,1);
Q = zeros(n+1);

if c~=0 
    Q(1:n,1:n) = c*P;
end
if a~=0
    Q(2:n+1,1:n) = Q(2:n+1,1:n) + a*P;
end
if b~=0
    Q(1:n,2:n+1) = Q(1:n,2:n+1) + b*P;
end
