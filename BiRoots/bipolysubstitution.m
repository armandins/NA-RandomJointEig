 function Q = bipolysubstitution(P,M)

%bipolysubstitution  Substitution of homogeneous bivariate polynomial
%
% Q = bipolysubstitution(P,M) returns coefficients of a homogeneous bivariate
% polynomial q(x,y,z), that we obtain by a linear substitution defined by
% 3x3 transition matrix M on polynomial p(x,y,z). Polynomial is given with 
% coefficients in an (n+1)x(n+1) square matrix P so that P(j,k) is a 
% coefficient of p(x,y,z) at x^(j-1)*y^(k-1)*z^(n-k-j+2). The elements 
% below the antidiagonal are not important.
%
% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 04.10.2016 Bor Plestenjak

n = size(P,1) - 1;
% this could be optimized
Q = zeros(n+1);
for j = 0:n
    for k = 0:n-j
        T = 1;
        for i = 1:j
            T = bipolymultiply(T,M(1,1),M(1,2),M(1,3));  % new x
        end 
        for i = 1:k
            T = bipolymultiply(T,M(2,1),M(2,2),M(2,3));  % new y
        end
        for i = 1:n-j-k
            T = bipolymultiply(T,M(3,1),M(3,2),M(3,3));  % new z
        end
        Q = Q + P(j+1,k+1)*T;
    end
end


