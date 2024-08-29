function P = bipolyshift(P,a,ver)

%bipolyshift  Shift a bivariate polynomial
%
% P = bipolyshift(P,a,ver) returns coefficients of a shifted bivariate
% polynomial. Polynomial is given with coefficients in a square matrix P so 
% that P(j,k) is a coefficient of p(x,y) at x^(j-1)*y^(k-1). The elements 
% below the antidiagonal are not important.
%
% ver: 1: y -> y + ax
%      2: y -> y + a
%      3: x -> x + ay
%      4: x -> x + a

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 24.11.2015 Bor Plestenjak

n = size(P,1);

len = 1;
c = zeros(1,n);
c(1) = 1;

if ver == 1 || ver == 3
  if ver == 3
    P = P.';
  end
  Q = rot90(P);
  for m = 2:n
    c(:,2:m) = c(:,2:m) + a*c(:,1:m-1);
    for i = 1:n-m+1
      v = P(i,m) * c(:,2:m);
      v = [zeros(1,i) v];
      Q = Q + diag(v, m-n+i-1);
    end
  end
  P = rot90(Q, -1);
  if ver == 3
    P = P.';
  end
else
  if ver == 2
    P = P.';
  end
  for i = 2:n
    c(:,2:i) = c(:,2:i) + a*c(:,1:i-1);
    for j = 1:n+1-i
      v = (P(i,j) * c(:,i:-1:2)).';
      P(1:i-1,j) = P(1:i-1,j) + v; 
    end
  end
  if ver == 2
    P = P.';
  end
end
