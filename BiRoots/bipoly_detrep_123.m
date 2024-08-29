function [A,B,C] = bipoly_detrep_123(P,method,opts)

%bipoly_detrep_123  Determinantal representation of a bivariate polynomial
%
% [A,B,C] = bipoly_detrep_123(P,method,opts)
%
% returns matrices A,B,C such that det(A + x*B + y*C) = p(x,y), where
% polynomial p is given with coefficients in a square matrix P so that
% P(j,k) is a coefficient of p(x,y) at x^(j-1)*y^(k-1). The elements 
% below the antidiagonal are not important.
%
% Example: For a determinantal representation of the polynomial
%
% p(x,y) = 1 + 2x + 3y + 4x^2 + 5xy + 6y^2 + 7x^3 + 8x^2y + 9xy^2 + 10y^3
%
% use
%
% P = [ 1 3 6 10; 2 5 9 0; 4 8 0 0; 7 0 0 0 ];
% [A,B,C] = bipoly_detrep_123(P)
%
% method: 1: linearization 1  (matrices of size n^2/4) - no computation
%         2: (default) linearization 2 (matrices of size n^2/6) - computation
%         3: linearization 2 with special construction for degree 3 or 4 
%            (matrices of size n^2/6) - computation
%
% options in opts:
%    - threshold: (1e-14) for linerization 2
%    - defshift: (0.5) shift used if leading coeffient is small in linearization 2
%
% See also: biroots

% Reference: B. Plestenjak, M. E. Hochstenbach: Roots of bivariate 
% polynomial systems via determinantal representations, SIAM J. Sci. 
% Comput. 38 (2016) A765-A788

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 29.09.2016 Bor Plestenjak

narginchk(1,3)

if nargin < 2 || isempty(method)
    method = 2;
end
if nargin < 3
    opts = [];
end
if method==3
    opts.reduce34 = 1;
end

n = size(P,2);
b = 0;                      % Beginning index
B = [];
C = [];

if method == 1 % Linearization 1
    while n > 0
        A(1,b+1) = P(1,1);      % Coefficients polynomial
        for i = 1:n-1
            B(1,b+i) = P(i+1,1);
            C(1,b+i) = P(i,2);
            if i < n-1            % Coupling x
                B(b+i+1,b+i)   =  1;
                A(b+i+1,b+i+1) = -1;
            end
        end
        P = P(1:n-2, 3:n);      % Extracting y^2
        if n > 2                % Coupling y
            C(b+n,b+1) =  1;
            A(b+n,b+n) = -1;
        end
        if n == 3
            C(1,b+n) = P(1,1);
            break
        end
        if n > 3                % Coupling y^2
            C(b+n+1,b+n)   =  1;
            A(b+n+1,b+n+1) = -1;
        end
        b = b+n;
        n = n-2;
    end
    % Final update
    m  = max([size(A) size(B) size(C)]);
    A0 = zeros(m); A0(1:size(A,1),1:size(A,2)) = A; A = A0;
    B0 = zeros(m); B0(1:size(B,1),1:size(B,2)) = B; B = B0;
    C0 = zeros(m); C0(1:size(C,1),1:size(C,2)) = C; C = C0;
    if mod(size(A,1),2)==0
        A(1,:) = -A(1,:); 
        B(1,:) = -B(1,:); 
        C(1,:) = -C(1,:);
    end
else  % Linearization 2
    [A, B, C] = bipolylin2(P,opts);
end

%-----------------------------------------------------------------------
%  Linearization 2 for bivariate polyomial
%-----------------------------------------------------------------------

function [A, B, C] = bipolylin2(P,opts)

% Determinantal representation of a bivariate polynomial using linearization 2 
%
% Last revision 24.11.2015 Bor Plestenjak

if nargin<2
    opts=[];
end

if isfield(opts,'threshold'), threshold = opts.threshold;  else threshold = 1e-14; end
if isfield(opts,'reduce34'),  reduce34 = opts.reduce34;    else reduce34 = 0;     end
if isfield(opts,'defshift'),  defshift = opts.defshift;    else defshift = 0.5;   end

n = size(P,2);

% we reduce the order of the polynomial if all leading coefficients are small
c = diag(rot90(P,-1)); % Coefficients of terms with max degree on antidiagonal
while (n>1) && (max(abs(double(c))) < threshold)
    P = P(1:n-1,1:n-1);
    n = n-1;
    c = diag(rot90(P,-1));
end

% special case for constant polynomial
if n==1
    A = P(1,1);
    B = 0;
    C = 0;
    return
end

% special case for linear polynomial
if n==2
    A = P(1,1);
    B = P(2,1);
    C = P(1,2);
    return
end

% special case (3x3 matrices) for cubic polynomial
if (n==4) && reduce34
    Pini = P;
    % if leading coefficient is zero or small, we have to swap
    q = diag(rot90(P, -1));  % Cubic terms
    swap = abs(q(1))<1e-1*norm(q(2:end));
    if swap
        P = P.';
    end
    q = diag(rot90(P, -1));                       % Cubic terms
    r = roots(q);                                 % Root(s) of cubic terms
    derv = polyval(polyder(q),r);                 % values of derivative in roots 
    [maxder, index] = sort(-abs(derv));
    r = r(index);                                 % Smallest real root
    P = bipolyshift(P, r(1), 3);                  % x -> x + ry to eliminate y^3
    if abs(maxder(1)) < 1e-4
        % multiple zero, shift does not exist
        opts.reduce34=0;
        % disp('Shift is multiple zero, reverting to linerization 2..')
        [A, B, C] = bipolylin2(Pini,opts);
        return
    end
    c = -P(1,3) / P(2,3);
    P = bipolyshift(P, c, 4);                    % x -> x + c  to eliminate y^2
    
    % Compute determinantal expression of P
    A(1,1) = P(1,1);                              % Constant and linear terms: a + bx + cy
    B(1,1) = P(2,1);
    C(1,1) = P(1,2);
    B(1,2) = P(3,1);                              % x^2 term
    C(1,2) = P(2,2);                              % xy term
    B(1,3) = P(4,1);                              % x^3 term
    C(1,3) = P(3,2) - B(1,3)*(-r(2)+r(1));        % x^2y term
    B(2,1) = -1;                                  % For x in vector
    A(2,2) =  1;
    B(3,2) = -1;                                  % For x-r in vector
    C(3,2) =  r(2)-r(1);
    A(3,3) =  1;
    
    % Resize
    A0 = zeros(3); A0(1:size(A,1),1:size(A,2)) = A; A = A0;
    B0 = zeros(3); B0(1:size(B,1),1:size(B,2)) = B; B = B0;
    C0 = zeros(3); C0(1:size(C,1),1:size(C,2)) = C; C = C0;
    
    % Undo shifts
    A = A-c*B;
    C = C-r(1)*B;
    % Undo swap
    if swap
        tmp = B; B = C; C = tmp;
    end
    return
end

% special case (4x4 matrices) for quartic polynomial
if (n==5) && reduce34
    Pini = P;
    % if leading coefficient is zero, we have to swap
    q = diag(rot90(P, -1));                       % Quartic terms
    swap = abs(q(1))<1e-1*norm(q(2:end));
    if swap
        P = P.';
    end
    q = diag(rot90(P, -1));             % Quartic terms
    r = roots(q);                       % roots of quartic terms
    derv = polyval(polyder(q),r);       % values of derivative in roots 
    [maxder, index] = sort(-abs(derv));
    r = r(index);                       % Root with the largest value of derivative
    P = bipolyshift(P, r(1), 3);        % x -> x + ry to eliminate y^4
    if abs(maxder) < 1e-4
        % multiple zero, shift does not exist
        opts.reduce34=0;
        % disp('Shift is multiple zero, reverting to linerization 2..')
        [A, B, C] = bipolylin2(Pini,opts);
        return
    end
    c = -P(1,4) / P(2,4);
    P = bipolyshift(P, c, 4);                     % x -> x + c  to eliminate y^3
    
    q = diag(rot90(P, -1)).';                     % Quartic terms
    q = fliplr(q(1:4));
    r2 = roots(q);                                % Roots of quartic terms
    derv = polyval(polyder(q),r2);                % values of derivative in roots 
    [maxder, index] = sort(-abs(derv));
    r2 = r2(index);                               % Smallest real root
    P = bipolyshift(P, r2(1), 1);                 % y -> y + rx to eliminate x^4
    if abs(maxder) < 1e-4
        % multiple zero, shift does not exist
        opts.reduce34=0;
        % disp('Shift is multiple zero, reverting to linerization 2..')
        [A, B, C] = bipolylin2(Pini,opts);
        return
    end
    c2 = -P(4,1) / P(4,2);
    P = bipolyshift(P, c2, 2);                   % y -> y + c  to eliminate x^3
    q = diag(rot90(P, -1)).';                    % Quartic terms
    r3 = roots(q(2:4));                          % Roots of quartic terms
    [~, index] = sort(abs(r3));
    r3 = r3(index);
    
    % Compute determinantal representation of P
    A(1,1) = P(1,1);                             % Constant and linear terms: a + bx + cy
    B(1,1) = P(2,1);
    C(1,1) = P(1,2);
    B(1,2) = P(3,1);                              % x^2 term
    C(1,2) = P(2,2);                              % xy term
    B(1,3) = P(3,2);                              % x^2y term
    C(1,3) = P(2,3);                              % xy^2 term
    C(1,3) = P(2,3);                              % x^2y term
    B(2,1) = -1;                                  % For x in vector
    A(2,2) =  1;
    C(3,2) = -1;                                  % For xy in vector
    A(3,3) =  1;
    B(4,3) = -1;                                  % For xy(x-r) in vector
    C(4,3) =  r3(1);
    A(4,4) =  1;
    A(5,5) =  1;                                  % For y in vector
    C(5,1) = -1;
    C(1,5) = P(1,3);                              % For quintic terms
    q1 = deconv(q(2:4), [1 -r3(1)]);
    B(1,4) = q1(1);
    C(1,4) = q1(2);
    
    % Resize
    A0 = zeros(5); A0(1:size(A,1),1:size(A,2)) = A; A = A0;
    B0 = zeros(5); B0(1:size(B,1),1:size(B,2)) = B; B = B0;
    C0 = zeros(5); C0(1:size(C,1),1:size(C,2)) = C; C = C0;
    
    % Undo shifts
    A = A-c2*C;
    B = B-r2(1)*C;
    A = A-c*B;
    C = C-r(1)*B;
    % Undo swap
    if swap
        tmp = B; B = C; C = tmp;
    end
    return
end

% we know that polynomial is at least of order 2
A = zeros(n);
B = zeros(n);
C = zeros(n);
shift = 0;
r = roots(c);
% if we need a shift we do it and update the roots

if abs(c(1)) < 1e-2*max(abs(c(2:end)))
    shift = defshift;  
    if shift ~= 0
        P = bipolyshift(P, shift, 1);
        % the roots are adjusted to the shift
        if ~isempty(r)
            r = [r./(1-shift*r); -1/shift*ones(n-1-length(r),1)];
        else
            r = -1/shift*ones(n,1);
        end
    end
end

A(1,1) = P(1,1);        % Constant term P00
[~, index] = sort(abs(r));
r = r(index);           % Reorder the roots
for k = 2:n-1           % Couplings for the first branch
   B(k,k-1) = -1;
   C(k,k-1) =  r(k-1);
   A(k,k) =  1;
end
B(1,1) = P(2,1);        % Linear terms P10*x + P01*y
C(1,1) = P(1,2);
q = [1 -r(1)];          % Match terms of the form a x^j + b x^(j-1) y
for k = 2:n-2
   B(1,k) = P(k+1,1) / q(1);
   C(1,k) = (P(k,2) - B(1,k)*q(2)) / q(1);
   for j = 3:length(q)
      P(k+2-j,j) = P(k+2-j,j) - q(j)*B(1,k);
   end
   for j = 3:length(q)+1
      P(k+2-j,j) = P(k+2-j,j) - q(j-1)*C(1,k);
   end
   q = conv(q, [1 -r(k)]);
end
B(1,n-1) = P(n,1);
C(1,n-1) = -P(n,1)*r(n-1);

% The remainder is in P(1:n-3, 3:n-1)
% if the remainder is not zero, we need additional node y to link 
% the respresentation tree of the remainder to the main tree

remP = P(1:n-3, 3:n-1);
remP = rot90(triu(rot90(remP,-1)),1); % explicitly put last diagonal to zero

if norm(remP,'fro')>threshold*norm(P,'fro')
    % build a subtree for the remainder recursively
    [AS, BS, CS] = bipolylin2(remP,opts);
    ns = size(AS,1);
    % add link node y 
    C(n,1) = -1;
    A(n,n) =  1;
    % if the remainder is just a constant we put it to y node
    if (ns==1) && (BS(1,1)==0) && (CS(1,1)==0)
        C(1,n) = AS(1,1);
    else
       % we join the tree
       % add y^2 node
       C(n+1,n)   = -1;
       A(n+1,n+1) =  1;
       B(n+1,n+1) =  0;
       % add the subtree
       if ns > 1
         A(n+2:n+ns,n+1:n+ns) = AS(2:ns,:);
         B(n+2:n+ns,n+1:n+ns) = BS(2:ns,:);
         C(n+2:n+ns,n+1:n+ns) = CS(2:ns,:);
       end
       A(1,n+1:n+ns) = AS(1,:);
       B(1,n+1:n+ns) = BS(1,:);
       C(1,n+1:n+ns) = CS(1,:);
    end
else
  % remainder is zero
  A = A(1:n-1,1:n-1);
  B = B(1:n-1,1:n-1);
  C = C(1:n-1,1:n-1);
end

% shift back
if shift~=0
  B = B - shift*C;
end