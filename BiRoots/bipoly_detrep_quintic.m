function [A,B,C,iter] = bipoly_detrep_quintic(P,opts)

%bipoly_detrep_quintic  Determinantal representation of a bivariate polynomial
%
% [A,B,C] = bipoly_detrep_quintic(P,opts)
%
% returns matrices A,B,C such that det(A + x*B + y*C) = p(x,y), where
% polynomial p is given with coefficients in a square matrix P so that
% P(j,k) is a coefficient of p(x,y) at x^(j-1)*y^(k-1). The elements 
% below the antidiagonal are not important. The degree of the polynomial
% must be 5 or less, i.e., P must be 6 x 6 or smaller.
%
% The algorithm returns an n x n representation for any polynomial of degree 
% n <= 5. However, due to numerical computation, it might fail to return a 
% representation good enough in the prescribed number of attempts. 
%
% Example: For a determinantal representation of the polynomial
%
% p(x,y) = 1 + 2x + 3y + 4x^2 + 5xy + 6y^2 + 7x^3 + 8x^2y + 9xy^2 + 10y^3
%
% use
%
% P = [ 1 3 6 10; 2 5 9 0; 4 8 0 0; 7 0 0 0 ];
% [A,B,C] = bipoly_detrep_quintic(P)
%
% options in opts:
%    - max_tries: (3) maximum number of attempts with a random orthogonal
%                     transformation to compute a stable representation 
%
% Output:
%    - A, B, C: matrices of the representation
%    - iter: -(max_tries-1): when the method fails, 
%             0,1,...,max_tries-1: number of additional attemtps required 
%
% See also: bipoly_detrep

% Reference: A. Buckley, B. Plestenjak: Simple determinantal 
% representations of up to quintic bivariate polynomials, arXiv.1609.00498

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% Last revision 07.10.2016 Bor Plestenjak

narginchk(1,2)

if nargin<2, opts = [];  end
if isfield(opts,'max_tries'),  max_tries = opts.max_tries; else max_tries = 3; end

n = size(P,1)-1; % n is degree, size is n+1
if (n<0) || (n>5) 
    error('Method works only for polynomials of degree <= 5')
end

QS = eye(3);
P1 = P;
[A1,B1,C1,flag] = bipoly_detrep_quintic_aux(P1);
iter = 1;
% if any of the bounds that should preserve numerical stability is
% exceeded, a random orthogonal transformation is applied up to max_tries 
% times before the method fails
while (flag~=0) && (iter<=max_tries)
    QS = orth(rand(3)); % random orthogonal transformation of (x,y,z)
    P1 = bipolysubstitution(P,QS');
    [A1,B1,C1,flag] = bipoly_detrep_quintic_aux(P1);
    iter = iter + 1;
end
if flag~=0
    A = [];
    B = [];
    C = [];
    iter = -(max_tries);
else
    B = QS(1,1)*B1 + QS(2,1)*C1 + QS(3,1)*A1;
    C = QS(1,2)*B1 + QS(2,2)*C1 + QS(3,2)*A1;
    A = QS(1,3)*B1 + QS(2,3)*C1 + QS(3,3)*A1;
    iter = iter-1;
end

%-----------------------------------------------------------------------
%  Determinantal representation for up to quintic bivariate polyomial
%-----------------------------------------------------------------------

function [A,B,C,fail] = bipoly_detrep_quintic_aux(P,fixed)

% fixed=0 is default. fixed=1 is used for recursive call where random
% orthogonal substitutions are not allowed and linearization is always 
% computed, no errors are reported

if nargin<2, fixed = 0; end

n = size(P,1)-1; % n is degree, size is n+1
Pnorm = norm(P,'fro');

A = []; B = []; C = []; fail = 0;

if n==0
    A = P(1,1);
    B = 0;
    C = 0;
    return
end

if n==1
    A = P(1,1); B = P(2,1); C = P(1,2);
    return
end

% (1) check if p = (ax+by+cz)^n
% ------------------------------------------------------------------
if n==5
    a = P(1,1)^(1/n);  b = P(n+1,1)^(1/n);  c = P(1,n+1)^(1/n);
    Q = [a c;b 0];
    for k = 1:n-1
        Q = bipolymultiply(Q,a,b,c);
    end
    if norm(P-Q)<1e-12*Pnorm
        A = a*eye(n); B = b*eye(n); C = c*eye(n);
        return
    end
end

% (2) if P(n+1,0) is very small, report fail=1 (if fixed = 0)
% ------------------------------------------------------------------
if (fixed==0) && (abs(P(n+1,1))<1e-6*Pnorm)
    fail = 1;
    return
end

% (3) compute alpha and beta
% ------------------------------------------------------------------
pxy = diag(rot90(P,-1));
alfa = roots(pxy);
[tmp,ord] = sort(-abs(alfa));
alfa = alfa(ord);
pxz = P(n+1:-1:1,1);
beta = roots(pxz);
[tmp,ord] = sort(-abs(beta));
beta = beta(ord);

% for n=4 we order alfa beta so that alfa(3) and beta(3) are minimal
if n==4
     alfa([3 end]) = alfa([3 end]);
     beta([3 end]) = beta([3 end]);
end
    
% for n=5 we make sure that (alfa3<>alfa4), (beta3<>beta4), and the
% intersection of lines is not on the polynomial curve
if n==5
    % first we take the most distant pairs alfa and beta
    [ai,aj] = max_dist(alfa);
    alfa([3 ai]) = alfa([ai 3]);
    alfa([4 aj]) = alfa([aj 4]);
    [bi,bj] = max_dist(beta);
    beta([3 bi]) = beta([bi 3]);
    beta([4 bj]) = beta([bj 4]);
    imen = alfa(3)-alfa(4);
    % intersection (x1,y1)
    x1 = (alfa(3)*beta(4)-alfa(4)*beta(3))/imen;
    y1 = (beta(4)-beta(3))/imen;
    absPxy1 = abs(bipolyval(P,x1,y1));
    % the second option is to interchange beta(3) and beta(4)
    x2 = (alfa(3)*beta(3)-alfa(4)*beta(4))/imen;
    y2 = -y1;
    absPxy2 = abs(bipolyval(P,x2,y2));
    if absPxy2 > absPxy1
        beta([3 4]) = beta([4 3]);
    end
    
    for iter = 1:3
        bull = 1;
        if (abs(alfa(3)-alfa(4))<1e-6*abs(alfa(3))) || (abs(beta(3)-beta(4))<1e-6*abs(beta(3)))
            bull = 0;
        end
        if bull
            imen = alfa(3)-alfa(4);
            x = (alfa(3)*beta(4)-alfa(4)*beta(3))/imen;
            y = (beta(4)-beta(3))/imen;
            absPxy = abs(bipolyval(P,x,y));
            if absPxy < 1e-8*Pnorm
                bull = 0;
            end
        end
        if bull
            break
        else % we try new random order (we try this up to three times)
            [tilda,ord1] = sort(rand(5,1));
            alfa = alfa(ord1);
            [tilda,ord2] = sort(rand(5,1));
            beta = beta(ord2);
        end
    end
    if bull==0
        fail = 3;
        return
    end
    Q5 = [alfa(3)*beta(4)-alfa(4)*beta(3), beta(4)-beta(3)  alfa(3)-alfa(4); ...
            1 -alfa(3) -beta(3); 1 -alfa(4) -beta(4)];
    if (fixed==0) && (cond(Q5)>200)
        fail = 4;
        return
    end
end

% (4) compute residual
% ------------------------------------------------------------------
T = P(n+1,1);
for k=1:n
    T = bipolymultiply(T,1,-alfa(k),-beta(k));
end

ost = P - T;
if (fixed==0) && (norm([ost(:,1) diag(rot90(ost))])>1e-12*Pnorm)
    fail = 5;
    return
end

R = ost(1:n-1,2:n);
if (fixed==0) && (norm(R)>1e4*norm(Pnorm))
    fail = 6;
    return
end

% (5) n = 2
% ------------------------------------------------------------------
if n==2
    B1 = [P(3,1) 0; 0 1];
    C1 = [-P(3,1)*alfa(1) -R(1,1); 0 -alfa(2)];
    A1 = [-P(3,1)*beta(1) 0; 1 -beta(2)];
end

% (6) n = 3
% ------------------------------------------------------------------
if n==3
    lead = P(4,1)^(1/3);
    B1 = [lead 0 R(2,1); 0 lead 0; 0 0 lead];
    C1 = [-lead*alfa(1) 0 R(1,2); 1 -lead*alfa(2) 0; 0 0 -lead*alfa(3)];
    A1 = [-lead*beta(1) 0 R(1,1); 0 -lead*beta(2) 0; 0 1 -lead*beta(3)];
end

% (7) n = 4
% ------------------------------------------------------------------
if n==4
    % (7a) substitution
    R1 = bipolysubstitution(R,[1 alfa(3) beta(3); 0 1 0; 0 0 1]);
    % (7b) apply Lemma 4.1 to R1
    [AT,BT,CT] = bipoly_detrep_quadric(R1);
    % (7c) substitution back
    AS = AT - beta(3)*BT;
    BS = BT;
    CS = CT - alfa(3)*BT;
    % ansatz
    lead = P(5,1)^(1/4);
    nor1 = norm([AS(1,1) BS(1,1) CS(1,1)]);
    nor2 = norm([AS(2,2) BS(2,2) CS(2,2)]);
    equa = sqrt(nor1*nor2);
    if equa==0
        nor1 = 1;
        nor2 = 1;
    end
    B1 = [lead 0 0 0; 0 lead BS(1,1)*equa/nor1 BS(1,2)/lead; 0 0 lead BS(2,2)*equa/nor2; 0 0 0 lead];
    C1 = [-lead*alfa(1) -1 0 0; 0 -lead*alfa(2) CS(1,1)*equa/nor1 CS(1,2)/lead; 0 0 -lead*alfa(3) CS(2,2)*equa/nor2; 0 0 0 -lead*alfa(4)];
    A1 = [-lead*beta(1)  0 0 0; 0 -lead*beta(2) AS(1,1)*equa/nor1 AS(1,2)/lead; 0 0 -lead*beta(3) AS(2,2)*equa/nor2; 1 0 0 -lead*beta(4)];
end

% (8) n = 5
% ------------------------------------------------------------------
if n==5
    % (8a) change of variables
    R1 = bipolysubstitution(R,inv(Q5));
    % (8b) recursive call
    [AT,BT,CT] = bipoly_detrep_quintic_aux(R1,1);
    % (8c) substitution back
    BS = Q5(1,1)*BT + Q5(2,1)*CT + Q5(3,1)*AT;
    CS = Q5(1,2)*BT + Q5(2,2)*CT + Q5(3,2)*AT;
    AS = Q5(1,3)*BT + Q5(2,3)*CT + Q5(3,3)*AT;
    % ansatz
    lead = P(6,1)^(1/5);
    B1 = [lead 0 0 0 0; 0 lead BS(1,1) 0 BS(1,3)/lead^2; 0 0 lead BS(2,2) 0; ...
           0 0 0 lead BS(3,3); 0 0 0 0 lead];
    C1 = [-lead*alfa(1) 1 0 0 0; 0 -lead*alfa(2) CS(1,1) 0 CS(1,3)/lead^2; ...
           0 0 -lead*alfa(3) CS(2,2) 0; 0 0 0 -lead*alfa(4) CS(3,3); 0 0 0 0 -lead*alfa(5)];
    A1 = [-lead*beta(1) 0 0 0 0; 0 -lead*beta(2) AS(1,1) 0 AS(1,3)/lead^2; ...
           0 0 -lead*beta(3) AS(2,2) 0; 0 0 0 -lead*beta(4) AS(3,3); 1 0 0 0 -lead*beta(5)];
end

% back substitution
B = B1; C = C1; A = A1; 
fail = 0;

if (fixed==0) && (max([norm(A) norm(B) norm(C)])>1e10*norm(P))
   fail = 10;   
end

%-----------------------------------------------------------------------
% Determinantal representation of a quadric from Lemma 4.1
%-----------------------------------------------------------------------

function [A,B,C] = bipoly_detrep_quadric(P)

% P is a 3x3 skew upper triangular matrix

normP = norm(P,'fro');

if max(abs(P(1,1)),abs(P(1,3))) < 1e-12*normP
    % we neglect coefficients P(1,1) and P(1,3)
    A = [0 -P(2,1); 0 P(1,2)];
    B = [0 -P(3,1); 1 P(2,2)];
    C = [1 0; 0 0];
elseif abs(P(1,1))>abs(P(1,3))   
    pyx = [P(1,1) P(2,1) P(3,1)];
    alfa = roots(pyx);
    pyz = [P(1,1) P(1,2) P(1,3)];
    beta = roots(pyz);
    T = P(1,1);
    T = bipolymultiply(T,-alfa(1),-beta(1),1);
    T = bipolymultiply(T,-alfa(2),-beta(2),1);
    R = P - T;
    A = [P(1,1) 0;0 1];
    B = [-P(1,1)*alfa(1) 0; 1 -alfa(2)];
    C = [-P(1,1)*beta(1) -R(2,2);0 -beta(2)];
else
    if normP==0
        A = [0 0;0 0];
        B = [0 0;1 0];
        C = [0 0;1 0];
    else
        pyx = [P(1,3) P(2,2) P(3,1)];
        alfa = roots(pyx);
        pyz = [P(1,3) P(1,2) P(1,1)];
        beta = roots(pyz);
        T = P(1,3);
        T = bipolymultiply(T,-alfa(1),1,-beta(1));
        T = bipolymultiply(T,-alfa(2),1,-beta(2));
        R = P - T;
        A = [-P(1,3)*beta(1) -R(2,1);0 -beta(2)];
        B = [-P(1,3)*alfa(1) 0; 1 -alfa(2)];
        C = [P(1,3) 0;0 1];
    end
end

%-----------------------------------------------------------------------
% Auxiliary function that returns indices of the two most distinct 
% elements in a given vector
%-----------------------------------------------------------------------

function [i,j] = max_dist(x)

x = x(:);
n = length(x);
M = abs(x*ones(1,n)-ones(n,1)*x.');
[tmp, indi] = max(M);
[tilda, j] = max(tmp);
i = indi(j);