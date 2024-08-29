% BiRoots Toolbox 
% Version 2.0               05-Dec-2016
% 
% Toolbox for roots of systems of bivariate polynomials
% 
% Roots of a polynomial can be computed as eigenvalues of the companion 
% matrix. In this toolbox this approach is generalized to systems of
% bivariate polynomials and the roots are computed as eigenvalues of 
% appropriate two-parameter eigenvalue problems.
% 
% ------------------------------------------------------------------------
% Determinantal representation of a bivariate polynomial
% ------------------------------------------------------------------------
% 
% Square matrices A,B,C form a determinantal representation (or 
% linearization) of a polynomial p(x,y) if det(A+x*B+y*C) = p(x,y)
% 
% Function bipoly_detrep in the toolbox returns a determinantal 
% representation of a given bivariate polynomial. Several methods are 
% available that give linearizations of different sizes.
% 
% ------------------------------------------------------------------------
% Supported determinantal representations of bivariate polynomials
% ------------------------------------------------------------------------
% 
%   - bipoly_detrep_nxn : representation of size n for a bivariate 
%     polynomial of degree n from [3]. 
%   - bipoly_detrep_quintic : representation of size n for a bivariate 
%     polynomial of degree 2<=n<=5 from [2], faster than bipoly_detrep_nxn
%   - bipoly_detrep_unif : representation of size 2n-1, where no 
%     computation is involved from [4]
%   - bipoly_detrep_123 : representations from [1]. One is of uniform type 
%     of asymptotic size  n^2/4, the second one involves computation and is 
%     of asymptotic size n^2/6. These were the representations in version 
%     1.0, the new representations from [2], [3], [4] are more efficient.
% 
% ------------------------------------------------------------------------
% Two-parameter eigenvalue problem
% ------------------------------------------------------------------------
% 
% A matrix two-parameter eigenvalue problem has the form
% 
% A1*x + lambda*B1*x + mu*C1*x = 0, 
% A2*y + lambda*B2*y + mu*C2*y = 0.
% 
% We are looking for an eigenvalue (lambda,mu) and nonzero eigenvectors 
% x,y such that the above system is satisfied. 
% 
% ------------------------------------------------------------------------
% Finding roots of a systems of bivariate polynomials
% ------------------------------------------------------------------------
% 
% To find solutions for the system p1(x,y)=0, p2(x,y)=0, where p1 and p2
% are bivariate polynomials, we find determinantal representations for 
% p1 and p2 such that
% 
% det(A1 + x*B1 + y*C1) = p1(x,y),
% det(A2 + x*B2 + y*C2) = p2(x,y).
% 
% The eigenvalues of the two-parameter eigenvalue problem 
% 
%    A1*x + lambda*B1*x + mu*C1*x = 0, 
%    A2*y + lambda*B2*y + mu*C2*y = 0.
% 
% are then the roots of p1(x,y)=0, p2(x,y)=0. 
% 
% If the size of A1 (A2) is larger than the degree of p1 (p2), then the 
% above two-parameter eigenvalue problems is singular. It can be solved 
% with a staircase-type algorithm, but due to the large size of the 
% matrices from the linearization and the difficulties of detecting the 
% correct numerical rank in the staircase algorithm, this works efficiently 
% and accurately only for polynomials of small order. It is not recommended 
% to use this toolbox for polynomials of degree 13 or more. It can also
% fail for some polynomials of smaller degree, where double precision is
% not sufficient for this approach.
%
% ------------------------------------------------------------------------
% Change log from version 1.0:
%  - Faster computation using three new representations: 
%       - (2n-1) x (2n-1) uniform representation, 
%       - n x n representation for square-free polynomials, 
%       - n x n representation for polynomials of degree <=5. 
%  - Biroots contains new heuristic approaches
%  - Contents.m is changes so that BiRoots shows in the list of installed 
%    packages
% 
% ------------------------------------------------------------------------
% Main Matlab functions in the toolbox 
% ------------------------------------------------------------------------
% 
%   - biroots: roots of a system of bivariate polynomials
%   - bipoly_detrep: determinantal representation of a bivariate polynomial
%   
% ------------------------------------------------------------------------
% References
% ------------------------------------------------------------------------
% 
% [1] B. Plestenjak, M. E. Hochstenbach: Roots of bivariate polynomial 
%     systems via determinantal representations, SIAM J. Sci. Comput. 38 
%     (2016) A765-A788
% [2] A. Buckley, B. Plestenjak: Simple determinantal representations of 
%     up to quintic bivariate polynomials, arXiv.1609.00498
% [3] B. Plestenjak: Minimal determinantal representations of 
%     bivariate polynomials, arXiv.1607.03969
% [4] A. Boralevi, J. van Doornmalen, J. Draisma, M. E. Hochstenbach, 
%     B. Plestenjak: Uniform determinantal representations, 
%     arXiv.1607.04873
% 
% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt