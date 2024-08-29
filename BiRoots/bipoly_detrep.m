function [A, B, C, method, hist] = bipoly_detrep(P, method, opts)

%bipoly_detrep  Determinantal representation of a bivariate polynomial
%
% [A,B,C,meth,iter] = bipoly_detrep(P,method,opts)
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
% [A,B,C] = bipoly_detrep(P)
%
% method: 1: linearization 1 (matrices of size n^2/4) - no computation [1]
%         2: linearization 2 (matrices of size n^2/6) - computation [1]
%         3: linearization 2 with special construction for degree 3 or 4 
%            (matrices of size n^2/6) - computation [1]
%         4: bipoly_detrep_quantic - size n for degree 5 or less [2]
%         5: bipoly_detret_nxn - size n for a square-free polynomial [3]
%         6: bipoly_detrep_unif - matrices of size 2n-1 - no computation [4] 
%         0: (default) not specified, default method is 4 for degree 5 or 
%            less, 5 for degree 9 or less, and 6 otherwise 
%         
% options in opts:
%    - options available in bipoly_detrep_123, bipoly_detrep_quintic, 
%      bipoly_detrep_nxn, and bipoly_detrep_unif
% 
% Output:
%    - A,B,C: matrices of the representation
%    - meth: method that was used 
%    - iter: number of additional iterations required (0 if none) or 
%            (-1)*(number of iterations) if the method fails
%
% See also: biroots

% References: 
%   [1] B. Plestenjak, M. E. Hochstenbach: Roots of bivariate 
%       polynomial systems via determinantal representations, SIAM J. Sci. 
%       Comput. 38 (2016) A765-A788
%   [2] A. Buckley, B. Plestenjak: Simple determinantal representations of 
%       up to quintic bivariate polynomials, arXiv.1609.00498
%   [3] B. Plestenjak: Minimal determinantal representations of 
%       bivariate polynomials, arXiv.1607.03969
%   [4] A. Boralevi, J. van Doornmalen, J. Draisma, M. E. Hochstenbach, 
%       B. Plestenjak: Uniform determinantal representations, 
%       arXiv.1607.04873

% BiRoots toolbox
% B. Plestenjak, University of Ljubljana
% M. E. Hochstenbach, TU Eindhoven
% FreeBSD License, see LICENSE.txt

% BP 16.10.2016 : new linearizatons added, default choice depends on degree
% Last revision 05.12.2016 Bor Plestenjak

narginchk(1,3)

n = size(P,2) - 1;
defmethod = 0;

if nargin < 3
    opts = [];
end

if nargin < 2 || isempty(method) || (method==0)
    defmethod = 1;
    if n<=5 
        method = 4; % quintic (n x n representation for degree<=5) [2]
        if ~isfield(opts,'try_best'),  opts.try_best=0;  end
    elseif n<=9
        method = 5; % n x n representation for square-free polynomials [3]
    else
        method = 6; % uniform representation [4]
    end
end

switch method
    case {1,2,3} % representation from [1]
        [A,B,C] = bipoly_detrep_123(P,method,opts);
        hist = 0;
    case 4
        [A,B,C,hist] = bipoly_detrep_quintic(P,opts); % linearization in [2]
        if defmethod && (hist<0)
            method = 6;
            hist = -hist;
            [A,B,C] = bipoly_detrep_unif(P);
        end
    case 5
        [A,B,C,hist] = bipoly_detrep_nxn(P,opts);
        if defmethod && (hist<0)
            method = 6;
            hist = -hist;
            [A,B,C] = bipoly_detrep_unif(P);
        end
    case 6
        [A,B,C] = bipoly_detrep_unif(P);
        method = 6;
        hist = 0;
    otherwise
        error('Supported methods are 1,2,3,4,5, and 6')
end

