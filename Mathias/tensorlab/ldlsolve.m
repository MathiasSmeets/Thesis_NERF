function x = ldlsolve(A,b,varargin)
%LDLSOLVE Solve linear system via LDL factorization.
%   X = LDLSOLVE(A,B) solves the linear system A*x = b for x in the case A is
%   a positive semidefinite symmetric/Hermitian matrix, using the LDL
%   factorization of A.  

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/02/19   NV      Initial version

    if ~isvector(b) && nargin >= 4
        L = A;
        D = b;
        p = varargin{1};
        b = varargin{2};
        varargin = varargin(3:end);
    else 
        % Compute LDL factorization
        [L,D,p] = ldl(A,'vector');
    end 

    
    ip = inputParser;
    ip.addOptional('tol', eps(class(A))*max(size(A,1),10));
    ip.parse(varargin{:});
    options = ip.Results;
    
 
    % Check D is diagonal and not block diagonal
    if any(diag(D,-1)~=0) || any(diag(D,1)~=0)
        if any(abs(diag(D,1)) > max(abs(diag(D)))*options.tol)
            error('A is not symmetric or not PSD.');
        end
    end 

    % pseudoinversion of diagonal of D to take zeros into account
    d = diag(D);
    dinv = 1./d;
    tol = max(abs(d))*options.tol;
    dinv(abs(d) < tol) = 0;
    dinv = spdiags(dinv, 0, numel(d), numel(d));
    
    % pseudoinversion of blocks on diagonal of D
    blkidx = find(abs(diag(D,1)) > eps(class(A)));
    for k = 1:numel(blkidx)
        idx = blkidx(k):blkidx(k)+1;
        [u,s,v] = svd(D(idx,idx));
        s = diag(s);
        s = (1./s) .* (s > tol);
        dinv(idx,idx) = v*diag(s)*u';
    end 
        
    % Solve system (with permutation)
    x = L' \ dinv * (L \ b(p,:));
    x(p,:) = x;
end
