function kappa = cpd_cond(U,absolute,varargin)
%CPDCONDNICO short discription
%   long description
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/22   NV      Initial version
    
    % Compress the factors (condition number is invariant under orthogonal projection)
    [~,U] = cellfun(@(u) qr(u,0), U, 'UniformOutput', false);
    
    if nargin < 2, absolute = false; end
    p = inputParser();
    p.addOptional('LargeScale', sum(cellfun(@numel,U)) > 10000);
    p.parse(varargin{:});
    options = p.Results;

    N = numel(U);
    R = size(U{1},2);
    
    % Balance norms
    lnrm = cellfun(@(u) sqrt(sum(abs(u).^2,1)), U, 'UniformOutput', false);
    tnrm = nthroot(prod(vertcat(lnrm{:}),1), N);
    U    = cellfun(@(u,n) bsxfun(@times, u, tnrm./n), U, lnrm, 'UniformOutput', false);
    
    % Compute Gramian for CPD
    kernel = CPD_Kernel(U);
    kernel.initialize(U);
    kernel.state(U);
    G = kernel.JHDJ(U);
    
    if ~options.LargeScale
        sv = sqrt(abs(eig(G)));
        kappa = 1./sv((N-1)*R+1);
    else
        % Compute largest singular value of Jacobian (=sqrt of largest eig of Gramian)
        [L,D,p] = ldl(G, 'vector');
        afun = @(x) L*(D*(L'*x));
        svl  = sqrt(abs(eigs(afun, sum(cellfun(@numel,U)), 1, 'LM')));
        
        % Compute smallest nonzero singular value of Jacobian (= sqrt of smallest nonzero eig of Gramian)
        afun = @(x) ldlsolve(L,D,p,x,'Tol', 1e-16);
        svs  = sqrt(abs(eigs(afun, sum(cellfun(@numel,U)), (N-1)*R+1, 'SM')));
        
        kappa = 1./svs(end);
    end 
    if nargin == 1 && absolute == false
        kappa = kappa * frob(U) / sqrt(N) / norm(tnrm);
    end 
end
