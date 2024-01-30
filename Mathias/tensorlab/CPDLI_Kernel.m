classdef CPDLI_Kernel < handle & matlab.mixin.CustomDisplay
%CPDLI_KERNEL Computational routines for CPD.
%   CPDLI_KERNEL(T, X, Z) is a class collecting computational routines for the
%   canonical polyadic decomposition of an Nth-order incomplete tensor T into
%   factors U{1}, ..., U{N} with linear constraints, i.e., U{n} = X{n}*Z{n} in
%   which X{n} are known matrices and Z{n} are variables, as discussed in
%   [1]. Z is a cell with initial values.
%
%   CPDLI_Kernel is an internal class. cpdli_nls and cpdli_minf are
%   interfaces to the user. Using CPDLI_Kernel is only recommended for
%   advanced users who want to have fine-grained control over all settings
%   and functions. To access more advanced information, type
%
%       edit CPDLI_Kernel
%
%   See also: cpdli_nls, cpdli_minf
%
%   References: 
%   [1] Vervliet, N., Debals, O., De Lathauwer, L., "Canonical polyadic
%       decomposition of incomplete tensors with linearly constrained factors",
%       Technical Report 16-172, ESAT-STADIUS, KU Leuven, Leuven, Belgium, 2017.

%   CPDLI_Kernel defines the following functions. Let kernel be an instance of
%   the CPDLI_KERNEL class, constructed for a tensor T, known matrices X and
%   unknown coefficients Z (an initial guess suffices):
%
%       kernel = CPDLI_Kernel(T, X, Z)
%
%   The following objective function routines can be used:
%    
%       kernel.objfun(z)
%
%   The following gradient functions can be used:
%  
%       kernel.grad(z)
%
%   The following Gramian functions can be used:
%
%       kernel.JHJ(z)
%
%   The following Gramian-vector product functions can be used:
%
%       kernel.JHJx(z, x)      
%
%   The following preconditioners can be used:
%
%       kernel.M_blockJacobi(~, b)  
%
%   The initialize and state functions cache computations which have to be done
%   once or once per iteration, respectively:
%
%       kernel.initialize()
%       kernel.state(Z)
%
%   The functions above have data-dependent and data-independent
%   implementations. The data-dependent is selected by setting the following
%   options to true:
%      
%       useDataObjFun        % use data-dependent objective function 
%       useDataGrad          % use data-dependent gradient 
%       useDataJHJ           % use data-dependent Gramian
%       useDataJHJx          % use data-dependent Gramian-vector product
%       useSparseGrad        % use sparse mtkrprod in gradient if true 
%
%   For large-scale problems, the Gramian-vector products are used, and a
%   preconditioner is used to improve convergence of the CG algorithm. For
%   reproducibility reasons in [1], the following options are available:
%
%       usePreconditioner    % Disable preconditioner if false
%       useOldPreconditioner % Use non-stochastic preconditioner 
%
%   As explained in [1], a solution with a machine precision accuracy may not
%   be attainable using the data-independent implementation. When this
%   accuracy loss is detected, the algorithm switches to the data-dependent
%   implementation automatically. This can be disabled by setting the
%   following option to false:
%
%       allowChangeToDependent 
%
%   See doc CPDLI_Kernel for more advanced parameters.
%
%   To use CPDLI_Kernel, a kernel is initialized, options are set, state is
%   called and an optimization algorithm is run. For example:
%
%       % Construct kernel
%       kernel = CPDLI_Kernel(T,X,R);
%       % Set options
%       kernel.useDataObjFun = true 
%       kernel.useDataJHJx   = true;
%       kernel.useDataGrad   = true;
%       kernel.useSparseGrad = true;
%       kernel.useDataObjFun = true;
%       % Initialize kernel
%       kernel.initialize(Z0);
%       % Call optimization routine
%       fval    = @(z) kernel.objfun(z);
%       dF.JHF  = @(z) kernel.grad(z);
%       dF.JHJx = @(z,x) kernel.JHJx(z,x);
%       dF.M    = @(z,b) kernel.M_blockJacobi(z,b);
%       [Zres, output] = nls_gndl(fval, dF, Z0);
    
% Author(s):  Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%             Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2017/04/02   NV      Initial documented and tested version for all routines.

    properties (SetAccess = protected)
        % Properties of the CPDLI_Kernel that cannot be changed after
        % construction, i.e., are read-only. A new kernel should be created
        % for other tensors, ranks or basis matrices. 
        
        T;                    % data
        X;                    % known (basis) matrices
        R;                    % rank
        N;                    % order
        size_tens;            % size of data
        size_x;               % size of variables
    end 
    
    properties
        
        %% Behavioral settings for CPDLI. 
        % Data-dependent or data-independent versions of the various computational
        % routines can be selected, as well as some implementation specific
        % properties.
        
        useDataObjFun         % Use data dependent implementation for objective function
        useDataGrad           % Use data dependent implementation for gradient
        useGramian            % Use Gramian information
        useDataJHJ            % Use data dependent implementation for Gramian
        useDataJHJx           % Use data dependent implementation for
                              % Gramian-vector product
        useSparseGrad         % Use sparse gradient implementation
        usePreconditioner     % Use preconditioner
        
        %% Advanced settings
        
        % Use non-stochastic preconditioner (without f)
        useOldPreconditioner  
        % In JHJx, first multiply then extend if true. If false, first
        % extend, then multiply.
        useJHJxMultiplyFirst  
        % Automatically change to data-dependent if needed. (If true.)
        allowChangeToDependent; 
        % Use dense tensor instead sparse tensor in mtkronprod computations,
        % if true.
        useSparseToFull       
        
        %% Cached variables
        
        residual;          % residual
        XE;                % extended basis matrices
        J;                 % (partial) Jacobian 
        W;                 % inner product of row-wise Khatri-Rao product of
                           % extended basis matrix
        XHX;               % inner product of basis matrices
        invUHU;            % (pseudo) inverse of inner products
        invXHX;            % (pseudo) inverse of inner products
        TX;                % projected tensor on basis
        T2;                % squared norm of tensor
        offset;            % offsets for vectorized variables
        fval0;             % first function value
        isInitialized      % True if initialize has been run
        
        %% Regularization 
        
        % Regularization for each variable. Empty if no regularization is required.
        regL2;            
    end 
    
    methods
        
        function this = CPDLI_Kernel(T, X, R, varargin)
        %CPDLI_Kernel Construct new kernel for CPDLI.
        %   KERNEL = CPDLI_KERNEL(T,X,R) constructs a new CPDLI kernel containing
        %   the computational routines necessary to decompose an incomplete Nth
        %   order tensor T using a rank-R CPD with factor matrices U{1}, U{2},
        %   ..., U{N}. The factor matrices have linear constraints, i.e., U{n} =
        %   X{n} * Z{n} with Z{n} the optimization variables and X{n} known
        %   matrices.
        %   
        %   KERNEL = CPDLI_KERNEL(T,X,Zinit) uses the initial factor matrices
        %   to derive the chosen rank R. 
        %   
        %   CPDLI_Kernel(T,X,R,options) or CPDLI_Kernel(T,X,R,'key',value) can be
        %   used to set the properties of the kernel:
        %     - UseDataObjFun
        %     - UseDataGrad
        %     - UseSparseGrad
        %     - UseDataJHJ
        %     - UseDataJHJx
        %     - UseGramian
        %     - UsePreconditioner
        %     - UseOldPreconditioner
        %     - AllowChangeToDependent
        %     - UseJHJxMultiplyFirst
        %     - RegL2 
        %   See doc CPDLI_Kernel for more info and properties. 
        %
        %   See also: initialize

            % Check and store T
            if ~isstruct(T)
                if ~isnumeric(T)
                    error('CPDLI_Kernel:T', 'T should be an incomplete tensor');
                else 
                    T = fmt(T); 
                end 
            end
            if ~isstruct(T) || ~isfield(T, 'incomplete') || ~T.incomplete
                error('CPDLI_Kernel:T', 'T should be an incomplete tensor');
            end 
            this.T = T;
            this.N = getorder(T);
            this.size_tens = getsize(T);
            
            % Check and store X
            X = X(:).';
            if ~iscell(X) || numel(X) ~= this.N || any(~cellfun(@ismatrix, X))
                error('CPDLI_Kernel:X', ...
                      'X should be a cell with numel(X) =  getorder(T) matrices.');                                
            end
            if any(cellfun('size', X, 1) ~= this.size_tens)
                error('CPDLI_Kernel:X', ...
                      'size(X{n},1) should be getsize(T,n) for all n.');                
            end 
            this.X = X;
            this.size_x = cellfun('size', X, 2);
            
            % Check and store R
            if iscell(R), R = size(R{1}, 2); end
            if ~isscalar(R)
                error('CPDLI_Kernel:R', ...
                      'R should be a scalar or a cell of coefficient matrices (Z).');
            end
            this.R = R;
            this.fval0 = inf;
            
            p = inputParser;
            p.addOptional('UseDataObjFun', true);
            p.addOptional('UseDataGrad', true);
            p.addOptional('UseSparseGrad', true);
            p.addOptional('UseGramian', true);
            p.addOptional('UseDataJHJ', true);
            p.addOptional('UseDataJHJx', true);
            p.addOptional('UsePreconditioner', true);
            p.addOptional('UseOldPreconditioner', false);
            p.addOptional('AllowChangeToDependent', true);
            p.addOptional('UseJHJxMultiplyFirst', false);
            p.addOptional('UseSparseToFull', nan);
            p.addOptional('RegL2', []);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            for f = fieldnames(p.Results)'
                fname = [lower(f{1}(1)) f{1}(2:end)];
                this.(fname) = p.Results.(f{1});
            end 
            
            % Check and reformat regularization terms to long vector
            this.convertRegL2();
                
            this.isInitialized = false;
            
            if ~islogical(p.Results.UseSparseToFull)
                scale = length(this.T.val)/prod(this.size_tens);
                this.useSparseToFull = prod(this.size_tens) < 1e9/8 && scale > 0.05;
            end 
        end 
        
        function initialize(this, ~)
        %INITIALIZE Initialize kernel before computations. 
        %  INITIALIZE initializes the kernel before iterations start and should be
        %  called before calling the optimization routine. Results indepedent of
        %  the variables are precomputed here. INITIALIZE should be called every
        %  time a setting is altered to ensure a correct state.
            
            N = this.N;
            R = this.R;
            
            % Compute data-independent matrix W if needed needed
            if ~this.useDataObjFun || ~this.useDataGrad || ...
                    (this.useGramian && (~this.useDataJHJ || ~this.useDataJHJx)) 
                this.W = ipkrt(this.X(end:-1:1), 'subs', this.T.sub(end:-1:1));
            end 
            % Compute extended known matrices XE if needed. 
            if this.useGramian && (this.useDataJHJ || this.useDataJHJx)
                this.XE = cellfun(@(x,i) x(i,:), this.X, this.T.sub, ...
                                  'UniformOutput', false);
                this.J = cell(1, N);
            end 
            % Precompute projected data TX if needed
            if ~this.useDataGrad || ~this.useDataObjFun
                this.T.incomplete = false;
                this.T.sparse     = true;
                if this.useSparseToFull
                    this.TX = tmprod(ful(this.T), this.X, 1:N, 'H');
                else 
                    tmp = mtkronprod(this.T, this.X, 0);
                    this.TX = reshape(full(tmp), this.size_x);
                end 
                this.T.incomplete = true;
                this.T.sparse     = false;
            end 
            % Precompute norm of data tensor if needed
            if ~this.useDataObjFun
                this.T2 = this.T.val'*this.T.val;
            end
            
            % Precompute inverses of known matrices for preconditioner if necessary
            if this.useGramian && this.usePreconditioner
                % For reproducible research reasons, the old non-stochastic
                % preconditioner can be selected, which ignores the distribution
                % of the known entries across the tensor. 
                if this.useOldPreconditioner
                    freq = repmat({1},1,N); 
                else
                    % Compute distribution of known entries. 
                    freq = cell(1,N);
                    for n = 1:N
                        ntot = prod(this.size_tens([1:n-1 n+1:N]));
                        freq{n} = histc(this.T.sub{n}, 1:this.size_tens(n))/ntot;
                    end 
                end 
                
                this.XHX = cellfun(@(x,f) x'*diag(f)*x, this.X, freq, ...
                                   'UniformOutput', false);
                this.invXHX = cellfun(@pinv, this.XHX, 'UniformOutput', false);
            end
            
            % Precompute offsets
            nbvar = this.size_x(:).' * R;
            this.offset = cumsum([1 nbvar]);
            
            % Reformat regularization terms to long vector
            this.convertRegL2();

            this.isInitialized = true;
        end 
        
        function state(this, z)
        %STATE Update the state of the kernel. 
        %  STATE(Z) updates the state for a regular iteration of the algorithm,
        %  which includes caching results that change every iteration but are
        %  use by multiple funtions.
            
            N = this.N;
            R = this.R;

            if this.useGramian && this.usePreconditioner
                % Cache the factor matrices' Gramians.
                UHU = zeros(N,R*R);
                for n = 1:N
                    tmp = conj(z{n}'*this.XHX{n}*z{n});
                    UHU(n,:) = tmp(:);
                end
            
                % Cache the inverses of the Gramians for the preconditioner.
                this.invUHU = cell(1,N);
                for n = 1:N
                    tmp = UHU([1:n-1 n+1:N],:);
                    if N > 2, tmp = prod(tmp,1); end
                    this.invUHU{n} = inv(reshape(tmp,[R R]));
                end
            end 
            
            % Precompute (partial) Jacobian matrices. 
            if this.useGramian && (this.useDataJHJx || this.useDataJHJ)
                size_x = this.size_x;
                xz = cellfun(@(x,z) x*z, this.X, z, 'uni', 0);
                for n = 1:N
                    dims = [1:n-1 n+1:N];
                    val = this.XE{dims(1)}*z{dims(1)};
                    for k = dims(2:end)
                        val = val .* (this.XE{k}*z{k});
                    end
                    this.J{n} = val;
                end
            end 
            
            % For relatively small tensors, treating a sparse residual as a
            % dense residual is more efficient when computing the gradient.
            if this.useDataGrad && isstruct(this.residual) && this.useSparseToFull
                this.residual.sparse = 1;
                this.residual.incomplete = 0;
                this.residual = ful(this.residual);
            end 
        end

        
        function fval = objfun(this, z)
        %OBJFUN Objective function value for CPDLI. 
            
            if ~isempty(this.regL2)
                % Add regularization terms
                zvec = cellfun(@(z) z(:), z, 'UniformOutput', false);
                zvec = vertcat(z{:});
                fval = 0.5 * this.regL2(:).' * (abs(zvec(:)).^2);
            else 
                fval = 0;
            end 
            if this.useDataObjFun
                xz = cellfun(@(x,z) x*z, this.X, z, 'UniformOutput', false);
                this.residual = cpdres(this.T, xz, 'Format', false);
                fval = fval + 0.5*(this.residual.val'*this.residual.val);
            else 
                Tz = cpdgen(z);
                fval = fval + real((0.5*Tz(:)'*this.W - this.TX(:)')*Tz(:)+0.5*this.T2);
                
                % Change to data-dependent implementation if allowed and needed.
                if this.allowChangeToDependent && ...
                        (isfinite(this.fval0) && fval <= this.fval0*1e-8)
                    this.useDataObjFun = true;
                    fval = this.objfun(z);
                end 
                % Store first objective value 
                if ~isfinite(this.fval0), this.fval0 = fval; end
            end 
        end 
             
        function grad = grad(this, z)
        %GRAD Gradient for CPDLI. 
            
            this.state(z);
            offset = this.offset;
            N = this.N;
            
            if ~isempty(this.regL2)
                % Gradient for regularization term
                zvec = cellfun(@(z) z(:), z, 'UniformOutput', false);
                grad = this.regL2 .* vertcat(zvec{:});
            else 
                grad = zeros(offset(end)-1,1);
            end 

            % Compute residual depending on chosen implementation
            if this.useDataGrad && ~this.useSparseGrad
                E = reshape(this.X{1}'*mtkronprod(this.residual, this.X, 1), ...
                            this.size_x);
            elseif this.useDataGrad && this.useSparseGrad
                E = this.residual;
                z = cellfun(@(u,v) u*v, this.X, z, 'UniformOutput', false);
            else 
                Tz = cpdgen(z);
                E = reshape(this.W*Tz(:), this.size_x) - this.TX;                
            end 
            % Remove matrix to avoid overhead in mtkrprod
            if isstruct(E), E.matrix = []; end
            % Compute gradient for CPD term
            for n = 1:N
                tmp = full(mtkrprod(E, z, n));
                if this.useDataGrad && this.useSparseGrad, 
                    tmp = this.X{n}'*tmp; 
                end 
                idx = offset(n):offset(n+1)-1;
                grad(idx) = grad(idx) + tmp(:);
            end 
        end
        
        function JHJ = JHJ(this, z)
        %JHJ Gramian for CPDLI.
            offset = this.offset;
            N = this.N;
            R = this.R;
            size_x = this.size_x;
            
            if ~isempty(this.regL2)
                % Gramian for regularization term
                JHJ = diag(this.regL2);
            else 
                JHJ = zeros(sum(size_x)*R);
            end 
            if this.useDataJHJ % Data-dependent implementation 
                for n = 1:N
                    idxn = offset(n):offset(n+1)-1;
                    JHJ(idxn,idxn) = ipkrt(this.J{n},this.XE{n});
                    for m = n+1:N
                        idxm = offset(m):offset(m+1)-1;
                        tmp = ipkrt2({this.J{m},this.XE{m}}, {this.J{n},this.XE{n}});
                        JHJ(idxm, idxn) = tmp;
                        JHJ(idxn, idxm) = tmp';
                    end               
                end
            else % Data-independent implementation
                W = reshape(this.W, [size_x, size_x]);
                for n = 1:N
                    idxn = offset(n):offset(n+1)-1;
                    JHJ(idxn,idxn) = reshape(krmtkrprod(W, z, z, n, n), ...
                                             [size_x(n)*R, size_x(n)*R]);
                    for m = n+1:N
                        idxm = offset(m):offset(m+1)-1;
                        tmp = reshape(krmtkrprod(W, z, z, m, n), ...
                                      [size_x(m)*R, size_x(n)*R]);
                        JHJ(idxm, idxn) = tmp;
                        JHJ(idxn, idxm) = tmp';
                    end               
                end
            end 
        end
        
        function y = JHJx(this, z, x)
        %JHJX Gramian-vector product for CPDLI.
            
            N  = this.N;
            R  = this.R;
            size_x = this.size_x;
            
            if this.useDataJHJx 
                % Data-dependent implementation
                nke = length(this.T.val);
                J = this.J;
                Jx = zeros(nke,1);
                offset = 0;
                for n = 1:N
                    xn = reshape(x(offset + (1:R*size_x(n))), size_x(n), R);
                    if this.useJHJxMultiplyFirst
                        tmp = this.X{n} * xn;
                        tmp = tmp(this.T.sub{n},:);
                    else 
                        tmp = this.XE{n} * xn;
                    end 
                    Jx = Jx + sum(J{n}.*tmp, 2);
                    offset = offset + R*size_x(n);
                end 
                offset = 0;
                y = zeros(size(x));
                for n = 1:N
                    idxn = offset + (1:R*size_x(n));
                    % Sparse matrix version of conj(bsxfun(@times, conj(this.XE{n}), Jx)'*J{n})
                    tmp = sparse(double(this.T.sub{n}), 1:nke, Jx, this.T.size(n), ...
                                 nke);
                    tmp = tmp * conj(J{n});
                    tmp = this.X{n}'*tmp;
                    y(idxn) = tmp(:);
                    offset = offset + R*size_x(n);
                end 
            else 
                % Data-independent implementation
                % Step 1: compute Jx = W*J*x
                Jx = zeros(prod(size_x),1);
                offset = 0;
                for n = 1:N
                    ztmp = z;
                    ztmp{n} = reshape(x(offset + (1:R*size_x(n))), size_x(n), R);
                    tmp = cpdgen(ztmp);
                    Jx = Jx + tmp(:);
                    offset = offset + R*size_x(n);
                end 
                Jx = this.W*Jx;
                Jx = reshape(Jx, size_x);
                % Step 2: compute y = J'*Jx
                offset = 0;
                y = zeros(size(x));
                for n = 1:N
                    idxn = offset + (1:R*size_x(n));
                    tmp = mtkrprod(Jx, z, n);
                    y(idxn) = tmp(:);
                    offset = offset + R*size_x(n);
                end 
            end 
            if ~isempty(this.regL2)
                y = y + this.regL2.*x;
            end 
        end

        function x = M_blockJacobi(this,~,b)
        %M_BLOCKJACOBI Stochastic Block-Jacobi preconditioner for CPDLI.
            N = this.N;
            R = this.R;
            offset = this.offset;

            x = nan(size(b));
            for n = 1:length(offset)-1
                idx = offset(n):offset(n+1)-1;
                B = reshape(b(idx),[],size(this.invUHU{n},1));
                tmp = this.invXHX{n} * B * this.invUHU{n};
                x(idx) = tmp(:);
            end
        end


        function isvalid = validate(this, z)
        %VALIDATE Check whether this kernel is in a valid.
        %  KERNEL.VALIDATE(Z) checks whether the combination of data and
        %  initial values Z is valid. This is a deep check, and may take some
        %  time, but it ensures that the computation will succeed.
            
            if ~isvalidtensor(this.T, false, 'incomplete');
                error('CPDLI_Kernel:T', ['T is not a valid incomplete tensor. ' ...
                                    'See isvalidtensor.']);
            end 
            if numel(this.X) ~= this.N 
                error('CPDLI_Kernel:X', ...
                      'numel(X) should be getorder(T).');                
            end 
            if numel(z) ~= this.N 
                error('CPDLI_Kernel:z', ...
                      'numel(Z) should be getorder(T).');                
            end                 
            if any(cellfun('size', this.X, 1) ~= this.size_tens)
                error('CPDLI_Kernel:X', ...
                      'size(X{n},1) should be getsize(T,n) for all n.');                
            end 
            if any(this.size_x ~= cellfun('size', z, 1))
                error('CPDLI_Kernel:z', ...
                      'size(z{n},1) should be size(X{n},n) for all n.');     
            end 
            if any(this.R ~= cellfun('size', z, 2))
                error('CPDLI_Kernel:z', ...
                      'size(z{n},2) should be R for all n.');     
            end 
            
            % Convert regL2 again in the case it is changed after construction 
            this.convertRegL2();  
            % Results should be 
            nbvar = sum(cellfun(@numel, z));
            if ~isempty(this.regL2) && numel(this.regL2) ~= nbvar
                error('CPDLI_Kernel:regL2', ['The number of regularization ' ...
                                    'parameters is incorrect.']);
            end 
            
            isvalid = true;
        end 
    end        
    
    methods (Access = private)
        
        function convertRegL2(this)
        %CONVERTREGL2 Convert different regL2 formats to standard.
            
            N      = this.N;
            R      = this.R;
            size_x = this.size_x;
            nbvar  = size_x(:).' * R;
            
            if iscell(this.regL2)
                if isempty(this.regL2) || all(cellfun(@isempty, this.regL2))
                    this.regL2 = [];
                elseif numel(this.regL2) ~= N
                    error('CPDLI_Kernel:regL2', ['The number of regularization ' ...
                                        'parameters is incorrect.']);
                else 
                    for n = 1:N
                        if isempty(this.regL2{n})
                            this.regL2{n} = zeros(size_x(n),R);
                        elseif numel(this.regL2{n}) == 1
                            this.regL2{n} = ones(size_x(n),R)* ...
                                this.regL2{n};
                        elseif numel(this.regL2{n}) == size_x(n)
                            this.regL2{n} = this.regL2{n}(:) * ones(1,R);
                        end 
                    end
                    this.regL2 = cellfun(@(x) x(:), this.regL2, ...
                                         'UniformOutput', false);
                    this.regL2 = vertcat(this.regL2{:});
                end 
            end 
            if numel(this.regL2) == 1
                this.regL2 = ones(sum(nbvar),1) * this.regL2;
            elseif numel(this.regL2) == N
                tmp = arrayfun(@(lambda,n) ones(n,1) * lambda, this.regL2(:), nbvar(:), ...
                               'UniformOutput', false);
                this.regL2 = vertcat(tmp{:});
            elseif numel(this.regL2) ~= sum(nbvar) && numel(this.regL2) > 0
                error('CPDLI_Kernel:regL2', ['The number of regularization ' ...
                                    'parameters is incorrect.']);
            end 
            this.regL2 = this.regL2(:);
        end 
    
    end 
    
    methods (Access = protected)
                
        function header = getHeader(this)
            if ~isscalar(this)
                header = getHeader@matlab.mixin.CustomDisplay(this);
            else
                className = matlab.mixin.CustomDisplay.getClassNameForHeader(this);
                newHeader = [className,[' CPD kernel for incomplete tensor ' ...
                                    'with linear constraints: ']];
                header = sprintf('%s\n',newHeader);
            end
        end
        
        function propgrp = getPropertyGroups(this)
            if ~isscalar(this)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(this);
            else
                dataprops = {'T', 'size_tens', 'X', 'size_x', 'R', 'regL2'};
                propgrp(1) = matlab.mixin.util.PropertyGroup(dataprops, ...
                                                             'Data properties');                
                settingprops = {'useDataObjFun', 'useDataGrad', 'useSparseGrad', ...
                               'useDataJHJ', 'useDataJHJx', 'usePreconditioner'};
                propgrp(2) = matlab.mixin.util.PropertyGroup(settingprops, ...
                                                             'Computational settings');
                specializedprops = {'useJHJxMultiplyFirst', 'useOldPreconditioner', ...
                                    'allowChangeToDependent', 'useSparseToFull'};
                propgrp(3) = matlab.mixin.util.PropertyGroup(specializedprops, ...
                                                             'Specialized settings');
            
            end
        end
        
        function footer = getFooter(this)
            if isscalar(this)
                footer = sprintf('%s\n','Cached variables hidden. ');
            else
                footer = '';
            end
        end
    end 
end 
