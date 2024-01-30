classdef CPDI_Kernel < TensorOptimizationKernel
%CPDI_KERNEL Computational routines for CPD.
%   CPDI_KERNEL(T, R) is a class collecting computational routines for the
%   rank-R canonical polyadic decomposition of an Nth-order incomplete tensor T
%   into factors U{1}, ..., U{N}, as discussed in [1]. U is a cell with factor
%   matrices.
%
%   CPDI_Kernel is an internal class. cpdi_nls and cpdi_minf are interfaces to
%   the user. Using CPDI_Kernel is only recommended for advanced users who want
%   to have fine-grained control over all settings and functions. To access more
%   advanced information, type
%
%       edit CPDI_Kernel
%
%   See also: cpdi_nls, cpdi_minf
%
%   References: 
%   [1] Vervliet, N., Debals, O., De Lathauwer, L., "Canonical polyadic
%       decomposition of incomplete tensors with linearly constrained factors",
%       Technical Report 16-172, ESAT-STADIUS, KU Leuven, Leuven, Belgium, 2017.

%   CPDI_Kernel defines the following functions. Let kernel be an instance of
%   the CPDI_KERNEL class, constructed for a tensor T and target rank-R (or
%   alternatively, initial guesses for the factor matrices U):
%
%       kernel = CPDI_Kernel(T, R) % or 
%       kernel = CPDI_Kernel(T, U) 
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
%   For large-scale problems, the Gramian-vector products are used, and a
%   preconditioner is used to improve convergence of the CG algorithm. For
%   reproducibility reasons in [1], the following options are available:
%
%       usePreconditioner    % Disable preconditioner if false
%       useOldPreconditioner % Use non-stochastic preconditioner 
%
%   See doc CPDLI_Kernel for more advanced parameters.
%
%   To use CPDI_Kernel directly, a kernel is initialized, options are set,
%   initialize is called and an optimization algorithm is run. For example:
%
%       % Construct kernel with options 
%       kernel = CPDI_Kernel(T, R, 'UseMex', true);
%       % Initialize kernel
%       kernel.initialize(U0);
%       % Call optimization routine
%       fval    = @(z) kernel.objfun(z);
%       dF.JHF  = @(z) kernel.grad(z);
%       dF.JHJx = @(z,x) kernel.JHJx(z,x);
%       dF.M    = @(z,b) kernel.M_blockJacobi(z,b);
%       [Ures, output] = nls_gndl(fval, dF, U0);
    
% Author(s):  Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%             Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2017/04/02   NV      Initial documented and tested version for all routines.

    properties (SetAccess = protected)
        % Properties of the CPDI_Kernel that cannot be changed after construction, i.e.,
        % are read-only. A new kernel should be created for other tensors or
        % ranks.
        
        T;                    % data
        R;                    % rank
        N;                    % order of data
        size_tens;            % size of data
    end 
    
    properties
        
        %% Behavioral settings for CPDI. 
        % Data-dependent or data-independent versions of the various computational
        % routines can be selected, as well as some implementation specific
        % properties.
        
        useGramian            % Use Gramian information
        usePreconditioner     % Use preconditioner
        
        %% Advanced settings
        
        % Use non-stochastic preconditioner (without f)
        useOldPreconditioner  
        % Use dense tensor instead sparse tensor in objfun and grad
        useSparseToFull
        % Use C/mex implementation (if compiled)
        useMex;
        % Has complex factor matrices (needed when using Mex)
        useComplexMex;
        
        %% Cached variables
        
        residual;          % residual
        J;                 % Jacobian 
        invUHU;            % (pseudo) inverse of inner products
        freq;              % Distribution of known entries accross T
        offset;            % offsets for vectorized variables
        isInitialized      % True if initialize has been run    
    end 
    
    methods
        
        function this = CPDI_Kernel(T, varargin)
        %CPDI_Kernel Construct new kernel for CPDI.
        %   KERNEL = CPDI_KERNEL(T,R) constructs a new CPDI kernel containing
        %   the computational routines necessary to decompose an incomplete Nth
        %   order tensor T using a rank-R CPD with factor matrices U{1}, U{2},
        %   ..., U{N}. 
        %   
        %   KERNEL = CPDI_KERNEL(T,Uinit) uses the initial factor matrices
        %   to derive the chosen rank R. 
        %   
        %   CPDI_Kernel(T,R,options) or CPDI_Kernel(T,R,'key',value) can be
        %   used to set the properties of the kernel:
        %     - UseGramian
        %     - UsePreconditioner
        %     - UseOldPreconditioner
        %     - UseMex
        %   See doc CPDI_Kernel for more info and properties. 
        %
        %   See also: initialize

            % Check and store T
            if ~isstruct(T)
                if ~isnumeric(T)
                    error('CPDI_Kernel:T', 'T should be an incomplete tensor');
                else 
                    T = fmt(T); 
                end 
            end
            if ~isstruct(T) || ~isfield(T, 'incomplete') || ~T.incomplete
                error('CPDI_Kernel:T', 'T should be an incomplete tensor');
            end 
            this.T = T;
            this.N = getorder(T);
            this.size_tens = getsize(T);
                       
            p = inputParser;
            p.addOptional('UseGramian', true);
            p.addOptional('UsePreconditioner', true);
            p.addOptional('UseOldPreconditioner', false);
            p.addOptional('TolLargeScale', 0.05);
            p.addOptional('UseMex', true);
            p.addOptional('UseComplexFactors', false);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            options = p.Results;
            
            this.useGramian           = options.UseGramian;
            this.usePreconditioner    = options.UsePreconditioner;
            this.useOldPreconditioner = options.UseOldPreconditioner;
            this.useMex               = options.UseMex;
            this.useMex = this.useMex && exist('cpd_jacobian_incomplete.mexa64', 'file');
            this.useComplexMex        = options.UseComplexFactors;
            
            this.useSparseToFull = length(T.ind)/prod(T.size) > ...
                options.TolLargeScale;
                
            this.isInitialized = false;
        end 
        
        function initialize(this, z)
        %INITIALIZE Initialize kernel before computations. 
        %  INITIALIZE(U) initializes the kernel before iterations start and should be
        %  called before calling the optimization routine. Results indepedent of
        %  the variables are precomputed here. U is a cell containing initial
        %  factor matrices and is used to derive properties of the data, e.g.
        %  whether the factors are complex or have higher order than the tensor.
        %  INITIALIZE should be called every time a setting is altered to ensure
        %  a correct state.
            
            R = size(z{1},2);
            N = this.N;
            T = this.T;
            size_tens = this.size_tens;
            this.R = R;
            this.offset = cumsum([1 cellfun(@numel,z)]);

            if this.useMex
                % Preallocation of sparse Jacobian when using mex
                this.useComplexMex = any(~cellfun(@isreal,z)) || any(~isreal(T.val));
                if this.useComplexMex, val = 1+1i; 
                else val = 1; end

                J = cell(1,length(z));
                i = repmat((1:numel(T.val)).', 1, R);
                for n = 1:N
                    j = bsxfun(@plus, double(T.sub{n}), size_tens(n)* ...
                               ones(length(T.val),1)*(0:R-1));
                    
                    J{n} = sparse(j.',i.',val);
                end
                if N < length(z)
                    j = bsxfun(@plus, ones(size(T.sub{n})), ...
                               ones(length(T.val),1)*(0:R-1));
                    for n = N+1:length(z)
                        J{n} = sparse(j.',i.',val);
                    end 
                end 
                this.J = vertcat(J{:});
            end 
            
            this.freq = cell(1,N);
            if this.useOldPreconditioner
                this.freq = repmat({1},1,this.N);
            else 
                for n = 1:N
                    ntot = prod(size_tens([1:n-1 n+1:N]));
                    this.freq{n} = histc(T.sub{n}, 1:size_tens(n))/ntot;
                end
            end 
            
            this.isInitialized = true;
        end 
        
        function state(this, z)
        %STATE Update the state of the kernel. 
        %  STATE(Z) updates the state for a regular iteration of the algorithm,
        %  which includes caching results that change every iteration but are
        %  use by multiple funtions.
            
            R = this.R;
            T = this.T;
            
            % Optionally cache the inverses of the Gramians for the preconditioner.
            % In a faster language, this should be the Cholesky factor instead.
            if this.usePreconditioner
                % Cache the factor matrices' Gramians.
                UHU = zeros(length(z),R*R);
                for n = 1:length(z)
                    tmp = conj(z{n}'*z{n});
                    UHU(n,:) = tmp(:);
                end
                
                this.invUHU = cell(1,length(z));
                for n = 1:length(z)
                    tmp = UHU([1:n-1 n+1:length(z)],:);
                    if length(z) > 2, tmp = prod(tmp,1); end
                    wrnstat = warning('query', 'MATLAB:nearlySingularMatrix');
                    warning('off', 'MATLAB:nearlySingularMatrix');
                    warning('');
                    this.invUHU{n} = inv(reshape(tmp,[R R]));
                    warning(wrnstat);                                            
                end
            end

            if this.useGramian || ~this.useSparseToFull
                if this.useMex
                    z = cellfun(@transpose, z, 'UniformOutput', false);
                    cpd_jacobian_incomplete(this.J, z, T.sub);
                else
                    size_tens = this.size_tens;
                    
                    J = cell(R,length(z));
                    i = (1:numel(T.val)).';
                    for n = 1:length(z);
                        dims = [1:n-1 n+1:length(z)];
                        for r = 1:R
                            val = z{dims(1)}(T.sub{dims(1)},r);
                            for k = dims(2:end)
                                if k <= this.N
                                    val = val .* z{k}(T.sub{k},r);
                                else 
                                    val = val .* z{k}(1,r);
                                end 
                            end
                            if n <= this.N
                                tmp = sparse(double(T.sub{n}),i,val, size_tens(n), ...
                                             numel(T.val));
                            else 
                                tmp = sparse(ones(size(T.sub{1})),i,val, 1, ...
                                             numel(T.val));
                            end 
                            J{r,n} = tmp;                
                        end
                    end
                    this.J = cat(1, J{:});
                end 
            end 
        end
        
        function fval = objfun(this, z)
        %OBJFUN Objective function value for CPDI. 
            T = this.T;
            N = this.N;
            R = this.R;
            size_tens = this.size_tens;
            
            if this.useSparseToFull
                fval = cpdgen(z); 
                fval = fval(T.ind)-T.val;
            else
                fval = -T.val;
                for r = 1:R
                    tmp = z{1}(T.sub{1},r);
                    for n = 2:N, tmp = tmp.*z{n}(T.sub{n},r); end
                    for n = N+1:length(z), tmp = tmp.*z{n}(1,r); end
                    fval = fval+tmp;
                end
            end
            E = T; 
            E.val = fval; 
            E.matrix = [];
            this.residual = E;

            fval = 0.5*(fval(:)'*fval(:));
        end 
             
        function grad = grad(this, z)
        %GRAD Gradient for CPDI. 
            this.state(z);
            if this.useSparseToFull
                offset = this.offset;
                grad = nan(offset(end)-1,1);
                E = this.residual;
                for n = 1:length(z)
                    tmp = full(mtkrprod(E,z,n));
                    grad(offset(n):offset(n+1)-1) = tmp(:);
                end
            else 
                grad = conj(this.J)*this.residual.val;
                grad = grad(:);
            end 
        end
        
        function JHJ = JHDJ(this, z)
        %JHJ Gramian for CPDI.
            J = this.J;
            JHJ = full(conj(J*J')); 
        end
        
        function y = JHDJx(this, z, y)
        %JHJX Gramian-vector product for CPDI.           
            y = y.'*this.J;
            y = conj(this.J*y');
            y = y(:);
        end

        function x = M_blockJacobi(this,~,b)
        %M_BLOCKJACOBI Stochastic Block-Jacobi preconditioner for CPDI.
            offset = this.offset;

            x = nan(size(b));
            for n = 1:length(offset)-1
                idx = offset(n):offset(n+1)-1;
                B = reshape(b(idx),[],size(this.invUHU{n},1));
                B = bsxfun(@rdivide, B, this.freq{n});
                tmp = B * this.invUHU{n};
                x(idx) = tmp(:);
            end
        end


        function isvalid = validate(this, z)
        %VALIDATE Check whether this kernel is in a valid.
        %  VALIDATE(U) checks whether the combination of data and initial factor
        %  matrices U is valid. This is a deep check, and may take some time,
        %  but it ensures that the computation will succeed.
            
            if ~isvalidtensor(this.T, false, 'incomplete');
                error('CPDI_Kernel:T', ['T is not a valid incomplete tensor. ' ...
                                    'See isvalidtensor.']);
            end 
            if ~iscell(z) || any(~cellfun(@ismatrix,z));
                error('CPDI_Kernel:U', 'U should be a cell of factor matrices.');
            end 
            sz = cellfun('size', z(:).', 1);
            if numel(z) ~= this.N && any(sz(this.N+1:end) ~= 1)
                error('CPDI_Kernel:U', ...
                      ['numel(U) should be getorder(T), or size(U{n},1) should ' ...
                       'be 1 for n > getorder(T).']);                
            end           
            if any(sz(1:this.N) ~= this.size_tens)
                error('CPDI_Kernel:U', ...
                      ['size(U{n},1) should be getsize(T,n) for all n.']);                
            end 
            R = size(z{1},2);
            if any(R ~= cellfun('size', z, 2))
                error('CPDI_Kernel:U', ...
                      'size(U{n},2) should be R for all n.');     
            end 
            
            isvalid = true;
        end 
    end        
    
    methods (Access = protected)
                
        function header = getHeader(this)
            if ~isscalar(this)
                header = getHeader@matlab.mixin.CustomDisplay(this);
            else
                className = matlab.mixin.CustomDisplay.getClassNameForHeader(this);
                newHeader = [className, ' CPD kernel for incomplete tensor'];
                header = sprintf('%s\n',newHeader);
            end
        end
        
        function propgrp = getPropertyGroups(this)
            if ~isscalar(this)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(this);
            else
                dataprops = {'T', 'size_tens', 'R'};
                propgrp(1) = matlab.mixin.util.PropertyGroup(dataprops, ...
                                                             'Data properties');                
                settingprops = {'useGramian', 'usePreconditioner'};
                propgrp(2) = matlab.mixin.util.PropertyGroup(settingprops, ...
                                                             'Computational settings');
                specializedprops = {'useOldPreconditioner', 'useSparseToFull', ...
                                   'useMex'};
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
