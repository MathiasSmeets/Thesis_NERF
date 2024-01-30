classdef CPDBasic_Kernel < TensorOptimizationKernel
%CPD_KERNEL Computational routines for CPD
%
%   CPD_KERNEL(T, z) is a class collecting computational routines for the
%   canonical polyadic decomposition of a Nth order tensor tensor T into factors
%   z{1}, ..., z{N}. This class stores intermediate results for improved
%   efficiency.
%
%   Let kernel be an instance of the CPD_KERNEL class. The following objective
%   function routines can be used:
%    
%       - kernel.objfun(z)
%
%   The following gradient functions can be used:
%  
%       - kernel.grad(z)
%
%   The following Gramian functions can be used:
%
%       - kernel.JHJ(z)
%
%   The following Gramian-vector product functions can be used:
%
%       - kernel.JHJx(z, x)      For dense, sparse, incomplete and structured
%                                tensors 
%       - kernel.JHJxi(z, x)     For incomplete tensors only.
%       - kernel.JHJxli(z, x)    For incomplete tensors with linearly
%                                constrained factors only.
%
%   The following preconditioners can be used:
%
%       - kernel.M_blockJacobi(~, b)  
%       - kernel.M_blockSSOR(z, b)  
%
%   The following line and plane search methods can be used
%   TODO
    
    properties
        %% Tensor properties
        
        T;                 % data
        R;                 % rank
        N;                 % order
        size_tens;         % size of data
        type;              % type 
        isstructured;      % true if type is not full, sparse or incomplete
        isincomplete;      % true if incomplete tensor
        
        %% Cached variables
        
        residual;          % residual
        offset;
        UHU;
        Wn                 % diagonal blocks
        invW;
        T2;                % squared Frobenius norm of tensor
        scale;             % fraction of known entries
        
        %% Settings 
        
        % Compute cached results for Gramians 
        useGramian = true; 
        % Compute cached results for preconditioner 
        usePreconditioner = true;
        
        P
        freq
        
        W
        widx
    end 
    
    
    methods
        
        function this = CPDBasic_Kernel(T, varargin)
            import tensorlab.auxiliary.*;
            
            p = inputParser;
            p.parse(varargin{:});
            p.KeepUnmatched = true;
            options = p.Results;

            this.handlesSymmetry  = false;
            this.handlesConstants = false;
            
            this.T = T;
            this.size_tens = getsize(T);
        end 
        
        function state(this, z)    
        % Cache the factor matrices' Gramians.
            N = this.N;
            R = this.R;
            
            if this.useGramian 
                for n = 1:N
                    tmp = conj(z{n}'*z{n});
                    this.UHU(:,:,n) = tmp;
                end
            end             
            if this.usePreconditioner
                this.invW = cell(1,N);
                for n = 1:N
                    Wn = prod(this.UHU(:,:,[1:n-1 n+1:N]), 3);
                    wrnstat1 = warning('off', 'MATLAB:nearlySingularMatrix');
                    wrnstat2 = warning('off', 'MATLAB:singularMatrix');
                    warning('');
                    this.invW{n} = inv(Wn);
                    [~,id] = lastwarn;
                    if ~isempty(id)
                        this.invW{n} = pinv(this.Wn{n});
                    end 
                    warning(wrnstat1);
                    warning(wrnstat2);
                end
            end 
        end 
        
        function initialize(this, z)
            this.R      = size(z{1},2);
            this.N      = numel(z);
            this.offset = cumsum([1 cellfun(@numel,z(:).')]);
            
            if this.useGramian
                this.UHU = zeros(this.R,this.R,this.N);
            end
        end
        
        function fval = objfun(this, z)
            fval = cpdres(this.T, z, 'Format', false);
            this.residual = reshape(fval, this.size_tens);
            fval = 0.5*(fval(:)'*fval(:));            
        end 
        
        function grad = grad(this, z)
            this.state(z); 
            offset = this.offset;
            grad   = zeros(offset(end)-1,1);

            % CPD scaled conjugate cogradient.
            for n = 1:this.N
                idxn = offset(n):offset(n+1)-1;
                tmp  = mtkrprod(this.residual, z, n);
                grad(idxn) = tmp(:);
            end
        end
        
        function JHJ = JHDJ(this, z)
            N = numel(z);
            R = this.R; 
            size_z = this.size_tens;
            offset = this.offset;
            
            % CPD Jacobian's Gramian.
            UHU = conj(this.UHU);
            JHJ = zeros(offset(end)-1);
            for n = 1:N
                % diagonal blocks
                idxn = offset(n):offset(n+1)-1;
                Wn   = prod(UHU(:,:,[1:n-1 n+1:N]), 3);
                JHJ(idxn,idxn) = kron(Wn,eye(size_z(n)));
                for m = n+1:N
                    % off-diagonal blocks
                    idxm  = offset(m):offset(m+1)-1;
                    Wnm   = prod(UHU(:,:,[1:n-1 n+1:m-1 m+1:N]), 3);
                    JHJnm = bsxfun(@times,reshape(z{n},[size_z(n) 1 1 R]), ...
                                   reshape(z{m}',[1 R size_z(m) 1]));
                    JHJnm = bsxfun(@times,JHJnm,reshape(Wnm,[1 R 1 R]));
                    JHJnm = reshape(JHJnm, size_z(n)*R, size_z(m)*R);
                    JHJ(idxn,idxm) = JHJnm;
                    JHJ(idxm,idxn) = JHJnm';
                end
            end
        end 
        
        function y = JHDJx(this, z, x)
            N = this.N;
            R = this.R;
            offset = this.offset;
            UHU = this.UHU;
            XHU = zeros(R,R,N);
            y = zeros(offset(end)-1,1);
            
            for n = 1:N
                Wn = prod(UHU(:,:,[1:n-1 n+1:N]),3);
                if iscell(x), tmp = x{n};
                else tmp = reshape(x(offset(n):offset(n+1)-1),size(z{n})); end                
                XHU(:,:,n) = conj(tmp'*z{n});
                idxn    = offset(n):offset(n+1)-1;
                y(idxn) = tmp * Wn;
            end
            for n = 1:N
                idxn = offset(n):offset(n+1)-1;
                Wn   = zeros(R,R);
                for m = n+1:N                    
                    idxm = offset(m):offset(m+1)-1;
                    if numel(this.size_tens) == 2
                        Wn = Wn+XHU(:,:,m);
                        JHJmnx = z{m}*XHU(:,:,n);
                    else
                        Wnm = prod(UHU(:,:,[1:n-1 n+1:m-1 m+1:N]), 3);
                        Wn  = Wn+Wnm.*XHU(:,:,m);
                        Wnm = Wnm.*XHU(:,:,n);
                        JHJmnx = z{m}*Wnm;
                    end
                    y(idxm) = y(idxm) + JHJmnx(:);
                end
                JHJnx = z{n}*Wn;
                y(idxn) = y(idxn) + JHJnx(:);
            end
        end

        function x = M_blockJacobi(this,~,b)
            x = zeros(size(b));
            for n = 1:this.N
                idx = this.offset(n):this.offset(n+1)-1;
                tmp = reshape(b(idx),[],this.R)*this.invW{n};
                x(idx) = tmp(:);
            end
        end
        
        function isvalid = validate(this, z)
            import tensorlab.auxiliary.*
            
            this.initialize(z);
            if ~isnumeric(fmt(this.T))
                error('tok:invaliddata', 'The given tensor should be a dense tensor');
            end 
            if ~isvalidtensor(this.T) 
                error('tok:invaliddata', 'The given tensor is not valid.');
            end 
            if ~iscellmat(z) 
                error('tok:invalidmodel', 'U should be a cell of matrices');
            end 

            sz1 = getsize(this.T);
            sz2 = getsize(z);
            if ~sizeequals(sz1, sz2)
                error('tok:invalidmodel', ...
                      ['The dimensions of the model (%s) do not match the ' ...
                       'dimensions of the data (%s).'], size2str(sz2), size2str(sz1));
            end
            
            if ~all(cellfun('size', z, 2) == size(z{1},2))
                error('tok:invalidmodel', 'size(U{n},2) should be equal for all n.');
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
                newheader = ['Basic CPD kernel routines for simple CPD of a dense tensor.'];
                header = sprintf('%s %s\n', className, newheader);
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
                specializedprops = {'handlesConstants', 'handlesSymmetry'};
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