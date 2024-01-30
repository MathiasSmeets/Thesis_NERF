classdef CPD_Kernel < TensorOptimizationKernel
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
    
    properties (SetAccess = protected)
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
        offset;            % offset of each variable in vector
        UHU;               % Cached inner products of variables
        Wn                 % diagonal blocks of Gramian in compressed form
        invW;              % inverse of Gramian blocks, or L factor
        T2;                % squared Frobenius norm of tensor

        %% Sumport for right multiplication constraints
        % Each factor matrix can have a right matrix multiplication constraint, e.g., to model an
        % LL1 decomposition. 
        
        % Cell with matrix contraints (one for each variable). Scalar 1 if unused
        P
        % If no linear constraints are used, disable computions
        useLinearConstraints = false;
        
        %% Support for symmetry constraints
        
        % Auxiliary vector indicating the number of times a variable is used. 
        freq
        
        %% Support for rank-1 weight tensor
        % A rank-1 weight, given as in the CPD format can be used. Virtual variables are used to
        % handle the (specific) case of symmetry constraints with unsymmetric weights.
        
        % cell of length order(T) with weight vectors; empty if no weight
        W
        % vector of length order(T) indicating a weight is used 
        widx
        % Heuristic weights in the case of incomlete tensors
        Winc
        
        %% Definition of virtual variables
        % Virtual variables are required to implement rank-1 weights when using symmetry and when
        % using constants that are not given as variables. 
        
        % Mapping of each (virtual) variable to mode 
        var2mode 
        % Mapping of each mode to a (virtual) variable
        mode2var
        % Mapping of real to virtual variables
        var2virt 
        
        %% Support for constants
        % Constants factor matrices are supported in two ways: given as a variable in z, and in a
        % separate cell, which is defined below.
        
        % Constants that are not given as variables (cell of length order(T))
        C;  
        % Mapping of each mode to a constant (0 otherwise)
        constantsInMode;
    end 
    
    properties
        %% Computational Settings 
        
        % Compute cached results for Gramians 
        useGramian = true; 
        % Compute cached results for preconditioner 
        usePreconditioner = true;
        % Use LDL factorization based preconditioner (only for unsymmetric weights)
        useLDLPrec = false;
        % Use structured approach as long as possible for full tensors
        useFastUpdate = false;
    end 
    
    methods
        
        function this = CPD_Kernel(T, varargin)
            import tensorlab.auxiliary.*;
            
            p = inputParser;
            p.addOptional('P', {});
            p.addOptional('L', nan);
            p.addOptional('Lmodes','last');
            p.addOptional('R1Weights', []);
            p.addOptional('VariablesInMode',[]);
            p.addOptional('Constants', []);
            p.addOptional('ConstantsInMode',[]);
            p.addOptional('IsConstant', []);
            p.addOptional('IsDataSymmetric', 'auto');
            p.addOptional('IsStructured', 'auto');
            p.parse(varargin{:});
            p.KeepUnmatched = true;
            options = p.Results;

            this.handlesSymmetry  = true;
            this.handlesConstants = true;
            
            N                 = getorder(T);
            this.T            = T;
            this.size_tens    = getsize(T);
            this.type         = getstructure(T);
            this.isstructured = ~any(strcmpi(this.type, {'full', 'incomplete'}));    
            if ~this.isstructured
                this.T = fmt(T);
                if ~isnumeric(this.T)
                    this.type = getstructure(this.T);
                end 
                this.isstructured = ~any(strcmpi(this.type, {'full', 'incomplete'}));
            end 
            % Override if given
            if ~ischar(options.IsStructured) || ~strcmpi(options.IsStructured, 'auto')
                this.isstructured = options.IsStructured;
            end 
            
            this.isincomplete = strcmpi(this.type, 'incomplete');            
            if this.isincomplete
                this.Winc = cell(1, N);
                for n = 1:N
                    ntot = prod(this.size_tens([1:n-1 n+1:N]));
                    this.Winc{n} = histc(this.T.sub{n}, 1:this.size_tens(n))/ntot;
                end
            end 
            
            vidx = options.VariablesInMode(:).';
            if isempty(vidx), vidx = 1:N; end 
            if any(strcmpi('VariablesInMode', p.UsingDefaults))
                vidx = 1:N;
            else 
                this.variablesInMode = vidx(:).';
                if any(~isnonnegint(vidx))
                    error('tok:invalidOptions', 'VariablesInMode should contain nonnegative integers.');
                end
                if numel(vidx) < N
                    error('tok:invalidOptions', ['numel(VariablesInMode)=%d should contain at least ' ...
                                        'getorder(T)=%d entries.'], numel(vidx), N);
                end 
            end 
            
            %% Constants that are not variables 
            if ~any(strcmpi('ConstantsInMode', p.UsingDefaults)) && any(options.ConstantsInMode)
                cidx = options.ConstantsInMode(:).';
                if any(~isnonnegint(cidx))
                    error('tok:invalidOptions', 'ConstantsInMode should contain nonnegative integers.');
                end
                if any(strcmpi('VariablesInMode', p.UsingDefaults))
                    error('tok:invalidOptions', ['VariablesInMode should be specified if ' ...
                                        'ConstantsInMode is given']);
                end 
                this.C = options.Constants;
                if ~iscell(this.C)
                    error('tok:invalidOptions', ['The Constants option should be given if ' ...
                                        'ConstantsInMode is defined']);
                end 
                if ~iscellmat(this.C)
                    error('tok:invalidOptions', 'Constants should be a cell of matrices.');
                end 
                if numel(unique(cidx(cidx>0))) ~= numel(this.C) || max(cidx) ~= numel(this.C)
                    error('tok:invalidOptions', 'numel(Constants) should match the max(ConstantsInMode).');
                end 
                if numel(vidx) > numel(cidx)
                    cidx = [cidx(:).', zeros(1,numel(vidx)-numel(cidx))];
                end 
                this.constantsInMode = cidx;
                if numel(cidx) > numel(vidx)
                    error('tok:invalidOptions', ['numel(ConstantsInMode) should be equal to ' ...
                                        'numel(VariablesInMode).']);
                end 

                % Check sizes
                sz1 = getsize(this.T, find(cidx > 0));
                sz2 = cellfun('size', this.C, 1);
                sz2 = sz2(cidx(cidx>0));
                if ~sizeequals(sz1,sz2)
                    error('tok:invalidOptions', ['size(Constants{ConstantsInMode(n)},1) should equal ' ...
                                        'getsize(T,n) for all n with ConstantsInMode(n) ~= 0.']);
                end
                if any(cellfun('size', this.C, 2) ~= size(this.C{1}, 2))
                    error('tok:invalidOptions', ['size(Constants{n},2) = R should hold for all n ' ...
                                        '= 1,...,numel(Constants)']);
                end
                % Check if all modes have a variable
                if ~all(cidx | vidx)
                    error('tok:invalidOptions', ['A variable or a constant should be defined for ' ...
                                        'each mode: VariablesInMode(n) > 0 or ConstantsInMode(n) ' ...
                                        '> 0']);
                end 
                if any(cidx & vidx) 
                    error('tok:invalidOptions', ['Either a variable or a constant should be defined ' ...
                                        'for each mode: VariablesInMode(n) > 0 if ConstantsInMode(n) ' ...
                                        '== 0 and vice versa.']);
                end 
            end 
            
            %% Symmetry constraints
            if ischar(options.IsDataSymmetric) && strcmpi(options.IsDataSymmetric, 'auto')
                if ~strcmpi(this.type, 'full')
                    this.isDataSymmetric = false; 
                else 
                    issym = true;
                    freq = histc(vidx(vidx>0),1:max(vidx));
                    for k = find(freq > 1)
                        ind = find(vidx == k);
                        T1 = tens2mat(this.T, ind(1));
                        for l = 2:numel(ind)
                            dims = 1:N;
                            dims(ind(l)) = ind(1);
                            dims(ind(1)) = [];
                            T2 = tens2mat(this.T, ind(l), dims);
                            if any(abs(T1(:) - T2(:)) > 1e-13)
                                issym = false;
                                break;
                            end 
                        end 
                        if ~issym, break; end
                    end 
                    this.isDataSymmetric = issym;
                end 
            end
            
            %% LL1 type constraints 
            if ~any(strcmpi('L', p.UsingDefaults))
                if ~isempty(options.P)
                    error('tok:invalidOptions', 'The P and L options cannot be specified simultaneously.');
                end 
                if ischar(options.Lmodes)
                    if strcmpi(options.Lmodes, 'first'),    options.Lmodes = 1; 
                    elseif strcmpi(options.Lmodes, 'last'), options.Lmodes = N;
                    else 
                        error('tok:invalidOptions', ...
                              'Unknown L mode ''%s'' (Lmodes can be ''first'', ''last'' or an array).', ...
                              options.Lmodes);
                    end 
                else
                    if any(~isnonnegint(options.Lmodes, true)) || any(options.Lmodes>N)
                        error('tok:invalidOptions', ['Lmodes should contain integers ' ...
                                            'between 1 and getorder(T), or should ''first'' or ' ...
                                            '''last''.']);
                    end 
                    if any(diff(sort(options.Lmodes)) == 0)
                        error('tok:invalidOptions', 'Lmodes should contain unique integers.');
                    end 
                end 
                if ~iscell(options.L), options.L = {options.L}; end 
                if numel(options.L) ~= numel(options.Lmodes)
                    error('tok:invalidOptions', 'numel(L) should equal numel(Lmodes).');
                end 
                if any(sum(options.L{1}) ~= cellfun(@sum,options.L))
                    error('tok:invalidOptions', 'sum(L{n}) should be equal for all n.');
                end 
                if ~isempty(intersect(find(options.ConstantsInMode), options.Lmodes))
                    error('tok:invalidOptions', ['Linear constraints are not allowed for constants ' ...
                                        'defined in C. (find(ConstantsInMode) and Lmodes should ' ...
                                        'not intersect.)']);
                end 
                P = repmat({1}, 1, N);
                for k = 1:numel(options.L)
                    ind = zeros(1,sum(options.L{k}));
                    ind(cumsum([1 options.L{k}(1:end-1)])) = 1;
                    P{options.Lmodes(k)} = sparse(cumsum(ind), 1:sum(options.L{k}), 1);
                end 
                % Check values in case of symmetry
                [v,i,~] = unique(vidx);
                i = i(1+(v(1)==0):end); % remove ref to first if vidx contains a zero.
                for k = v(:).'
                    ind = find(vidx == k);
                    for l = 2:numel(ind)
                        if any(size(P{ind(k)}) ~= size(P{ind(1)})) || any(P{ind(k)}(:) ~= P{ind(1)}(:))
                            error('tok:invalidOptions', ['Conflicting values for L are found because ' ...
                                                'of symmetry as defined by VariablesInMode.']);
                        end
                    end 
                end 
                P = P(i);
                this.P = P;
            elseif ~any(strcmpi('Lmodes', p.UsingDefaults))
                error('tok:invalidOptions', 'Lmodes defined, but no L given.');
            elseif ~isempty(options.P) 
                P = options.P;
                if numel(P) < max(vidx)
                    P = [P, repmat({1},1,numel(P)-max(vidx))];
                elseif numel(P) > max(vidx)
                    error('tok:invalidOptions', ['One constraint per variable should be provided: ' ...
                                        'numel(P) should be numel(z) or max(VariablesInMode).']);
                end 
                ind = cellfun(@(p) ~isempty(p) && ~isscalar(p), P);
                if numel(ind) > 0 && any(diff(cellfun('size',P(ind),2)))
                    error('tok:invalidOptions', ['size(P{n},2) should be equal for all nonscalar ' ...
                                        'P{n}, n=1,...,numel(z).']);
                end 
                this.P = P;
            else
                this.P = {};
            end
            
            %% Constant variables 
            if ~isempty(options.IsConstant)
                isconst = options.IsConstant;
                if numel(isconst) < max(vidx)
                    isconst = [isconst(:).' zeros(1,max(vidx) - numel(isconst))];
                elseif numel(isconst) > max(vidx) && ~any(strcmpi('VariablesInMode', p.UsingDefaults))
                    error('tok:invalidOptions', ['numel(IsConstant)=%d should equal ' ...
                                        'max(VariablesInMode)=%d.'], numel(isconst), max(vidx));
                end
                this.isConstant = isconst;
            end
            
            %% Weights
            this.W = options.R1Weights(:).';
            this.prepare();
        end 

        function prepare(this)
            if isempty(this.W), this.W = repmat({[]},1,getorder(this.T)); end
            N = getorder(this.T);
            if numel(this.W) < N
                this.W = [this.W, repmat({[]}, 1, N - numel(this.W))];
            end 
            if numel(this.W) ~= N
                error('tok:invalidweight', ['numel(W) should equal getorder(T).']);   
            end 
            this.widx = ~cellfun(@(w) isempty(w) || (isscalar(w) && w == 1), this.W);
            ind     = this.widx;
            weights = this.W(ind);
            sz1     = this.size_tens(ind);
            if any(~cellfun(@isvector,weights))
                error('tok:invalidweight', ['W{n} should be a column vector for n=1,...,' ...
                                    'getorder(T) and weights W.']);
            end 
            this.W = cellfun(@(w) w(:), this.W, 'UniformOutput', false);
            if any(sz1 ~= cellfun('size',weights,1))
                error('tok:invalidweight', ['numel(W{n}) should equal getsize(T,n) for n=1,...,' ...
                                    'getorder(T) and weights W.']);
            end           
            if any(cellfun(@(w) any(w<0), weights))
                error('tok:invalidweight', ['W{n} should be nonnegative for n=1,...,' ...
                                    'getorder(T) and weights W.']);
            end 
            
            if ~this.isstructured && any(this.widx)
                if isstruct(this.T) % incomplete 
                    for k = find(this.widx)
                        this.T.val = this.T.val .* this.W{k}(this.T.sub{k});
                    end 
                else 
                    tomode = @(v,n) reshape(v, [ones(1,n-1), numel(v) 1]);
                    for k = find(this.widx)
                        this.T = bsxfun(@times, this.T, tomode(this.W{k},k));
                    end 
                end 
            end
        end 
        
        function initialize(this, z)
            import tensorlab.auxiliary.vecequals;

            if isempty(this.P)
                this.P = repmat({1}, 1, numel(z)); 
            elseif ~iscell(this.P)
                error('tok:P', 'P should be empty or a cell of length numel(U).');                
            end 
            this.P(cellfun(@isempty, this.P)) = {1};
            if isscalar(this.P{1});
                this.R = size(z{1},2);
            else 
                this.R = size(this.P{1},2);
            end 
            this.useLinearConstraints = any(~cellfun(@isscalar, this.P));
            this.N      = numel(z);
            this.offset = cumsum([1 cellfun(@numel,z(:).')]);
            
            if isempty(this.variablesInMode)
                this.variablesInMode = 1:numel(z);
            else
                this.variablesInMode = this.variablesInMode(:).';
            end 
                        
            % Correct weights for 1 modes
            if numel(this.widx) < numel(this.variablesInMode)
                this.widx = [this.widx false(1, numel(this.variablesInMode) - numel(this.widx))];
            end 
            
            %% Create virtual variables if needed
            vidx = this.variablesInMode;
            [nidx,var2mode_old,i] = unique(vidx + (max(vidx) + 1) .* (vidx == 0));
            newvaridx = numel(z)+1;
            var2virt = [];
            var2mode = [];
            mode2var = vidx;
            for k = 1:numel(nidx)
                if nidx(k) == max(vidx) + 1, continue; end
                idx = find(i == nidx(k));
                if numel(idx) > 1 && any(cellfun(@(w) ~vecequals(this.W{idx(1)},w), this.W(idx(2:end))))
                    % if this variable is used in multiple modes, but rank-1 weights are not
                    % identical, create virtual variables 
                    var2mode = [var2mode idx(:).'];
                    var2virt = [var2virt ones(1,numel(idx))*nidx(k)];
                    mode2var(idx) = numel(var2virt) + (-numel(idx)+1:0);
                    this.useLDLPrec = true;
                else 
                    var2virt = [var2virt nidx(k)];
                    var2mode = [var2mode idx(1)];              
                    mode2var(idx) = numel(var2virt);
                end 
            end
            % Remove symmetry in constants if weights are not symmetric
            cidx = this.constantsInMode;
            for k = 1:max(cidx)
                idx = find(cidx == k);
                if numel(idx)>1 && any(cellfun(@(w) ~vecequals(this.W{idx(1)},w), this.W(idx(2:end))))
                    this.C = [this.C, repmat(this.C(k), 1, numel(idx)-1)];
                    this.constantsInMode(idx) = [cidx(idx(1)) max(this.constantsInMode) + (1:numel(idx)-1)];
                end 
            end 
            % Add virtual variables for constants that are not variables
            nbvar         = numel(var2virt);
            var2virt      = [var2virt, max(var2virt)+(1:numel(this.C))];
            cidx          = this.constantsInMode;
            ind           = cidx > 0;
            mode2var(ind) = cidx(ind) + nbvar;
            [v,i] = unique(cidx + (numel(this.C) + 1).*(cidx==0));
            var2mode      = [var2mode i(v<numel(this.C)+1).'];
            if numel(this.P) ~= numel(z) + numel(this.C)
                this.widx     = [this.widx, zeros(1, sum(ind))];
                this.P        = [this.P, repmat({1}, 1, numel(this.C))];
            end 
            % Store results
            this.var2virt = var2virt;
            this.var2mode = var2mode;
            this.mode2var = mode2var;
            this.freq     = histc(mode2var, unique(mode2var));
                        
            if this.useGramian
                this.UHU = zeros(numel(var2virt),this.R^2);
            end
            
            if isempty(this.isConstant)
                this.isConstant = [false(1, numel(z)), true(1,numel(this.C))];
            elseif numel(this.isConstant) == numel(z)
                this.isConstant = [this.isConstant, true(1,numel(this.C))];
            end
            
            if sum(this.widx) > 0
                this.T2 = frob(this.T, 'squared', this.W);
            else 
                this.T2 = frob(this.T, 'squared');
            end

        end

        function state(this, z)    
        % Cache the factor matrices' Gramians.
            N = this.N;
            R = this.R;
            var2virt     = this.var2virt;
            var2mode     = this.var2mode;
            size_tens    = this.size_tens;
            widx         = this.widx(var2mode);
            if this.useGramian && this.isincomplete
                tmp  = this.Winc;
                ind  = this.widx;
                tmp(ind) = cellfun(@times, tmp(ind), this.W(ind), 'UniformOutput', false);
            else 
                tmp = this.W;
            end 
            weights      = cell(1,N); 
            ind          = var2mode <= numel(this.W);
            weights(ind) = this.W(var2mode(ind));
            
            if this.useGramian 
                if numel(this.C) > 0, z = [z, this.C]; end
                z = z(var2virt);
                P = this.P(var2virt);
                for n = 1:numel(z)
                    if widx(n)
                        z{n} = bsxfun(@times, weights{n}, z{n}); 
                    end
                    Wn = conj(z{n}'*z{n});
                    if this.useLinearConstraints
                        Wn = P{n}.'*Wn*conj(P{n});
                    end 
                    this.UHU(n,:) = Wn(:);
                end
            end             
            if this.usePreconditioner
                this.invW = cell(1,N);
                size_p = cellfun('size', z, 2);
                for n = 1:N
                    if this.isConstant(n), continue, end;
                    if this.useLDLPrec
                        modes = find(var2virt == n);
                        M = 0;
                        for k = modes 
                            mult = this.freq;
                            mult(k) = mult(k) - 1;
                            if any(mult > 1), Wn = this.UHU .^ (mult(:) * ones(1, R^2));
                            else, Wn = this.UHU(mult==1,:); end 
                            Wn = reshape(prod(Wn, 1), [R R]);
                            Wn = this.freq(n) * conj(P{k}*Wn*P{k}');
                            if widx(n)
                                M = M + kr(Wn, weights{k}.^2*ones(1,size_p(k)));
                            else 
                                M = M + kr(Wn, ones(size(z{k})));
                            end 
                        end 
                        this.Wn{n}   = M;
                        this.invW{n} = ldldblk(M);
                    else 
                        mult = this.freq;
                        mult(n) = mult(n) - 1;
                        if any(mult > 1), Wn = this.UHU .^ (mult(:) * ones(1, R^2));
                        else, Wn = this.UHU(mult==1, :); end
                        Wn = reshape(prod(Wn, 1), [R R]);
                        this.Wn{n} = this.freq(n) * (P{n}*Wn*P{n}');
                        % Try to use inv (faster) unless nearly singular.
                        wrnstat1 = warning('off', 'MATLAB:nearlySingularMatrix');
                        wrnstat2 = warning('off', 'MATLAB:singularMatrix');
                        warning('');
                        this.invW{n} = inv(this.Wn{n});
                        [~,id] = lastwarn;
                        if ~isempty(id)
                            this.invW{n} = pinv(this.Wn{n});
                        end 
                        warning(wrnstat1);
                        warning(wrnstat2);
                    end 
                end
            end 
        end 
        
        function fval = objfun(this, z)
            import tensorlab.auxiliary.vecequals;
            
            if numel(this.C) > 0, z = [z this.C]; end 
            z = cellfun(@mtimes, z, this.P, 'UniformOutput', false);
            if ~vecequals(this.variablesInMode, 1:numel(z))
                z = z(this.var2virt(this.mode2var));
            end
            for k = find(this.widx)
                z{k} = bsxfun(@times, this.W{k}, z{k});
            end 
            
            if this.isstructured
                ztmp = z;
                for k = find(this.widx)
                    ztmp{k} = bsxfun(@times, this.W{k}, z{k});
                end 
                fval = abs(0.5*this.T2 - real(inprod(this.T, ztmp)) + 0.5*frob(z,'squared'));
                this.isstructured = ~this.useFastUpdate || abs(fval) > 1e-8*abs(this.T2);
                if this.isstructured, return; end
            end 
            fval = cpdres(this.T, z, 'Format', false);
            if isstruct(fval)
                this.residual = fval;
                fval = fval.val;
                fval = 0.5*(fval(:)'*fval(:));
            else 
                this.residual = reshape(fval, this.size_tens);
                fval = 0.5*(fval(:)'*fval(:));            
            end 
        end 
        
        function grad = grad(this, z)
            import tensorlab.auxiliary.vecequals;

            this.state(z); 
            offset   = this.offset;
            grad     = zeros(offset(end)-1,1);
            updated  = false(1, this.N);
            var2mode = this.var2mode;
            var2virt = this.var2virt;
            freq     = this.freq;

            if this.useLinearConstraints
                z = cellfun(@mtimes, z, this.P, 'UniformOutput', false);
            end 
            if numel(this.C) > 0, z = [z this.C]; end 
            z = z(var2virt);
            if ~vecequals(this.mode2var, 1:numel(z))
                z = z(this.mode2var); 
            end
            for k = find(this.widx)
                z{k} = bsxfun(@times, this.W{k}, z{k});
            end 
            if this.isstructured
                z2 = z;
                for k = find(this.widx)
                    z2{k} = bsxfun(@times, this.W{k}, z{k});
                end 
            end
            
            % If the data is not symmetric, all blocks need to be computed and summed, otherwise
            % multiplications can be used to save some computations.
            if ~this.isDataSymmetric
                var2virt = this.variablesInMode;
                var2mode = 1:numel(var2virt);
                freq = ones(1, numel(var2virt));
            end 
            % Correct P
            P = cell(1, length(var2virt));
            ind = var2virt > 0;
            P(ind) = this.P(var2virt(ind));
            
            % CPD scaled conjugate cogradient.
            for n = find(ind)
                if this.isConstant(var2virt(n)), continue; end 
                if this.isstructured 
                    tmp = mtkrprod(this.T,z2,var2mode(n));
                    if this.widx(var2mode(n)), tmp = bsxfun(@times,this.W{var2mode(n)}, tmp); end 
                    tmp = mtkrprod(z,z,var2mode(n)) - tmp;
                else 
                    tmp = full(mtkrprod(this.residual, z, var2mode(n)));
                end 
                if this.widx(var2mode(n))
                    tmp = bsxfun(@times, this.W{var2mode(n)}, tmp);
                end 
                if this.useLinearConstraints
                    tmp = tmp * P{n}';
                end 
                idxn = offset(var2virt(n)):offset(var2virt(n)+1)-1;
                if ~updated(var2virt(n))
                    grad(idxn) = freq(n) * tmp(:);
                else 
                    grad(idxn) = grad(idxn) + freq(n) * tmp(:);
                end 
                updated(var2virt(n)) = true;
            end
        end
        
        function JHJ = JHDJ(this, z)
            freq         = this.freq;
            offset       = this.offset;
            var2mode     = this.var2mode;
            var2virt     = this.var2virt;
            N            = numel(var2virt);
            R            = this.R;
            updated      = false(this.N, this.N);
            const        = this.isConstant(var2virt);
            widx         = this.widx(var2mode);
            P            = this.P(var2virt);
            uselincon    = this.useLinearConstraints;
            if this.isincomplete
                tmp  = this.Winc;
                ind  = this.widx;
                tmp(ind) = cellfun(@times, tmp(ind), this.W(ind), 'UniformOutput', false);
            else 
                tmp = this.W;
            end 
            weights      = cell(1,N); 
            ind          = var2mode <= numel(this.W);
            weights(ind) = this.W(var2mode(ind));

            if numel(this.C) > 0, z = [z this.C]; end 
            size_p = cellfun('size', z, 2);
            if uselincon
                z = cellfun(@mtimes, z, this.P, 'UniformOutput', false);
            end 
            z = z(var2virt);
            size_z = cellfun('size', z, 1);
            size_p = size_p(var2virt);
            for k = find(widx)
                z{k} = bsxfun(@times, weights{k}.^2, z{k});
            end 
            
            % CPD Jacobian's Gramian.
            UHU = conj(this.UHU);
            JHJ = zeros(offset(end)-1);
            for n = find(~const)
                idxn = offset(var2virt(n)):offset(var2virt(n)+1)-1;
                mult = [freq(1:n-1) freq(n)-1 freq(n+1:end)]; 
                
                if any(mult > 1), Wn = UHU .^ (mult(:)*ones(1,R^2));
                else Wn = UHU(mult==1,:); end 
                Wn = reshape(prod(Wn, 1), [R, R]);
                Wn = freq(n) * (conj(P{n}) * Wn * P{n}.');
                if widx(n), tmp = kron(Wn,diag(weights{n}.^2));
                else, tmp = kron(Wn,eye(size_z(n))); end 
                if updated(var2virt(n),var2virt(n))
                    JHJ(idxn,idxn) = JHJ(idxn,idxn) + tmp;
                else 
                    JHJ(idxn,idxn) = tmp;
                end 
                updated(var2virt(n),var2virt(n)) = true;
                for m = n+(freq(n)==1):N
                    if const(m), continue; end
                    mult = freq; 
                    mult([n m]) = mult([n m]) - 1 - (m==n);
                    
                    idxm = offset(var2virt(m)):offset(var2virt(m)+1)-1;
                    if any(mult > 1), Wnm = UHU .^ (mult(:)*ones(1,R^2));
                    else Wnm = UHU(mult==1,:); end 
                    Wnm = reshape(prod(Wnm, 1), [R, R]);
                    Wnm = Wnm * freq(n) * (freq(m) - (n==m));
                    JHJnm = bsxfun(@times,reshape(z{n},[size_z(n) 1 1 R]), ...
                                   reshape(z{m}',[1 R size_z(m) 1]));
                    JHJnm = bsxfun(@times,JHJnm,reshape(Wnm,[1 R 1 R]));
                    if uselincon
                        if ~isscalar(P{m})
                            JHJnm = reshape(JHJnm, [], R) * P{m}.';
                            JHJnm = reshape(JHJnm, [size_z(n), R, size_z(m), size_p(m)]);
                        end 
                        if ~isscalar(P{n})
                            JHJnm = reshape(permute(JHJnm, [1 3 4 2]), [], R);
                            JHJnm = JHJnm * P{n}';
                            JHJnm = reshape(JHJnm, [size_z(n),size_z(m),size_p(m),size_p(n)]);
                            JHJnm = permute(JHJnm, [1 4 2 3]);
                        end 
                    end 
                    JHJnm = reshape(JHJnm, size_z(n)*size_p(n), size_z(m)*size_p(m));
                    if n == m
                        JHJ(idxn,idxn) = JHJ(idxn,idxn) + JHJnm;
                    elseif updated(var2virt(n),var2virt(m)) 
                        JHJ(idxn,idxm) = JHJ(idxn,idxm) + JHJnm; 
                        JHJ(idxm,idxn) = JHJ(idxm,idxn) + JHJnm';
                    else 
                        JHJ(idxn,idxm) = JHJnm;  
                        JHJ(idxm,idxn) = JHJnm'; 
                    end 
                    updated(var2virt(n),var2virt(m)) = true;
                end
            end
        end 
        
        function y = JHDJx(this, z, x)
            freq         = this.freq;
            offset       = this.offset;
            UHU          = this.UHU;
            var2mode     = this.var2mode;
            var2virt     = this.var2virt;
            N            = numel(this.var2virt);
            R            = this.R;
            const        = this.isConstant(var2virt);
            widx         = this.widx(var2mode);
            P            = this.P(var2virt);
            uselincon    = this.useLinearConstraints;
            % Handle incompleteness
            if this.isincomplete
                tmp  = this.Winc;
                ind  = this.widx;
                tmp(ind) = cellfun(@times, tmp(ind), this.W(ind), 'UniformOutput', false);
            else 
                tmp = this.W;
            end 
            weights      = cell(1,N); 
            ind          = var2mode <= numel(this.W);
            weights(ind) = this.W(var2mode(ind));
            % Preallocation
            XHU          = zeros(R,R,N);
            y            = zeros(offset(end)-1,1);
            Wod          = cell(1,N); % off-diagonal contributions
            updated      = false(1,numel(z));

            % Convert x to a cell
            if ~iscell(x)
                tmp = x;
                x = cell(1,numel(z));
                for n = 1:numel(z)
                    x{n} = reshape(tmp(offset(n):offset(n+1)-1),size(z{n}));
                end 
            end 
            % Expand to virtual variables (less efficient, but more readible)
            x = x(var2virt(var2virt <= numel(z)));
            if numel(this.C) > 0, z = [z this.C]; end                         
            z = z(var2virt);

            for n = find(~const)
                if widx(n)
                    z{n} = bsxfun(@times, weights{n}, z{n});
                    x{n} = bsxfun(@times, weights{n}, x{n});
                end
                XHZ = x{n}'*z{n};
                if uselincon, XHZ = P{n}'*XHZ*P{n}; end
                XHU(:,:,n) = freq(n) * conj(XHZ);
                Wod{n} = 0;
            end
            
            % Off-diagonal contributions
            for n = find(~const)
                Wn = zeros(R);
                for m = n+(freq(n)==1):N
                    if const(m), continue; end;
                    mult = freq; 
                    mult([n m]) = mult([n m]) - 1 - (m==n);
                    
                    if sum(mult) == 0
                        if m ~= n, Wn = Wn+XHU(:,:,m); end 
                        Wnm = P{m}*XHU(:,:,n)*P{m}';
                    else
                        if any(mult > 1), Wnm = UHU .^ (mult(:)*ones(1,R^2));
                        else Wnm = UHU(mult==1,:); end 
                        Wnm = reshape(prod(Wnm, 1), [R, R]);
                        if m ~= n, Wn = Wn+Wnm.*XHU(:,:,m); end 
                        Wnm = Wnm.*XHU(:,:,n);
                        if uselincon, Wnm = P{m}*Wnm*P{m}'; end 
                    end
                    Wod{m} = Wod{m} + (freq(m)-(n==m)) * Wnm;
                end
                if uselincon, Wn = P{n}*Wn*P{n}'; end
                Wod{n} = Wod{n} + freq(n) * Wn;
            end

            for n = find(~const)
                % Diagonal block contributions
                mult = [freq(1:n-1) freq(n)-1 freq(n+1:end)]; 
                if any(mult > 1), Wn = UHU .^ (mult(:)*ones(1,R^2));
                else Wn = UHU(mult==1,:); end 
                Wn = freq(n)*reshape(prod(Wn, 1), [R, R]);
                if uselincon, Wn = P{n}*Wn*P{n}'; end
                % Full contribution
                tmp = x{n}*Wn + z{n}*Wod{n};
                if widx(n), tmp = bsxfun(@times, weights{n}, tmp); end
                % Store in y
                idxn = offset(var2virt(n)):offset(var2virt(n)+1)-1;
                if ~updated(var2virt(n)), y(idxn) = tmp(:);
                else, y(idxn) = y(idxn) + tmp(:); end
                updated(var2virt(n)) = true;
            end 
        end

        function x = M_blockJacobi(this,~,b)
            N = this.N;
            R = this.R;
            offset = this.offset;
            var2virt = this.var2virt;
            var2mode = this.var2mode;
            
            if this.isincomplete
                weights = this.Winc;
                ind  = this.widx;
                weights(ind) = cellfun(@times, weights(ind), this.W(ind), 'UniformOutput', false);
            else 
                weights = this.W;
            end 
            
            x = zeros(size(b));
            for n = 1:N
                if this.isConstant(n), continue; end
                idx = offset(n):offset(n+1)-1;
                if this.useLDLPrec
                    tmp = ldldblksolve(this.invW{n}, b(idx));
                else 
                    tmp = reshape(b(idx),[],size(this.invW{n},1))*this.invW{n};
                    if this.widx(var2mode(n)) || this.isincomplete
                        w     = weights{var2mode(n)};
                        nzidx = abs(w.^2) > max(abs(w.^2))*numel(w)*eps(class(w))*10;
                        w     = 1./w.^2;
                        w(~nzidx) = 0;
                        tmp = bsxfun(@times, w, tmp);
                    end 
                end 
                x(idx) = tmp(:);
            end
        end
        
        function z = postaction(this, z)
        % Normalize the data 
            
            if this.useLinearConstraints, return; end
            if numel(z) <= 1, return; end
            
            ind = find(~this.isConstant);
            
            nrm    = cellfun(@(u) sqrt(sum(abs(u).^2,1)), z(ind), 'UniformOutput', false);
            z(ind) = cellfun(@(u,n) bsxfun(@rdivide, u, n), z(ind), nrm, 'UniformOutput', false);
            freq   = accumarray(this.var2virt(:), this.freq(:));
            freq   = freq(ind);
            
            N = getorder(this.T);
            if numel(this.variablesInMode) > N && ~this.isConstant(this.variablesInMode(N+1))
                nrm = prod(vertcat(nrm{:}).^(freq*ones(1,this.R)),1);
                idx = this.variablesInMode(getorder(this.T)+1);
                z{idx} = z{idx} .* nrm;
            else 
                nrm = nthroot(prod(vertcat(nrm{:}).^(freq*ones(1,this.R)),1), sum(freq));
                z(ind) = cellfun(@(u) bsxfun(@times, u, nrm), z(ind), 'UniformOutput', false);
            end 
        end         
        
        function isvalid = validate(this, z)
            import tensorlab.auxiliary.*
            
            isvalid = false;
            if ~isvalidtensor(this.T, false, this.type) 
                error('tok:invaliddata', 'The given tensor is not valid.');
            end 
            if ~iscellmat(z) 
                error('tok:invalidmodel', 'U should be a cell of matrices.');
            end 
            if ~isempty(this.C) && ~iscellmat(this.C)
                error('tok:invalidconstant', 'C should be a (possibly empty) cell of matrices.');
            end 

            if isempty(this.variablesInMode), vidx = 1:numel(z); 
            else vidx = this.variablesInMode; end
            if max(vidx) ~= numel(z)
                error('tok:invalidmodel', 'numel(z)=%d should equal max(variablesInMode)=%d.', ...
                      numel(z), max(vidx));
            end 
            if numel(this.C) > 0
                cidx = this.constantsInMode;
                if numel(cidx) ~= numel(vidx)
                    error('tok:invalidmodel', ...
                          'numel(variablesInMode) should equal numel(constantsInMode).');
                end
                if any((cidx == 0) & (vidx == 0)) || any((cidx ~= 0) & (vidx ~= 0))
                    error('tok:invalidmodel', ...
                          ['Either variablesInMode(n) or constantsInMode(n) should be nonzero, but ' ...
                           'not both, for all n.']);
                end 
                if ~isnonnegint(cidx(cidx>0), true) || max(cidx) > numel(this.C)
                    error('tok:cidx', ['constantsInMode should contain integers ' ...
                                       'between 1 and numel(C).']);
                end 
                sz2 = zeros(1,numel(vidx));
                sz2(vidx > 0) = cellfun('size', z(vidx(vidx>0)), 1);
                sz2(cidx > 0) = cellfun('size', this.C(cidx(cidx > 0)), 1);
            else 
                sz2 = getsize(z(vidx));
            end 
            if ~isnonnegint(vidx(vidx>0), true) || max(vidx) > numel(z)
                error('tok:vidx', ['variablesInMode should contain integers ' ...
                                   'between 1 and numel(z).']);
            end 
            sz1 = getsize(this.T);
            if ~sizeequals(sz1, sz2)
                error('tok:invalidmodel', ...
                      ['The dimensions of the model (%s) do not match the ' ...
                       'dimensions of the data (%s).'], size2str(sz2), size2str(sz1));
            end
            if ~isempty(this.P) 
                if numel(this.P) < numel(z)
                    this.P = [this.P(:).', repmat({1}, 1, numel(z)-numel(this.P))];
                end
                this.P(cellfun(@(p)isempty(p)||isscalar(p),this.P)) = {1};
                if ~iscellmat(this.P)
                    error('tok:invalidP', 'P should be a cell of scalars or matrices');
                end 
                if numel(z) ~= numel(this.P)
                    error('tok:invalidP', 'numel(P) should be equal to numel(z).');
                end 
                sz1 = cellfun('size', z(:).', 2);
                sz2 = cellfun('size', this.P(:).', 1);
                scal = cellfun(@isscalar, this.P(:).');
                if ~all( (sz1 == sz2) | scal )
                    error('tok:invalidP', 'size(P{n},1) should equal size(z{n},2) or be 1 for n=1,...,numel(z).');
                end 
            end
            if ~isempty(this.P)
                zext = cellfun(@(z,p) z*p, z, this.P, 'UniformOutput', false);
            else 
                zext = z;
            end
            if any(cellfun('size', zext, 2) ~= size(zext{1},2))
                if isempty(this.P) || all(cellfun(@isscalar, this.P(:).'))
                    msg = 'size(U{n},2) should be equal for all n.';
                    error('tok:invalidmodel', msg);
                else 
                    msg = 'size(U{n}*P{n},2)should be equal for all n.';
                    error('tok:invalidP', msg);                                    
                end 
            end 
            if ~isempty(this.C) && any(cellfun('size', this.C, 2) ~= size(zext{1},2))
                msg = 'size(Constants{n},2) should be equal to R for all n.';
                error('tok:invalidmodel', msg);
            end 
            if ~isempty(this.isConstant)
                isconst = this.isConstant;
                if numel(isconst) < numel(z)
                    this.isConstant =  [isconst(:).', zeros(1, numel(z)-numel(isconst))];
                elseif numel(isconst) > numel(z)
                    error('tok:invalidConstant', 'numel(isConstant) should equal numel(z)');
                end 
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
                newheader = ['Advanced CPD kernel routines including constants, symmetry and linear constraints.'];
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