function [U,output] = cpd_core(T,U0,varargin)
%CPD_CORE Computational routines for CPD decomposition
%   CPD_CORE should not be called directly. Use CPD_MINF or CPD_NLS instead. 
%
%   See also cpd_minf, cpd_nls.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Optimization-based
%       algorithms for tensor decompositions: canonical polyadic
%       decomposition, decomposition in rank-(Lr,Lr,1) terms and a new
%       generalization," SIAM J. Opt., 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

    
    
    
    
    
% unstructuredtypes = {'full', 'incomplete', 'sparse'};
% type = getstructure(T);
% isstructured = ~any(strcmpi(type, unstructuredtypes));
% if ~isstructured, 
%     T = fmt(T,true); 
%     type = getstructure(T);
% end
% size_tens = getsize(T);

% isincomplete = strcmp('incomplete',type);
% issparse = strcmp('sparse',type);
    
% % Check the initial factor matrices U0.
% U0 = U0(:).';
% N = length(U0);
% R = size(U0{1},2);
% if any(cellfun('size',U0,2) ~= R)
%     error('cpd_core:U0','size(U0{n},2) should be the same for all n.');
% end
% if length(size_tens) > N
%     error('cpd_core:U0', 'length(U0) should be >= getorder(T)');
% end
% sz = cellfun('size', U0, 1);
% if any(sz(1:length(size_tens)) ~= size_tens & sz(1:length(size_tens))>0)
%     error('cpd_core:U0', 'size(U0{n},1) should be size(T,n) for all n');
% end
% if any(sz(length(size_tens)+1:end)~=1) 
%     error('cpd_core:U0', 'size(U0{n},1) should be 1 for all n > getorder(T)');
% end

% Check the options structure.
p = inputParser;
p.addOptional('OptimizationType', 'nls');
p.addOptional('M', nan);
p.addOptional('MaxIter', 200);
p.addOptional('CGMaxIter', 15);
p.addOptional('TolFun', 1e-12);
p.addOptional('TolX', 1e-8);
p.addOptional('TolAbs', 0);
p.addOptional('LargeScale', 'auto');
p.addOptional('TolLargeScale', 0.02);
p.addOptional('JHasFullRank', false);
p.addOptional('Display', 0);
p.addOptional('UseCPDI', false);
p.addOptional('FastUpdate', false);
p.addOptional('ForceStructured', true);
p.addOptional('LowerBound', {});
p.addOptional('UpperBound', {});
p.addOptional('VariablesInMode', []);
p.KeepUnmatched = true;
p.parse(varargin{:});

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

% Test type of tensor
T = fmt(T);
type = getstructure(T);

% Fast update if full tensor
kernelopt = struct;
isincomplete = strcmpi(type, 'incomplete') && options.UseCPDI;
if isincomplete
    kernel = CPDI_Kernel(T, kernelopt);
else 
    if strcmpi(type, 'full')
        if options.FastUpdate, kernelOpt.IsStructured = true; end
    end 
    if ~any(strcmpi(p.UsingDefaults, 'ForceStructured'))
        kernelopt.IsStructured = options.ForceStructured;
    end
    if ~isempty(options.VariablesInMode)
        kernelopt.VariablesInMode = options.VariablesInMode(:).';
    end
    kernel = CPD_Kernel(T,kernelopt);    
    if options.FastUpdate
        kernel.useFastUpdate = true;
    end 
end 


kernel.validate(U0);
try 
    kernel.initialize(U0);
catch e
    if strcmpi(e.identifier, 'frob:notImplemented');
        error('cpd_core:notImplemented', ...
              ['cpd_core does not support the structured tensor type %s, yet. Use ' ...
               'ful(T) instead.'], kernel.type);
    else 
        rethrow(e)
    end
end

% Line and plane search settings
isfunc = @(f)isa(f,'function_handle');
% Adapt line/plane search if it is a CPD line/plane search.
if isfield(options,'LineSearch') && isfunc(options.LineSearch) && ...
   ~isempty(strfind(func2str(options.LineSearch),'cpd_'))
    linesearch = options.LineSearch;
    options.LineSearch = @ls;
end
if isfield(options,'PlaneSearch') && isfunc(options.PlaneSearch) && ...
   ~isempty(strfind(func2str(options.PlaneSearch),'cpd_'))
    planesearch = options.PlaneSearch;
    options.PlaneSearch = @ps;
end
% Set absolute tolerances
if any(strcmpi(p.UsingDefaults, 'TolAbs')) 
    if ~isincomplete
        if kernel.isstructured
            options.TolAbs = 0.5*1e-15*kernel.T2;
        else 
            options.TolAbs = 0.5*options.TolAbs*kernel.T2;
        end
    end 
end

% Compute bounds
useBoundSolver = false;
lb = options.LowerBound;
if ~isempty(lb)
    if ~iscell(lb), 
        lb = repmat({lb},1,numel(U0)); 
    elseif numel(lb) ~= numel(U0)
        error('cpd_core:lb', 'numel(LowerBound) should be equal to numel(U0).')
    end
    for n = 1:numel(U0)
        if isempty(lb{n})
            lb{n} = -inf(size(U0{n}))
        elseif ~ismatrix(lb{n}) || ~all((size(lb{n}) == size(U0{n})) | size(lb{n}) == [1 1])
            error('cpd_core:lb', 'Dimension mismatch for LowerBound');
        else 
            lb{n} = bsxfun(@times, ones(size(U0{n})), lb{n});
        end 
    end 
    useBoundSolver = true;
end
ub = options.UpperBound;
if ~isempty(ub)
    if ~iscell(ub), 
        ub = repmat({ub},1,numel(U0)); 
    elseif numel(lb) ~= numel(U0)
        error('cpd_core:ub', 'numel(UpperBound) should be equal to numel(U0).')
    end
    for n = 1:numel(U0)
        if isempty(ub{n}), 
            ub{n} = inf(size(U0{n})); 
        elseif  ~ismatrix(ub{n}) || ~all((size(ub{n}) == size(U0{n})) | size(ub{n}) == [1 1])
            error('cpd_core:ub', 'Dimension mismatch for UpperBound');
        else 
            ub{n} = bsxfun(@times, ones(size(U0{n})), ub{n});
        end 
    end 
    useBoundSolver = true;
    if isempty(lb)
        ub = cellfun(@(u) -inf(size(u)), U0, 'UniformOutput', false);    
    end
elseif useBoundSolver
    ub = cellfun(@(u) inf(size(u)), U0, 'UniformOutput', false);    
end
if useBoundSolver %&& any(strcmpi(p.UsingDefaults, 'Algorithm'))
    options.Algorithm = @(F,dF,z0,options) nlsb_gndl(F,dF,lb,ub,z0,options);
end 
        
% Test large scale
if ischar(options.LargeScale)
    options.LargeScale = sum(cellfun(@numel, U0)) > 500;
end 

% Validate the algorithm
if isfield(options, 'Algorithm')
    options.Algorithm = validate_algorithm(options.OptimizationType, ...
                                           options.Algorithm);
else 
    options.Algorithm = validate_algorithm(options.OptimizationType);    
end

% Call the optimization method.
if strcmpi(options.OptimizationType, 'nls')
    kernel.useGramian = true;
    if options.LargeScale, dF.JHJx = @kernel.JHDJx; 
    else dF.JHJ = @kernel.JHDJ; end
   
    dF.JHF = @kernel.grad;
    if isnan(options.M), options.M = 'block-Jacobi'; end
    switch options.M
      case 'block-SSOR',   dF.M = @M_blockSSOR;
      case 'block-Jacobi', dF.M = @kernel.M_blockJacobi;
      case 'Jacobi',       dF.M = @M_Jacobi;
      otherwise, if isa(options.M,'function_handle'), dF.M = options.M; end
    end
    if isfield(dF, 'M'), kernel.usePreconditioner = true; end
else 
    dF = @kernel.grad;
    if isfield(options, 'M') && isnan(options.M)
        options = rmfield(options, 'M'); 
    end
end 
options.PostAction = @kernel.postaction;

[U,output] = options.Algorithm(@kernel.objfun,dF,U0(:).',options);
output.Name = func2str(options.Algorithm);

if output.info == 4 && ~isincomplete && kernel.isstructured
    warning('cpd_core:accuracy', ...
            ['Maximal numerical accuracy for structured tensors reached. The ' ...
             'result may be improved using ful(T) instead of T.']);
end

% % Computational core
% function state(z,firstrun)

%     if nargin == 2 && firstrun
%         % Store the fraction of known elements.
%         if isincomplete
%             cache.scale = length(T.val)./prod(T.size);
%         end
%         if options.UseCPDI
%             [cache.j, cache.i] = cellfun(@sort, T.sub, 'UniformOutput', 0);
%             cache.i = cellfun(@double, cache.i, 'UniformOutput', false);
%             cache.j = cellfun(@double, cache.j, 'UniformOutput', false);
            
%             % new
%             J = cell(1,N);
%             for n = 1:N
%                 i = cache.i{n};
%                 j = cache.j{n};
%                 idx = sortrows([i,j]);
%                 i = repmat(idx(:,1), 1, R);
%                 j = bsxfun(@plus, idx(:,2), size_tens(n)* ...
%                            ones(length(T.val),1)*(0:R-1));
%                 J{n} = sparse(j.',i.',1);
%             end
%             cache.J = vertcat(J{:});
%             % end new
            
%             cache.freq = cell(1,N);
%             for n = 1:N
%                 ntot = prod(size_tens([1:n-1 n+1:N]));
%                 cache.freq{n} = histc(T.sub{n}, 1:size_tens(n))/ntot;
%             end             
%         end
%     end

%     % Cache the factor matrices' Gramians.
%     cache.UHU = zeros(N,R*R);
%     for n = 1:N
%         tmp = conj(z{n}'*z{n});
%         cache.UHU(n,:) = tmp(:);
%     end
    
%     % Optionally cache the inverses of the Gramians for the preconditioner.
%     % In a faster language, this should be the Cholesky factor instead.
%     if ischar(options.M) || isa(options.M,'function_handle')
%         cache.invW = cell(1,N);
%         for n = 1:N
%             tmp = cache.UHU([1:n-1 n+1:N],:);
%             if N > 2, tmp = prod(tmp,1); end
%             cache.invW{n} = inv(reshape(tmp,[R R]));
%         end
%     end
    
%     if options.UseCPDI && (nargin < 2 || ~firstrun)
%         J = cell(R,length(z));
%         i = (1:numel(T.val)).';
%         for n = 1:length(z);
%             dims = [1:n-1 n+1:length(z)];
%             for r = 1:R
%                 val = z{dims(1)}(T.sub{dims(1)},r);
%                 for k = dims(2:end)
%                     if k <= N
%                         val = val .* z{k}(T.sub{k},r);
%                     else 
%                         val = val .* z{k}(1,r);
%                     end 
%                 end
%                 if n <= N
%                     tmp = sparse(double(T.sub{n}),i,val, size_tens(n), ...
%                                  numel(T.val));
%                 else 
%                     tmp = sparse(ones(size(T.sub{1})),i,val, 1, ...
%                                  numel(T.val));
%                 end 
%                 J{r,n} = tmp;                
%             end
%         end
%         cache.J = cat(1, J{:});
%     end 
% end

% function fval = objfun(z)
% % CPD objective function.
%     if isstructured
%         fval = abs(0.5*cache.T2 - real(inprod(T, z)) + 0.5*frob(z,'squared'));
%         isstructured = ~options.FastUpdate || abs(fval/cache.T2) > 1e-8;
%         if isstructured, return; end
%     end
    
%     if ~isincomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
%         fval = cpdgen(z); %z{1}*kr(z(end:-1:2)).';
        
%         if isincomplete, fval = fval(T.ind)-T.val;
%         elseif issparse
%             if ~isempty(T.ind), fval(T.ind) = fval(T.ind)-T.val; end
%         else
%             fval = fval-reshape(T,size(fval));
%         end
%     else
%         fval = -T.val;
%         for r = 1:R
%             tmp = z{1}(T.sub{1},r);
%             for n = 2:length(z), tmp = tmp.*z{n}(T.sub{n},r); end
%             fval = fval+tmp;
%         end
%     end
%     if isincomplete
%         E = T; E.val = fval;
%         if ~isempty(E.matrix)
%             E.matrix = sparse(double(E.sub{1}), ...
%                 double(1+idivide(E.ind-1,int64(size(E.matrix,1)))), ...
%                 double(fval),size(E.matrix,1),size(E.matrix,2));
%         end
%         cache.residual = E;
%     else
%         cache.residual = reshape(fval,size_tens);
%     end
%     fval = 0.5*(fval(:)'*fval(:));
    
% end

% function grad = grad(z)
%     if usestate, state(z); end
%     offset = cache.offset;    
%     grad = nan(offset(end)-1,1);
    
%     if isstructured
%         for n = 1:length(z)
%             tmp = mtkrprod(z,z,n) - mtkrprod(T,z,n);
%             grad(offset(n):offset(n+1)-1) = tmp(:);
%         end
%     else     
%         % CPDn scaled conjugate cogradient.
%         E = cache.residual;
%         if options.UseCPDI
%             grad = conj(cache.J)*E.val;
%             grad = grad(:);
%         else 
%             for n = 1:length(z)
%                 tmp = full(mtkrprod(E,z,n));
%                 grad(offset(n):offset(n+1)-1) = tmp(:);
%             end
%         end 
%     end
% end

% function JHJ = JHJ(z)
    
%     % CPD Jacobian's Gramian.
%     UHU = conj(cache.UHU);
%     JHJ = zeros(cache.offset(end)-1);
%     for n = 1:N
%         idxn = cache.offset(n):cache.offset(n+1)-1;
%         Wn = reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]);
%         JHJ(idxn,idxn) = kron(Wn,eye(size_tens(n)));
%         for m = n+1:N
%             idxm = cache.offset(m):cache.offset(m+1)-1;
%             Wnm = reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),[R R]);
%             JHJnm = bsxfun(@times,reshape(z{n},[size_tens(n) 1 1 R]), ...
%                     reshape(conj(z{m}),[1 size_tens(m) R 1]));
%             JHJnm = bsxfun(@times,JHJnm,reshape(Wnm,[1 1 R R]));
%             JHJnm = permute(JHJnm,[1 3 2 4]);
%             JHJnm = reshape(JHJnm,[size_tens(n)*R size_tens(m)*R]);
%             JHJ(idxn,idxm) = JHJnm;
%             JHJ(idxm,idxn) = JHJnm';
%         end
%     end
    
%     % If incomplete, approximate the effect of missing entries.
%     if isincomplete
%         JHJ = JHJ*cache.scale;
%     end
    
% end

% function JHJ = JHJi(~)
%     J = cache.J;
%     JHJ = conj(J*J');
% end

% function y = JHJx(z,x)
    
%     % CPD fast Jacobian's Gramian vector product.
%     % Ignores the fact that the tensor might be incomplete.
%     offset = cache.offset;
%     UHU = cache.UHU;
%     XHU = zeros(R,R,N);
%     y = nan(offset(end)-1,1);
%     for n = 1:N
%         Wn = UHU([1:n-1 n+1:N],:);
%         if N > 2, Wn = prod(Wn,1); end
%         tmp = reshape(x(offset(n):offset(n+1)-1),[],R);
%         XHU(:,:,n) = conj(tmp'*z{n});
%         y(offset(n):offset(n+1)-1) = tmp*reshape(Wn,[R R]);
%     end
%     for n = 1:N-1
%         idxn = offset(n):offset(n+1)-1;
%         Wn = zeros(R);
%         for m = n+1:N
%             idxm = offset(m):offset(m+1)-1;
%             if N == 2
%                 Wn = Wn+XHU(:,:,m);
%                 JHJmnx = z{m}*XHU(:,:,n);
%             else
%                 Wnm = UHU([1:n-1 n+1:m-1 m+1:N],:);
%                 if N > 3, Wnm = prod(Wnm,1); end
%                 Wnm = reshape(Wnm,[R R]);
%                 Wn = Wn+Wnm.*XHU(:,:,m);
%                 JHJmnx = z{m}*(Wnm.*XHU(:,:,n));
%             end
%             y(idxm) = y(idxm)+JHJmnx(:);
%         end
%         JHJnx = z{n}*Wn;
%         y(idxn) = y(idxn)+JHJnx(:);
%     end
    
%     % If incomplete, approximate the effect of missing entries.
%     if isincomplete
%         y = y*cache.scale;
%     end
    
% end

% function y = JHJxi(z,x)
%     J = cache.J;
%     y = x.'*J;
%     y = conj(J*y');
%     y = y(:);
% end

% function x = M_blockJacobi(~,b)

%     % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
%     % Equivalent to simultaneous ALS updates for each of the factors.
%     x = nan(size(b));
%     for n = 1:length(cache.offset)-1
%         idx = cache.offset(n):cache.offset(n+1)-1;
%         tmp = reshape(b(idx),[],size(cache.invW{1},1));
%         if options.UseNewPreconditioner
%             tmp = bsxfun(@rdivide, tmp, cache.freq{n});
%         end 
%         tmp = tmp*cache.invW{n};
%         x(idx) = tmp(:);
%     end
    
%     % If incomplete, approximate the effect of missing entries.
%     if isincomplete && ~options.UseNewPreconditioner
%         x = x/cache.scale;
%     end
    
% end

% function x = M_blockSSOR(z,b)
    
%     % Solve Mx = b, where M is a block-Symmetric Successive Overrelaxation
%     % preconditioner.
%     % x = inv(D)*(U+D)*b
%     B = cell(size(z));
%     UHU = cache.UHU;
%     BHU = nan(R,R,N);
%     for n = 1:N
%         B{n} = reshape(b(cache.offset(n):cache.offset(n+1)-1), size_tens(n), R);
%         BHU(:,:,n) = B{n}'*z{n};
%     end
%     X = B;
%     for n = 1:N-1
%         Wsum = zeros(R);
%         for m = n+1:N
%             Wnm = conj(reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),R,[]));
%             Wsum = Wsum+Wnm.*conj(BHU(:,:,m));
%         end
%         Wn = conj(reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]));
%         X{n} = X{n}+(z{n}*Wsum)/Wn;
%     end
%     % x = (L+D)*x
%     B = X;
%     for n = 1:N
%         Wn = conj(reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]));
%         BHU(:,:,n) = B{n}'*z{n};
%         X{n} = B{n}*Wn;
%     end
%     for n = 2:N
%         Wsum = zeros(R);
%         for m = 1:n-1
%             Wnm = conj(reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),R,[]));
%             Wsum = Wsum+Wnm.*conj(BHU(:,:,m));
%         end
%         X{n} = X{n}+z{n}*Wsum;
%     end
%     x = cell2mat(cellfun(@(x)x(:),X(:),'UniformOutput',false));
    
%     % If incomplete, approximate the effect of missing entries.
%     if isincomplete
%         x = x/cache.scale;
%     end
    
% end

function [alpha,output] = ls(~,~,z,p,state,options)
    state.UHU = cache.UHU;
    state.T2 = cache.T2;
    [alpha,output] = linesearch(T,z,p,state,options);
end

function [alpha,output] = ps(~,~,z,p,q,state,options)
    state.UHU = cache.UHU;
    state.T2 = cache.T2;
    [alpha,output] = planesearch(T,z,p,q,state,options);
end

end
