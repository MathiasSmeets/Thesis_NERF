function [Z, output] = cpdli_minf(T, X, Z0, varargin)
%CPDLI_MINF CPD with linearly constrained factors by unconstrained minimization.
%   Z = CPDLI_MINF(T,X,Z0) computes the coefficients matrices Z{1}, Z{2}, ...,
%   Z{N} belonging to a canonical polyadic decomposition of the Nth-order
%   incomplete tensor T with factor matrices U{n} = X{n}*Z{n}, n = 1,...,N, by
%   minimizing 0.5*frob(cpdres(T, U))^2 in which only known entries in T are
%   taken into account. The algorithm is initialized with the coefficient
%   matrices Z0{n}.
%
%   Z = CPDLI_MINF(T,X,Z0, 'RegL2', lambda) can be used to add an addtional L2
%   regularization term 0.5*sum(lambda{n} .* Z{n}.^2) in which the sum goes
%   from n = 1,...,N. If lambda is a scalar, all variables have the same
%   regularization weight. If lambda is an array of length N, each variable
%   Z{n} has regularization weight lambda(n). If lambda is a cell of length
%   N, numel(lambda{n}) should be numel(Z{n}) and each entry in each variable
%   gets a separate weight. If numel(lambda{n}) is 0, 1 or size(X{n},2),
%   lambda{n} is automatically expanded to a matrix of size(Z{n}).  
%     
%   [Z,OUTPUT] = CPDLI_MINF(T,X,Z0) additionally returns the structure output
%   containing the following information:
%    
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   CPDLI_MINF(T,X,Z0,options), in which options is a struct, or
%   CPDLI_MINF(T,X,Z0,'key',value) may be used to set the following options:
%
%      options.RegL2 = []          - Add an L2 regularization term. (See above.)
%      options.UseDataDependent =  - Select the data-dependt or the
%      [true|false|'auto']           data-independent algorithm.
%      options.Algorithm =         - The desired optimization method.
%      [@minf_ncg|@minf_lbfgs|            
%       {@minf_lbfgsdl}|@minf_sr1cgs]       
%      options.<...>               - Parameters passed to the selected method,
%                                    e.g., options.TolFun, options.TolX and
%                                    options.CGMaxIter. See also help
%                                    [options.Algorithm].
%
%   See also cpdli_nls.

%   Authors: Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
%   Version History:
%   - 2017/04/03   NV      Initial version.
% 
%   References:
%   [1] Vervliet, N., Debals, O., De Lathauwer, L., "Canonical polyadic
%       decomposition of incomplete tensors with linearly constrained factors",
%       Technical Report 16-172, ESAT-STADIUS, KU Leuven, Leuven, Belgium, 2017.

    p = inputParser;
    p.addOptional('Algorithm', @minf_lbfgs);
    p.addOptional('MaxIter', 500);
    p.addOptional('TolFun', 1e-12);
    p.addOptional('TolX', 1e-8);
    p.addOptional('UseDataDependent', 'auto');
    p.addOptional('AllowChangeToDependent', true);
    p.addOptional('CheckInputs', true);
    p.addOptional('RegL2', []);
    p.addOptional('Display', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);

    % Check tensor
    if isnumeric(T), T = fmt(T); end
    if ~strcmpi(getstructure(T), 'incomplete')
        error('CPDLI_Kernel:T', 'T should be an incomplete tensor');
    end 
    if ~iscell(Z0)
        error('CPDLI_Kernel:Z', ...
              'Z should be a cell of coefficient matrices');
    end 
    
    % Data-dependence 
    if ischar(options.UseDataDependent) || isnan(options.UseDataDependent)
        % Select implementation based on complexity and memory if auto
        Nke = numel(T.val); 
        R   = size(Z0{1},2);
        N   = numel(Z0);
        szx = cellfun('size', Z0, 1);
        if prod(szx).^2 > 1e9/8 
            options.UseDataDependent = true;
        else 
            options.UseDataDependent = R*N^2*Nke < prod(szx)^2;
        end 
    end     

    % Construct kernel options
    fields = {'UseDataObjFun', 'UseDataGrad'};
    for f = fields
        if ~isfield(options, f),
            options.(f{1}) = options.UseDataDependent; 
        end
    end 

    % Construct kernel and check for errors
    newoptions = [fieldnames(options)'; struct2cell(options)'];
    kernel = CPDLI_Kernel(T, X, Z0, newoptions{:});
    if options.CheckInputs, kernel.validate(Z0); end 
    kernel.usePreconditioner = false;
    kernel.useGramian = false;
    
    % Initialize kernel
    kernel.initialize();
    
    % Construct optimization object
    fval = @kernel.objfun;
    dF = @kernel.grad;
    
    % Call optiomization routine
    [Z, output] = options.Algorithm(fval,dF,Z0(:).',newoptions{:});
    output.kernel = kernel;
    
end 