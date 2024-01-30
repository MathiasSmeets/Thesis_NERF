function [U, output] = cpdi_nls(T, U0, varargin)
%CPDI_NLS CPD of an incomplete tensor by nonlinear least squares.
%   Z = CPDI_NLS(T,U0) computes the factor matrices U{1}, U{2}, ..., U{N}
%   belonging to a canonical polyadic decomposition of the Nth-order incomplete
%   tensor T by minimizing 0.5*frob(cpdres(T, U))^2 in which only known entries
%   in T are taken into account. The algorithm is initialized with the
%   coefficient matrices U0{n}.
%     
%   [U,OUTPUT] = CPDI_NLS(T,U0) additionally returns the structure output
%   containing the following information:
%    
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   CPDI_NLS(T,U0,options), in which options is a struct, or
%   CPDI_NLS(T,U0,'key',value) may be used to set the following options:
%
%      options.Algorithm =         - The desired optimization method.
%      [@nls_gncgs| ...            
%       {@nls_gndl}|@nls_lm]       
%      options.LargeScale          - If true, the Gauss-Newton or Levenberg-
%      = [true|false|{'auto'}]       Marquardt steps are computed using a
%                                    preconditioned conjugate gradient algorithm.
%                                    Otherwise, a direct solver is used. When
%                                    'auto' is selected, the problem is
%                                    large-scale when the number of variables
%                                    exceeds 100, i.e., sum(cellfun(@numel,U0))
%                                    > 100.
%      options.M =                 - The preconditioner to use when
%      [{'block-Jacobi'}|...         options.LargeScale is true.
%       false]                     
%      options.<...>               - Parameters passed to the selected method,
%                                    e.g., options.TolFun, options.TolX and
%                                    options.CGMaxIter. See also help
%                                    [options.Algorithm].
%
%   See also cpdi_minf.

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
    p.addOptional('Algorithm', @nls_gndl);
    p.addOptional('MaxIter', 200);
    p.addOptional('TolFun', 1e-12);
    p.addOptional('TolX', 1e-8);
    p.addOptional('M', 'block-Jacobi');
    p.addOptional('CGMaxIter', 25);
    p.addOptional('LargeScale', 'auto');
    p.addOptional('CheckInputs', true);
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
        error('CPDI_Kernel:T', 'T should be an incomplete tensor');
    end 
    if ~iscell(U0)
        error('CPDI_Kernel:U', ...
              'U should be a cell of factor matrices');
    end 
    
    % Large scale test
    if ischar(options.LargeScale) || isnan(options.LargeScale) 
        options.LargeScale = sum(cellfun(@numel, U0)) > 100;
    end 
    
    % Construct kernel options
    if ~isfield(options, 'UsePreconditioner')
        options.UsePreconditioner = options.LargeScale && ...
            (~islogical(options.M) || options.M);
    end 

    % Construct kernel and check for errors
    newoptions = [fieldnames(options)'; struct2cell(options)'];
    kernel = CPDI_Kernel(T, U0, newoptions{:});
    if options.CheckInputs, kernel.validate(U0); end 
    kernel.useGramian = true;
    
    % Initialize kernel
    kernel.initialize(U0);
    
    % Construct optimization object
    fval = @kernel.objfun;
    dF = struct;
    dF.JHF = @kernel.grad;
    if options.LargeScale
        dF.JHJx = @kernel.JHJx;
    else 
        dF.JHJ = @kernel.JHJ;
    end
    if kernel.usePreconditioner
        dF.M = @kernel.M_blockJacobi;
    end 
    
    % Call optiomization routine
    [U, output] = options.Algorithm(fval,dF,U0(:).',newoptions{:});
    output.kernel = kernel;
end 