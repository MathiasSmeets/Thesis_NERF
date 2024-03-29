function [U,output] = cpd_als(T,U0,varargin)
%CPD_ALS CPD by alternating least squares.
%   [U,OUTPUT] = CPD_ALS(T,U0) computes the factor matrices U{1}, ..., U{N}
%   belonging to a canonical polyadic decomposition of the N-th order
%   tensor T. The algorithm is initialized with the factor matrices U0{n}.
%   The structure output returns additional information:
%
%      output.alpha      - The value of the line or plane search step
%                          length(s) in every iteration.
%      output.fval       - The value of the objective function
%                          0.5*frob(T-cpdgen(U))^2 in every iteration.
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Step size tolerance reached.
%                             3: Maximum number of iterations reached.
%                             4: Absolute objective function tolerance
%                                reached.
%                             5: Measure tolerance reached. 
%      output.iterations - The number of iterations.
%      output.relfval    - The difference in objective function value
%                          between every two successive iterates, relative
%                          to its initial value.
%      output.relstep    - The step size relative to the norm of the 
%                          current iterate in every iteration.
%
%   CPD_ALS(T,U0,OPTIONS) or CPD_ALS(T,U0,'key',value) may be used to set the
%   following options:
%
%      options.Measure            - Function handle accepting the optimization 
%                                   variables as only input and returning a 
%                                   'measure' of goodness. (Smaller is better.)
%                                   E.g., @(z) frob(cpdres(T,z));
%                                   Default: @(z) nan;
%      options.MeasureIter        - Computes options.Measure very 
%      = options.Display            options.MeasureIter iterations. 
%      options.MeasureTol         - Tolerance for the measure, i.e., the algorithm 
%      = nan                        stops if options.Measure(z)<options.MeasureTol. 
%                                   Only computed if options.MeasureIter > 0. 
%      options.Display = 0        - Displays output information each
%                                   options.Display iterations.
%      options.CustomDisplay      - Printer object used to print progress and
%      = ProgressPrinter()          termination messages. 
%      options.ShowCurves         - Show convergence curves in plot every 
%      = false                      options.Display iterations. Note that 
%                                   this can be slow.
%      options.LineSearch =       - A function handle to the desired line
%      [{false}|@cpd_aels|...       search algorithm.
%       @cpd_els|@cpd_lsb]
%      options.LineSearchOptions  - An options structure passed to the
%                                   selected line search algorithm.
%      options.PlaneSearch =      - A function handle to the desired plane
%      [{false}|@cpd_eps]           search algorithm. Searches in the plane
%                                   spanned by the last two ALS updates.
%      options.PlaneSearchOptions - An options structure passed to the
%                                   selected plane search algorithm.
%      options.MaxIter = 500      - The maximum number of iterations.
%      options.Order = 1:N        - The order in which to update the factor
%                                   matrices. Must be a permutation of 1:N.
%      options.TolFun = 1e-8      - The tolerance for output.relfval. Note
%                                   that because the objective function is
%                                   a squared norm, TolFun can be as small
%                                   as eps^2.
%      options.TolX = 1e-6        - The tolerance for output.relstep.
%      options.TolAbs = -inf      - The tolerance for output.fval.
%      options.Delta = 1          - Step length restriction for update. If
%                                   Delta = 1, no restriction is used.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2018/04/09    NV    Added measure and tolerance, and custom printer
% - 2016/01/30    NV    Support for structured tensors

% Format the tensor T.
type = getstructure(T);
isstructured = ~any(strcmp(type, {'full', 'incomplete', 'sparse'}));
if ~isstructured, T = fmt(T,true); end

% Check the initial factor matrices U0.
U = U0(:).'; U1 = {};
R = size(U{1},2);
N = length(U);
if any(cellfun('size',U,2) ~= R)
    error('cpd_als:U0','size(U0{n},2) should be the same for all n.');
end
sz = cellfun('size', U, 1);
size_tens = getsize(T);
if length(sz) < length(size_tens)
    error('cpd_als:U0','length(U0) should >= getorder(T).');
end
if any(sz(1:length(size_tens)) ~= size_tens)
    error('cpd_als:U0','size(U0{n},1) should be size(T,n).');
end
if any(sz(length(size_tens)+1:end) ~= 1)
    error('cpd_als:U0','size(U0{n},1) should be 1 for n > getorder(T).');
end

p = inputParser;
p.addOptional('Display', 0);
p.addOptional('FastUpdate', true);
p.addOptional('LineSearch', false);
p.addOptional('LineSearchOptions', struct);
p.addOptional('PlaneSearch', false);
p.addOptional('PlaneSearchOptions', struct);
p.addOptional('MaxIter', 500);
p.addOptional('Order', 1:N);
p.addOptional('TolFun', 1e-8);
p.addOptional('TolX', 1e-6);
p.addOptional('TolAbs', -inf);
p.addOptional('Delta', nan);
p.addOptional('Measure', []);
p.addOptional('MeasureIter', nan);
p.addOptional('MeasureTol', nan);
p.addOptional('CustomDisplay', ProgressPrinter());
p.addOptional('ShowCurves', false);
p.KeepUnmatched = true;
p.parse(varargin{:});
options = p.Results;

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if ~xsfunc(options.LineSearch) ...
   && (~isa(options.LineSearch,'function_handle') && options.LineSearch)
    error('cpd_als:LineSearch','Not a valid line search algorithm.');
end
if ~xsfunc(options.PlaneSearch) ...
   && (~isa(options.PlaneSearch,'function_handle') && options.PlaneSearch)
    error('cpd_als:PlaneSearch','Not a valid line search algorithm.');
end
if any(strcmpi(p.UsingDefaults, 'TolAbs')) && isstructured
    options.TolAbs = 1e-15;
end
if isempty(options.Measure)
    options.MeasureIter = nan;
elseif isnan(options.MeasureIter)
    options.MeasureIter = options.Display;
end
if options.MeasureIter == 0
    options.MeasureIter = nan;
end

% Cache some intermediate variables.
try 
    T2 = frob(T,'squared');
catch e
    if strcmpi(e.identifier, 'frob:notImplemented');
        error('cpd_als:notImplemented', ...
              ['cpd_als does not support the structured tensor type %s, yet. Use ' ...
               'ful(T) instead.'], type);
    end
end
K = cell(1,N);
UHU = zeros(R,R,N);
for n = 1:N, UHU(:,:,n) = U{n}'*U{n}; end

% Alternating least squares.
first = options.Order(1);
last = options.Order(end);
K{first} = mtkrprod(T,U,first);
if isstructured,
    D = 0.5*T2 - real(inprod(T,U)) + 0.5*frob(U,'squared');
else 
    D = cpdres(T,U);
    D = 0.5*(D(:)'*D(:));
end 
output.alpha = [];
output.fval = D;
output.info = false;
output.iterations = 0;
output.relgain = [];
output.relfval = [];
output.relstep = [];
output.measure = [];
if ~isnan(options.MeasureIter)
    output.measure(1) = options.Measure(U);
end

% Progress printer
printer = options.CustomDisplay;
printer.Display = options.Display;
printer.ShowConvergenceCurves = options.ShowCurves;
printer.addField('fval',    'fval',    '%-15s', '%14.8e', '', '=1/2*norm(F)^2', true);
printer.addField('relfval', 'relfval', '%-15s', '%14.8e', options.TolFun, 'TolFun = %4.e', true);
printer.addField('relstep', 'relstep', '%-15s', '%14.8e', options.TolX, 'TolX = %4.e', true);
if ~isnan(options.MeasureIter)
    printer.addField('measure', 'measure', '%-11s', '%12.6e', options.MeasureTol, 'TolM = %4.e', true);
end 
printer.addField('alpha', 'alpha', '%-11s', '%10.4e');
printer.print(output, 'fval');

while ~output.info
    
    % Save current state.
    U2 = U1;
    U1 = U;
    
    % Update factor matrices.
    for n = options.Order
        W = prod(UHU(:,:,[1:n-1 n+1:N]),3);
        if n ~= first
            K{n} = mtkrprod(T,U,n);
        end
        if isnan(options.Delta) || options.Delta == 1
            U{n} = K{n}/conj(W);
        else 
            U{n} = (1-options.Delta)*U{n} + options.Delta*K{n}/conj(W);
        end
        UHU(:,:,n) = U{n}'*U{n};
    end
    
    % Line/plane search.
    [U,alpha,outputsrch] = search();
	
    % Update the output structure.
    output.alpha(:,end+1) = alpha;
    if isfield(outputsrch,'relgain')
        output.relgain(end+1) = outputsrch.relgain;
    end
    output.fval(end+1) = outputsrch.fval;
    output.iterations = output.iterations+1;
    output.relfval(end+1) = ...
        abs(diff(output.fval(end:-1:end-1)))/abs(output.fval(1));
    output.relstep(end+1) = sqrt(sum(cellfun(@(u,v)(u(:)-v(:))'* ...
                                             (u(:)-v(:)),U,U1)))/sqrt(sum(cellfun(@(u)u(:)'*u(:),U)));
    if output.relfval(end) <= options.TolFun,   output.info = 1; end
    if output.relstep(end) <= options.TolX,     output.info = 2; end
    if output.iterations >= options.MaxIter,    output.info = 3; end
    if 2*output.fval(end) <= options.TolAbs*T2, output.info = 4; end
    if options.MeasureIter > 0 && (mod(output.iterations, options.MeasureIter) == 0 || ...
                                   output.info > 0) 
        output.measure(end+1:output.iterations+1) = nan;
        output.measure(end) = options.Measure(U);
        if output.measure(end) < options.MeasureTol, 
            output.info = 5; 
        end
    end 
    
    % Display progress.
    printer.print(output);
end

% Display termination message.
printer.printTermination(output);

function [Ua,alpha,outputsrch,dU,dU1,dU2,D,n] = search()
% Line/plane search wrapper that guarantees the return values alpha and
% outputsrch.fval, and prepares the next ALS iteration.
    
    % Attempt line/plane search, and fill in missing alpha.
    outputsrch = struct;
    state = output; state.T2 = T2; state.UHU = UHU;
    if options.FastUpdate
        state.K1 = first; state.KN = last; state.K = K;
    end
    if isfunc(options.LineSearch)
        dU = cellfun(@(u,v)u-v,U,U1,'UniformOutput',false);
        [alpha,outputsrch] = options.LineSearch(T,U1,dU, ...
            state,options.LineSearchOptions);
        if any(isnan(alpha)) || isempty(alpha), alpha = [1 1];
        elseif length(alpha) == 1, alpha(2) = 1; end
    elseif isfunc(options.PlaneSearch)
        if output.iterations >= 1
            dU1 = cellfun(@(u,v)u-v,U,U1,'UniformOutput',false);
            dU2 = cellfun(@(u,v)u-v,U1,U2,'UniformOutput',false);
            [alpha,outputsrch] = options.PlaneSearch(T,U1,dU1,dU2, ...
                state,options.PlaneSearchOptions);
        else
            alpha = nan;
        end
        if any(isnan(alpha)) || isempty(alpha), alpha = [1 0 1];
        elseif length(alpha) == 2, alpha(3) = 1; end
    else
        alpha = 1;
    end
    
    % Fill in missing fval.
    if ~isfield(outputsrch,'fval')
        if isfunc(options.LineSearch)
            Ua = cellfun(@(u,v)alpha(2)*(u+alpha(1)*v),U1,dU, ...
                'UniformOutput',false);
        elseif isfunc(options.PlaneSearch) && output.iterations >= 1
            Ua = cellfun(@(u,v,w)alpha(3)*(u+alpha(1)*v+alpha(2)*w), ...
                U1,dU1,dU2,'UniformOutput',false);
        else
            Ua = U;
        end
        if options.FastUpdate
            K1 = mtkrprod(T,Ua,first);
            UHU1 = zeros(size(UHU));
            for n = 1:N, UHU1(:,:,n) = Ua{n}'*Ua{n}; end
        end
        if options.FastUpdate && log10(output.fval(end)) > log10(T2)-16+2.5
            outputsrch.fval = abs(.5*(T2+sum(sum(real(prod(UHU1,3)))))- ...
                            real(sum(dot(K1,Ua{first}))));
        else
            if isstructured
                D = 0.5*T2 - real(inprod(T,Ua)) + 0.5*frob(Ua,'squared');
                outputsrch.fval = D;
            else
                D = cpdres(T,Ua);
                outputsrch.fval = 0.5*(D(:)'*D(:));
            end
        end
        if outputsrch.fval > output.fval(end), K1 = []; end
    end
    
    % If line search did not improve the objective function, don't take it.
    if outputsrch.fval > output.fval(end)
        if isfunc(options.LineSearch)
            alpha = [1 1];
        elseif isfunc(options.PlaneSearch) && output.iterations >= 1
            alpha = [1 0 1];
        else
            alpha = 1;
        end
        if options.FastUpdate && log10(output.fval(end)) > log10(T2)-16+2.5
            outputsrch.fval = abs(0.5*(T2+sum(sum(real( ...
                W.*UHU(:,:,last)))))-real(sum(dot(K{last},U{last}))));
        else
            if isstructured
                D = 0.5*T2 - real(inprod(T,U)) + 0.5*frob(U,'squared');
                outputsrch.fval = D;
            else
                D = cpdres(T,U);
                outputsrch.fval = 0.5*(D(:)'*D(:));
            end
        end
    end
    
    % Prepare next ALS iteration, if necessary.
    if options.FastUpdate
        if exist('K1','var') && ~isempty(K1)
            K{first} = K1;
            UHU = UHU1;
        else
            if isfunc(options.LineSearch)
                Ua = cellfun(@(u,v)alpha(2)*(u+alpha(1)*v),U1,dU, ...
                    'UniformOutput',false);
            elseif isfunc(options.PlaneSearch) && output.iterations >= 1
                Ua = cellfun(@(u,v,w)alpha(3)*(u+alpha(1)*v+alpha(2)*w),...
                    U1,dU1,dU2,'UniformOutput',false);
            else
                Ua = U;
            end
            K{first} = mtkrprod(T,Ua,first);
            for n = 1:N, UHU(:,:,n) = Ua{n}'*Ua{n}; end
        end
    end
    
end

end
