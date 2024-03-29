function [z,output] = minf_lbfgs(f,g,z0,varargin)
%MINF_LBFGS Minimize a function by L-BFGS with line search.
%   [z,output] = minf_lbfgs(f,g,z0) starts at z0 and attempts to find a
%   local minimizer of the real-valued function f(z). The input variables z
%   may be a scalar, vector, matrix, tensor or even a (nested) cell array
%   of tensors and its contents may be real or complex.
%
%   If f(x) is a function of real variables x, the function g(x) should
%   compute the partial derivatives of f with respect to the real variables
%   x, i.e. g(xk) := df(xk)/dx. If f(z) is a function of complex variables
%   z, the function g(z) should compute two times the partial derivative
%   of f with respect to conj(z) (treating z as constant), i.e. g(zk) :=
%   2*df(zk)/d(conj(z)) = 2*conj(df(zk)/dz). If g is the empty matrix [],
%   the real gradient or scaled conjugate cogradient is approximated with
%   finite differences. The output of the function g(z) may have the same
%   structure as z (although this is not necessary). The structure output
%   returns additional information:
%
%      output.alpha      - The line search step length in every iteration.
%      output.fevals     - The total number of function/gradient calls.
%      output.fval       - The value of the objective function f in every
%                          iteration.
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Step size tolerance reached.
%                             3: Maximum number of iterations reached.
%                             4: Absolute objective function tolerance
%                                reached.
%                             5: Measure tolerance reached. 
%      output.infols     - The circumstances under which the line search
%                          terminated in every iteration.
%      output.iterations - The number of iterations.
%      output.relfval    - The difference in objective function value
%                          between every two successive iterates, relative
%                          to its initial value.
%      output.relstep    - The step size relative to the norm of the 
%                          current iterate in every iteration.
%
%   minf_lbfgs(f,g,z0,options) may be used to set the following options:
%
%      options.Measure           - Function handle accepting the optimization 
%                                  variables as only input and returning a 
%                                  'measure' of goodness. (Smaller is better.)
%                                  E.g., @(z) frob(cpdres(T,z));
%                                  Default: @(z) nan;
%      options.MeasureIter       - Computes options.Measure very 
%      = options.Display           options.MeasureIter iterations. 
%      options.MeasureTol        - Tolerance for the measure, i.e., the algorithm 
%      = nan                       stops if options.Measure(z)<options.MeasureTol. 
%                                  Only computed if options.MeasureIter > 0. 
%      options.Display = 10      - Displays output information each
%                                  options.Display iterations. Set to 0 to
%                                  disable.
%      options.CustomDisplay     - Printer object used to print progress and
%      = ProgressPrinter()         termination messages. 
%      options.ShowCurves        - Show convergence curves in plot every 
%      = false                     options.Display iterations. Note that 
%                                  this can be slow.
%      options.LineSearch        - The line search used to minimize the
%      = @ls_mt                    objective function in the quasi-Newton
%                                  descent direction.
%      options.LineSearchOptions - The options structure passed to the line
%                                  search routine.
%      options.M                 - The number of L-BFGS updates to store.
%      = min(30,length(z0))
%      options.MaxIter = 500     - The maximum number of iterations.
%      options.TolFun = 1e-6     - The tolerance for output.relfval.
%      options.TolX = 1e-8       - The tolerance for output.relstep.
%      options.TolAbs = -inf     - The tolerance for output.fval.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables", SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.
%
% Version History:
% - 2018/04/09    NV    Added measure and tolerance, and custom printer
% - 2016/02/18    NV    Added absolute tolerance and option parser
    
% Evaluate the objective function and gradient.
dim = structure(z0);
fval = f(z0);
if ~isa(g,'function_handle') && isempty(g)
    grad = serialize(deriv(f,z0,fval));
else
    grad = serialize(g(z0));
end
z = z0;

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
p = inputParser;
p.addOptional('Display', 10);
p.addOptional('LineSearch', @ls_mt);
p.addOptional('LineSearchOptions', struct);
p.addOptional('M', min(30,length(grad))); 
p.addOptional('MaxIter', 500);
p.addOptional('TolFun', 1e-6);
p.addOptional('TolX', 1e-8);
p.addOptional('TolAbs', -inf);
p.addOptional('Measure', []);
p.addOptional('MeasureIter', nan);
p.addOptional('MeasureTol', nan);
p.addOptional('CustomDisplay', ProgressPrinter());
p.addOptional('ShowCurves', false);
p.KeepUnmatched = true;
p.parse(varargin{:});
options = p.Results;
if ~isfunc(options.LineSearch), options.LineSearch = @ls_mt; end
if ~isfield(options.LineSearchOptions,'alpha')
    options.LineSearchOptions.alpha = 1;
end
if ~isfield(options.LineSearchOptions,'c2')
    options.LineSearchOptions.c2 = 0.9;
end

if isempty(options.Measure)
    options.MeasureIter = nan;
elseif isnan(options.MeasureIter)
    options.MeasureIter = options.Display;
end
if options.MeasureIter == 0
    options.MeasureIter = nan;
end

% Initialize the algorithm.
S = zeros(numel(grad),options.M);
Y = zeros(numel(grad),options.M);
a = zeros(1,options.M);
r = zeros(1,options.M);
m = 0;
midx = [];

% L-BFGS with line search.
output.alpha = [];
output.fevals = 1;
output.fval = fval;
output.info = false;
output.infols = [];
output.iterations = 0;
output.relfval = [];
output.relstep = [];
output.measure = [];
if ~isnan(options.MeasureIter)
    output.measure(1) = options.Measure(z);
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

    % Compute the quasi-Newton step pqn = -H*grad.
    pqn = -grad;
    for i = 1:m
        a(i) = r(midx(i))*real(S(:,midx(i))'*pqn);
        pqn = pqn-a(i)*Y(:,midx(i));
    end
    if m > 0
        y1 = Y(:,midx(1));
        y1y1 = y1'*y1;
        gamma = y1y1*r(midx(1));
        pqn = 1/gamma*pqn;
    end
    for i = m:-1:1
        b = r(midx(i))*real(Y(:,midx(i))'*pqn);
        pqn = pqn+(a(i)-b)*S(:,midx(i));
    end
    
    % Minimize f along z+alpha*pqn.
    state = output; state.grad = grad;
    [alpha,outputls] = options.LineSearch( ...
        f,g,z,deserialize(pqn,dim),state,options.LineSearchOptions);
    output.alpha(:,end+1) = alpha;
    if length(alpha) < 2, alpha(2) = 1; end
    
    % Update iterate.
    z1 = serialize(z);
    z = alpha(2)*(z1+alpha(1)*pqn);
    if alpha(2) == 1
        s = alpha(1)*pqn;
    else
        s = z-z1;
    end
    z = deserialize(z,dim);
    
    % Update gradient and Hessian approximation.
    grad1 = grad;
    if isfield(outputls,'grad')
        grad = outputls.grad;
    else
        if ~isa(g,'function_handle') && isempty(g)
            grad = serialize(deriv(f,z,output.fval(end)));
        else
            grad = serialize(g(z));
        end
    end
    y = grad-grad1;
    sy = real(y'*s);
    if sy > 0
        m = min(m+1,options.M);
        if length(midx) < options.M, midx = m:-1:1;
        else midx = circshift(midx,[0 1]); end
        S(:,midx(1)) = s;
        Y(:,midx(1)) = y;
        r(:,midx(1)) = 1/sy;
    end
    
    % Update the output structure.
    if isfield(outputls,'fevals')
        output.fevals = output.fevals+outputls.fevals;
    end
    if isfield(outputls,'fval')
        output.fval(end+1) = outputls.fval;
    else
        output.fval(end+1) = f(z);
    end
    if isfield(outputls,'info')
        output.infols(end+1) = outputls.info;
    end
    output.iterations = output.iterations+1;
    output.relfval(end+1) = ...
        abs(diff(output.fval(end:-1:end-1)))/abs(output.fval(1));
    output.relstep(end+1) = norm(s)/norm(z1);
    if isnan(output.relstep(end)), output.relstep(end) = 0; end
    if output.relfval(end) <= options.TolFun, output.info = 1; end
    if output.relstep(end) <= options.TolX, output.info = 2; end
    if output.iterations >= options.MaxIter, output.info = 3; end
    if output.fval(end) < options.TolAbs, output.info = 4; end
    if options.MeasureIter > 0 && (mod(output.iterations, options.MeasureIter) == 0 || ...
                                   output.info > 0) 
        output.measure(end+1:output.iterations+1) = nan;
        output.measure(end) = options.Measure(z);
        if output.measure(end) < options.MeasureTol, 
            output.info = 5; 
        end
    end 
    
    % Display progress.
    printer.print(output);
end

% Display termination message.
printer.printTermination(output);

end

function [z,offset] = deserialize(z,dim,offset)
    if iscell(dim)
        v = z;
        z = cell(size(dim));
        if nargin < 3, offset = 0; end
        for i = 1:numel(z)
            if iscell(dim{i})
                [z{i},offset] = deserialize(v,dim{i},offset);
            else
                n = prod(dim{i}(:));
                z{i} = reshape(v(offset+(1:n)),dim{i});
                offset = offset+n;
            end
        end
    elseif ~isempty(dim)
        z = reshape(z,dim);
    end
end

function z = serialize(z)
    if iscell(z)
        for i = find(cellfun(@iscell,z(:).'))
            z{i} = serialize(z{i});
        end
        s = cellfun(@numel,z(:)); o = [0; cumsum(s)];
        c = z; z = zeros(o(end),1);
        for i = 1:length(s), z(o(i)+(1:s(i))) = c{i}(:); end
    else
        z = z(:);
    end
end

function dim = structure(z)
    if iscell(z)
        dim = cellfun(@size,z,'UniformOutput',false);
        for i = find(cellfun(@iscell,z(:).'))
            dim{i} = structure(z{i});
        end
    else
        dim = size(z);
        if numel(z) == dim(1), dim = []; end
    end
end
