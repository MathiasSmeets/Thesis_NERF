function [z,output] = minf_sr1cgs(f,g,z0,varargin)
%MINF_SR1CGS Minimize a function by SR1 with CG-Steihaug.
%   [z,output] = minf_sr1cgs(f,g,z0) starts at z0 and attempts to find a
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
%      output.cginfo       - The circumstances under which CG-Steihaug
%                            terminated:
%                               1: Maximum number of iterations reached.
%                               2: Trust region radius reached.
%                               3: Direction of negative curvature found.
%      output.cgiterations - The number of CG iterations to solve the
%                            trust-region subproblem.
%      output.cgrelres     - The relative residual norm of the computed
%                            Steihaug step.
%      output.delta        - The trust region radius at every step attempt.
%      output.fevals       - The total number of function calls.
%      output.fval         - The value of the objective function f in every
%                            iteration.
%      output.gevals       - The total number of gradient calls.
%      output.info         - The circumstances under which the procedure
%                            terminated:
%                               1: Objective function tolerance reached.
%                               2: Step size tolerance reached.
%                               3: Maximum number of iterations reached.
%                               4: Absolute objective function tolerance
%                                  reached.
%                               5: Measure tolerance reached. 
%      output.iterations   - The number of iterations.
%      output.relfval      - The difference in objective function value
%                            between every two successive iterates,
%                            relative to its initial value.
%      output.relstep      - The step size relative to the norm of the 
%                            current iterate in every iteration.
%      output.rho          - The trustworthiness at every step attempt.
%
%   minf_sr1cgs(f,g,z0,options) may be used to set the following options:
%
%      options.CGMaxIter = 30 - The maximum number of CG iterations for
%                               solving the trust-region subproblem.
%      options.CGTol = 1e-6   - The tolerance for the CG method for solving
%                               the trust-region subproblem.
%      options.Measure        - Function handle accepting the optimization 
%                               variables as only input and returning a 
%                               'measure' of goodness. (Smaller is better.)
%                               E.g., @(z) frob(cpdres(T,z));
%                               Default: @(z) nan;
%      options.MeasureIter    - Computes options.Measure very 
%      = options.Display        options.MeasureIter iterations. 
%      options.MeasureTol     - Tolerance for the measure, i.e., the algorithm 
%      = nan                    stops if options.Measure(z)<options.MeasureTol. 
%                               Only computed if options.MeasureIter > 0. 
%      options.Delta =        - The initial trust region radius.
%      0.3*max(1,norm(z0))
%      options.Display = 10   - Displays the objective function value, its
%                               difference with the previous iterate
%                               relative to the first iterate and the
%                               relative step size each options.Display
%                               iterations. Set to 0 to disable.
%      options.CustomDisplay  - Printer object used to print progress and
%      = ProgressPrinter()      termination messages. 
%      options.ShowCurves     - Show convergence curves in plot every 
%      = false                  options.Display iterations. Note that 
%                               this can be slow.
%      options.M =            - The number of updates to store.
%      min(30,length(z0))
%      options.MaxIter = 500  - The maximum number of iterations.
%      options.TolFun = 1e-6  - The tolerance for output.relfval.
%      options.TolX = 1e-8    - The tolerance for output.relstep.
%      options.TolAbs = -inf  - The tolerance for output.fval.

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

% Store the structure of the input space and evaluate the gradient.
dim = structure(z0);
fval = f(z0);
if ~isa(g,'function_handle') && isempty(g)
    grad = serialize(deriv(f,z0,fval));
else
    grad = serialize(g(z0));
end
z0 = serialize(z0);

% Define the Hessian-vector product.
function y = Bx(x)
    y = x;
    if m > 0
        V = Y(:,midx)-B0*S(:,midx);
        y = B0*y+V*(D*real(V'*x));
    end
end

% Check the options structure.
p = inputParser;
p.addOptional('CGMaxIter', 30);
p.addOptional('CGTol', 1e-6);
p.addOptional('Delta', 0.3*max(1,norm(z0)));
p.addOptional('Display', 10);
p.addOptional('M', min(30,length(z0)));
p.addOptional('MaxIter', 500);
p.addOptional('PlaneSearch', false);
p.addOptional('PlaneSearchOptions', struct);
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

if isempty(options.Measure)
    options.MeasureIter = nan;
elseif isnan(options.MeasureIter)
    options.MeasureIter = options.Display;
end
if options.MeasureIter == 0
    options.MeasureIter = nan;
end

% Initialize the algorithm.
S = zeros(numel(z0),options.M);
Y = zeros(numel(z0),options.M);
a = zeros(1,options.M);
r = zeros(1,options.M);
m = 0;
midx = [];

% SR1 with CG-Steihaug.
output.cginfo = [];
output.cgiterations = [];
output.cgrelres = [];
output.delta = options.Delta;
output.fevals = 1;
output.fval = fval;
output.gevals = 1;
output.info = false;
output.iterations = 0;
output.relfval = [];
output.relstep = [];
output.rho = [];
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
printer.addField('delta', 'delta', '%-11s', '%10.4e');
printer.addField('rho', 'rho', '%-11s', '%10.4e');
printer.addField('cginfo', 'cginfo', '%-11s', '%-11s', '', '', false, ...
                 {'Res tol','Max iter', 'Max step len','Neg curv'});
printer.print(output, 'fval');

while ~output.info

    % Solve trust region subproblem with CG-Steihaug.
    rho = -inf;
    while rho <= 0

        % Compute the CG-Steihaug step p and estimate objective function
        % improvement.
        delta = output.delta(end);
        [p,output.cginfo(end+1),output.cgrelres(end+1), ...
            output.cgiterations(end+1)] ...
            = pcgsh(@Bx,-grad,delta,options.CGTol,options.CGMaxIter);
        dfval = -real(p'*grad)-0.5*real(p'*Bx(p));
        if isnan(output.delta(end))
            delta = max(1,norm(p));
            output.delta(end) = delta;
        end

        % Compute the trustworthiness rho.
        if dfval > 0
            z = deserialize(z0+p,dim);
            fval = f(z);
            rho = (output.fval(end)-fval)/dfval;
            if isnan(rho), rho = -inf; end
            output.rho(end+1) = rho;
            output.fevals = output.fevals+1;
        end

        % Update trust region radius delta.
        if rho > 0.5
            output.delta(end+1) = max(delta,2*norm(p));
        else
            sigma = (1-0.25)/(1+exp(-14*(rho-0.25)))+0.25;
            output.delta(end+1) = sigma*delta;
        end
        
        % Check for convergence.
        relstep = norm(p)/norm(z0); if isnan(relstep), relstep = 0; end
        if rho <= 0 && relstep <= options.TolX
            output.rho(end+1) = rho;
            fval = output.fval(end);
            z = deserialize(z0,dim);
            break;
        end

    end
    
    % Save current state.
    if rho > 0
        z0 = serialize(z);
        grad1 = grad;
    end

    % Evaluate the gradient and update step information.
    if rho > 0
        if ~isa(g,'function_handle') && isempty(g)
            grad = serialize(deriv(f,z,fval));
        else
            grad = serialize(g(z));
        end
        s = p;
        y = grad-grad1;
        ymBs = y-Bx(s);
        if abs(real(ymBs'*s)) > 1e-6*norm(s)*norm(ymBs)
            m = min(m+1,options.M);
            if length(midx) < options.M, midx = m:-1:1;
            else midx = circshift(midx,[0 1]); end
            S(:,midx(1)) = s;
            Y(:,midx(1)) = y;
            B0 = real(Y(:,midx(1))'*S(:,midx(1)))/ ...
                 (Y(:,midx(1))'*Y(:,midx(1)));
            L = real(Y(:,midx)'*S(:,midx));
            D = tril(L,-1)+diag(diag(L))+tril(L,-1)';
            D = D-B0*real(S(:,midx)'*S(:,midx));
            D = pinv(D);
        end
    end
    
    % Update the output structure.
    output.fval(end+1) = fval;
    output.gevals = output.gevals+1;
    output.iterations = output.iterations+1;
    output.relfval(end+1) = ...
        abs(diff(output.fval(end:-1:end-1)))/abs(output.fval(1));
    output.relstep(end+1) = relstep;
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

function [x,flag,relres,iter] = pcgsh(A,b,delta,tol,maxit,M)

% Check the options.
if nargin < 3, delta = nan; end
if nargin < 4 || isempty(tol), tol = 1e-6; end
if nargin < 5 || isempty(maxit), maxit = min(20,length(b)); end
PC = nargin > 5 && (isa(M,'function_handle') || ...
     (isnumeric(M) && all(size(M) == length(b))));

% Initialize PCG-Steihaug.
x = zeros(size(b));
r = -b;
if PC
    if isnumeric(M), y = M\r;
    else y = M(r); end
    d = -y;
    rr = r'*y;
else
    d = -r;
    rr = r'*r;
end
normb = sqrt(rr);
flag = 1;

% PCG-Steihaug.
for iter = 1:maxit

    if isnumeric(A), Ad = A*d;
    else Ad = A(d); end
    dAd = d'*Ad;
    x1 = x;
    
    if ~isnan(delta) && dAd <= 0
        xx = x1'*x1;
        dd = d'*d;
        c = real(x1'*d);
        alpha = (delta^2-xx)/(c+sqrt(c^2+dd*(delta^2-xx)));
        x = x1+alpha*d;
        flag = 3;
    end
    
    alpha = rr/dAd;
    x = x+alpha*d;
    
    if ~isnan(delta) && norm(x) >= delta
        xx = x1'*x1;
        dd = d'*d;
        c = real(x1'*d);
        alpha = (delta^2-xx)/(c+sqrt(c^2+dd*(delta^2-xx)));
        x = x1+alpha*d;
        flag = 2;
    end
    
    r = r+alpha*Ad;
    rr1 = rr;
    if PC
        if isnumeric(M), y = M\r;
        else y = M(r); end
        rr = r'*y;
    else
        rr = r'*r;
    end
    
    if PC, relres = norm(r)/normb;
    else relres = sqrt(rr)/normb; end
    if flag ~= 1, break; end
    if relres < tol, flag = 0; break; end

    beta = rr/rr1;
    if PC, d = -y+beta*d;
    else d = -r+beta*d; end

end

end
