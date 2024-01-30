function [x,state] = struct_genvander(z,task,deg,mult)
%STRUCT_GENVANDER Generalized Vandermonde matrix.
%   [x,state] = struct_vander(z) uses the generator vector z to compute a
%   square Vandermonde matrix x as
%
%      [z(:).^0 z(:).^1 z(:).^2 ... z(:).^(length(z)-1)]
%
%   The structure state stores information which is reused in computing the
%   right and left Jacobian-vector products.
%
%   [x,state] = struct_vander(z,[],deg) generates x as the column vector 
%   z(:) raised to the powers deg(1):deg(2). The vector deg may contain
%   negative integers.
%
%   [x,state] = struct_vander(z,[],deg,transpose) transposes the Vandermonde
%   matrix if transpose is true. 
%
%   struct_vander(z,task,deg,transpose) computes the right or left
%   Jacobian-vector product of this transformation, depending on the structure
%   task. Use the structure state and add the field 'r' of the same shape as z
%   or the field 'l' of the same shape as x to obtain the structure task for
%   computing the right and left Jacobian-vector products
%   
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%   
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.
%   
%   See also struct_hankel, struct_toeplitz.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

if nargin < 2, task = []; end
if nargin < 3 || isempty(deg), deg = length(z)-1; end

sz_z = size(z);
z = z(:).';
if isempty(task) || (isempty(task.l) && isempty(task.r))
    d    = max(mult);
    n    = deg;

    pasc = zeros(n,d+1);
    
    % loop through each row
    baseVector = 1;
    for i = 1:n
        % store the actual base vector in the matrix
        pasc(i, 1:min(i,d+1)) = baseVector(1:min(length(baseVector),d+1));
        % update the base vector: generate the next row by convolution
        baseVector = conv(baseVector, [1 1]);
        baseVector = baseVector(1:min(length(baseVector),d+1));
    end
    %pasc = bsxfun(@rdivide, pasc, sqrt(n.^(0:d)));

    tmp  = bsxfun(@power, z, (-1:deg-1).');
    x    = zeros(deg, sum(mult+1));
    k    = 1;
    der  = zeros(deg, sum(mult+1));
    for i = 1:length(z)
        for d = 0:mult(i)
            x(d+1:end,k) = pasc(d+1:end,d+1) .* tmp(2:end-d,i);
            der(d+1:end,k) = pasc(d+1:end,d+1) .* tmp(1:end-d-1,i) .* (0:deg-d-1).';
            k = k + 1;
        end
    end
    state.deriv = der;
    ind = arrayfun(@(i) ones(mult(i)+1,1)*i, 1:length(mult), 'UniformOutput', false);
    state.P = sparse(1:sum(mult+1), vertcat(ind{:}),1);
elseif ~isempty(task.r)
    x = bsxfun(@times,task.deriv,task.r(:).'*task.P.');
    state = [];
elseif ~isempty(task.l)
    x = sum((conj(task.deriv).*task.l)*task.P,1);
    state = [];
end

end

function V = vander(g,deg)

if deg(1) > 0
    V = ones(size(g,1),deg(2)+1);
    V(:,2) = g;
    for i = 3:deg(2)+1, V(:,i) = V(:,i-1).*g; end
    V = V(:,deg(1)+1:deg(2)+1);
elseif deg(2) < 0
    V = ones(size(g,1),-deg(1)+1);
    V(:,end-1) = 1./g;
    for i = -deg(1)-1:-1:1, V(:,i) = V(:,i+1)./g; end
    V = V(:,1:-deg(1)+deg(2)+1);
else
    V = ones(size(g,1),deg(2)-deg(1)+1);
    for i = -deg(1):-1:1, V(:,i) = V(:,i+1)./g; end
    V(:,-deg(1)+2,:) = g;
    for i = -deg(1)+3:deg(2)-deg(1)+1, V(:,i) = V(:,i-1).*g; end
end

end
