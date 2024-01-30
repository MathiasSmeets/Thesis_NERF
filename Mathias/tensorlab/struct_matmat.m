function [x,state] = struct_matmat(z,task,cnst)
%STRUCT_MATMAT Matrix-matrix multiplication.

%   [x,state] = struct_matmat(z,[]) computes x as z{1}*z{2} if numel(z) = 2 or
%   z{1}*z{2}*...*z{N} if numel(z) = N > 2. The structure state stores
%   information which is reused in computing the right and left Jacobian-vector
%   products.
%
%   struct_matmat(z,task) computes the right or left Jacobian-vector product of
%   this transformation, depending on the structure task. Use the structure
%   state and add the field 'r' of the same shape as z or the field 'l' of the
%   same shape as x to obtain the structure task for computing the right and
%   left Jacobian-vector products
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
%   See also struct_plus, struct_times, struct_matvec.

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%
% Version History:
% - 2017/02/09   NV      Extension to arbitrary number of variables.
% - 2017/02/08   NV      Initial version for two variables.

if nargin < 2, task = []; end
state = [];

right = isstruct(task) && isfield(task,'r') && ~isempty(task.r);
left  = isstruct(task) && isfield(task,'l') && ~isempty(task.l);

if ~left && ~right
    x = z{1};
    for k = 2:numel(z)
        x = x * z{k};
    end    
    % Cumulative products from the left
    state.lefts = cell(1,numel(z));
    tmp = 1;
    state.lefts{1} = tmp;
    for k = 2:numel(z)
        tmp = tmp * z{k-1}; 
        state.lefts{k} = tmp;
    end 
    % Cumulative products from the rights
    state.rights = cell(1,numel(z));
    tmp = 1;
    state.rights{end} = tmp;
    for k = numel(z)-1:-1:1
        tmp = z{k+1}*tmp;
        state.rights{k} = tmp;
    end     
elseif right
    x = 0;
    for k = 1:numel(z)
        x = x + task.lefts{k} * task.r{k} * task.rights{k};
    end   
elseif left
    x = cell(size(z));
    for k = 1:numel(z)
        x{k} = task.lefts{k}'*task.l*task.rights{k}';
    end 
end
