function res = cellequals(A,B)
%CELLEQUALS Check if two cells have identical content.
%   CELLEQUALS(A,B) checks if both A and B are cells with the same dimension and identical content.
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/07/05   NV      Initial version
    
    res = false;
    if ~iscell(A) || ~iscell(B), return; end
    if numel(A) ~= numel(B), return; end
    if ~all(cellfun(@numel, A) == cellfun(@numel, B)), return; end
    if ~all(cellfun(@ndims, A) == cellfun(@ndims, B)), return; end
    if ~all(cellfun(@(a,b) all(size(a)==size(b)), A, B)), return; end
    res = all(cellfun(@(a,b) all(a(:) == b(:)), A, B));
end
