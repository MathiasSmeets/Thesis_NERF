function res = vecequals(a,b)
%VECEQUALS Test if two vectors have identical entries
%   VECEQUALS(A,B) tests if two vectors A and B contain identical entries.

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/07/05   NV      Initial version
    
    res = (isempty(a) && isempty(b)) || ...
          (isvector(a) && isvector(b) && numel(a) == numel(b) && all(a(:) == b(:)));
end
