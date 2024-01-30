function res = iscellmat(x)
%ISCELLMAT Check if a cell of matrices.
%   ISCELLMAT(X) checks if X is a cell of (numerical) matrices (including scalars and vectors).
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/07/05   NV      Initial version
    
    res = iscell(x) && all(cellfun(@isnumeric, x)) && all(cellfun(@ismatrix, x));
end
