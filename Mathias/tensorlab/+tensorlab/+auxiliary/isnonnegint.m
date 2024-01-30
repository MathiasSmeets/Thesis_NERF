function res = isnonnegint(x, strict)
%ISNONNEGINT Check if all entries are nonnegative integers.
%   ISNONNEGINT(X) checks if all entries are nonnegative integers.
%
%   ISNONNEGINT(X,true) checks if all entries are strictly nonnegative integers.

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/07/05   NV      Initial version
    
    res = all(round(x) == x);
    if nargin > 1 && strict 
        res = res & all(x > 0);
    else 
        res = res & all(x >= 0);
    end 
end
