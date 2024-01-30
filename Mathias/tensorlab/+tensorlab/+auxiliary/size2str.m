function s = size2str(sz, delim, brackets)
%SIZE2STR Convert size vector to string.
%   SIZE2STR(SZ) converts a size vector SZ to a string with 'x' between each
%   dimension. For example, SIZE2STR([2 3 4]) results in '[2x3x4]'.
%
%   SIZE2STR(SZ,DELIM) uses the DELIM as delimiter, which defaults to 'x'.
%
%   SIZE2STR(SZ,DELIM,BRACKETS) disables the brackets if BRACKETS = false.

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/02/08   NV      Initial version
    
    if nargin < 2 || isempty(delim), delim = 'x'; end
    if nargin < 3, brackets = true; end
    s = strrep(strrep(mat2str(sz(:).'), '  ', ' '), ' ', delim);
    if ~brackets, s = regexprep(s, '[\[\]]', ''); end
end
