function isequal = sizeequals(sz1, sz2, strict)
%SIZEEQUALS Check if two size vectors are equal.
%   SIZEEQUALS(SZ1,SZ2) checks if two size vectors are equal regardless of
%   squeezed dimensions at the end.
%
%   SIZEEQUALS(SZ1,SZ2,TRUE) checks if the order and the dimensions are equal. 

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/02/08   NV      Initial version
    
    if nargin < 3, strict = false; end
    sz1 = sz1(:).'; 
    sz2 = sz2(:).';
    if strict && numel(sz1) ~= numel(sz2)
        isequal = false;
        return;
    end 
    if length(sz1) > length(sz2)
        sz2 = [sz2 ones(1, length(sz1)-length(sz2))];
    elseif length(sz1) < length(sz2)
        sz1 = [sz1 ones(1, length(sz2)-length(sz1))];
    end
    isequal = all(sz1 == sz2);
end
