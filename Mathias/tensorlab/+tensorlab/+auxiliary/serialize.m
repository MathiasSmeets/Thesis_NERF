function z = serialize(z)
%SERIALIZE Serialize a (nested) cell of arrays. 
%   SERIALIZE(Z) returns a vector which containing all entries in the given, potentially nested,
%   cell of arrays. 

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/07/06   NV      Initial version

    import tensorlab.auxiliary.serialize;
    if iscell(z)
        for i = 1:numel(z)
            if iscell(z{i}), z{i} = serialize(z{i}); end
            z{i} = z{i}(:);
        end
        z = vertcat(z{:});
    else
        z = z(:);
    end
end
