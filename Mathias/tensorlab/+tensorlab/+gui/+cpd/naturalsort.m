function [res,ind] = naturalsort(str)
%NATURALSORT Sort strings using natural order. 
%   NATURALSORT(STRINGS) sorts a cell of strings using natural ordering, i.e., 'example10' comes
%   after 'example2'. 

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/18   NV      Initial version
    
    res = str;
    [numbers, str] = regexp(str, '[0-9]+', 'match', 'split');
    for k = 1:numel(str)
        num = cellfun(@(n) sprintf('%010d', str2num(n)), numbers{k}, 'UniformOutput', false);
        tmp = [reshape([str{k}(1:end-1); num], 1, []), str{k}(end)];
        str{k} = horzcat(tmp{:});
    end 
    [~,ind] = sort(str);
    res = res(ind);
end
