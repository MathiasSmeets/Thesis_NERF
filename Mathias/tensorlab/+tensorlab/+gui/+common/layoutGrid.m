function layoutGrid(handles, offset, widths, vertsep, vcenter, hmargin)
    if nargin < 5, vcenter = false; end
    if nargin < 6, hmargin = 0; end
    if numel(hmargin) < numel(widths)
        hmargin = [hmargin ones(1, numel(widths)-numel(hmargin))*hmargin(end)];
    end 
    pos = get(handles, 'position');
    if iscell(pos)
        pos = cell2mat(reshape(pos,size(handles)));
        heights = max(pos(:,4:4:end),[],2);
    else 
        heights = pos(4);
    end 
    if nargin >= 4, heights = heights + vertsep; end
    for r = 1:size(handles,1)
        height = get(handles(r,:), 'pos');
        if iscell(height), height = vertcat(height{:}); end
        height = height(:,4);
        if vcenter
            loff = 0.5*(max(height) - height);
            loff(height == 0) = 0;
        else 
            loff = zeros(1,size(handles,2));
            height = ones(1,size(handles,2))*max(height);
        end 
        for c = 1:size(handles,2)
            pos = [offset(1) + sum(widths(1:c-1)), ...
                   offset(2) + loff(c) + sum(heights(r+1:end)), ...
                   widths(c)-hmargin(c), ...
                   height(c)];
            set(handles(r,c), 'Units', 'pixels', 'Position', pos);
        end 
    end 
end 