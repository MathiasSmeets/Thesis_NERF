function varargout = selectTensor(path, autoselect, varargin)
%SELECTTENSOR Select tensor from file or workspace.
%   NAME = SELECTTENSOR() creates a GUI to select a (valid) tensor from the workspace and returns
%   its NAME.
%
%   NAME = SELECTTENSOR(PATH) used the MAT-file pointed to by PATH to select a valid tensor. An
%   empty path selects the workspace. 
%  
%   NAME = SELECTTENSOR(PATH, AUTOSELECT) immediatly returns the variable name if only one
%   (valid) variable is found and AUTOSELECT is true. Default: AUTOSELECT = false.

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/08/09   NV      Initial version

        
    p = inputParser();
    p.addOptional('SelectionFcn', @(T) isvalidtensor(T) && getorder(T) >= 3);
    p.addOptional('Name', 'tensor');
    p.addOptional('NameArticle', []); 
    p.addOptional('NamePlural', []);
    p.addOptional('RankProperty', false);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    if isempty(options.NameArticle)
        if ~isempty(regexp(options.Name, '^[aeiou]'))
            options.NameArticle = 'an';
        else 
            options.NameArticle = 'a';
        end 
    end
    if isempty(options.NamePlural)
        options.NamePlural = [options.Name 's'];
    end
    
    if nargin < 1 || isempty(path)
        tensorsrc = 'Workspace';
        variables = populateFromWorkspace();
    else
        [~,file,ext] = fileparts(fullfile(path));
        tensorsrc = [file, ext];
        variables = populateFromFile(path);
    end 
    if nargin >= 2 && autoselect
        if numel(variables) == 1
            varargout = {variables{1}.varname};
        end 
    end 

    % Create figure
    fh = figure();
    create();
    cbResize();
    fillList(variables);
    % Wait for user
    varargout = {[]};
    uiwait(fh);
    try 
        close(fh);
    end
    
    function create()
        import tensorlab.gui.common.*;
        set(fh, 'WindowStyle', 'Modal', ...
                'NumberTitle', 'off', ...
                'Units', 'pixels', ...
                'Name', 'Select a tensor', ...
                'Visible', 'off');
        
        styles = getDefaultStyles();
        
        ph1 = uipanel(fh, styles.panel, ...
                     'Title', 'Select a variable');
        uicontrol(ph1, styles.comment, ...
                  'Tag', 'lblcomment', ...
                  'String', sprintf(['Select %s %s from %s.\r\nOnly valid %s are shown.'], ...
                                    options.NameArticle, options.Name, tensorsrc, options.NamePlural), ...
                  'Position', [0 0 0 40]);
        uicontrol(ph1, styles.list, ...
                  'Tag', 'lstvariables', ...
                  'Callback', @cbSelectionChanged, ...
                  'Position', [0 0 0 40]);
                  
        ph2 = uipanel(fh, styles.panel, ...
                     'Title', sprintf('Properties of selected %s', options.Name), ...
                      'Position', [0 0 0 110]);
        props = [];
        props(1,1) = uicontrol(ph2, styles.label, 'String', 'Variable');
        props(1,2) = uicontrol(ph2, styles.value, 'Tag',    'txtvariable');
        props(2,1) = uicontrol(ph2, styles.label, 'String', 'Order');
        props(2,2) = uicontrol(ph2, styles.value, 'Tag',    'txtorder');
        props(3,1) = uicontrol(ph2, styles.label, 'String', 'Size');
        props(3,2) = uicontrol(ph2, styles.value, 'Tag',    'txtsize');
        if options.RankProperty
            props(4,1) = uicontrol(ph2, styles.label, 'String', 'Rank');
            props(4,2) = uicontrol(ph2, styles.value, 'Tag',    'txtrank');
        else 
            props(4,1) = uicontrol(ph2, styles.label, 'String', 'Type');
            props(4,2) = uicontrol(ph2, styles.value, 'Tag',    'txttype');
        end 
        layoutGrid(props, [10 5], [150 180], 0);
        
        ph3 = uipanel(fh, styles.panel, ...
                     'Title', 'Finalize', ...
                     'Position', [0 0 0 75]);
        uicontrol(ph3, styles.accentbutton, ...
                  'String', 'Select', ...
                  'Tag', 'btnselect', ...
                  'Callback', @cbSelect, ...
                  'Position', [10 10 140 40]);
        uicontrol(ph3, styles.button, ...
                  'String', 'Cancel', ...
                  'Callback', @cbCancel, ...
                  'Position', [10 10 140 40]);

        % Position figure elements
        pos = get(fh, 'Position');
        %% Set minimum size
        pos(3:4) = [350 500];
        set(fh, 'Position', pos);
        movegui(fh, 'center');
        set(fh, 'Visible', 'On');
        setminsize(fh);
        set(fh, 'SizeChangedFcn', @cbResize);
    end

    function variables = populateFromWorkspace()
        variables = evalin('base', 'who');
        selfunname = 'selfun2315634';
        assignin('base', selfunname, options.SelectionFcn);
        for k = 1:numel(variables)
            try 
                if evalin('base', sprintf('%s(%s)', selfunname, variables{k}))
                    name = variables{k};
                    variables{k} = struct();
                    variables{k}.varname = name;
                    variables{k}.order   = evalin('base', sprintf('getorder(%s)', name));
                    variables{k}.size    = evalin('base', sprintf('getsize(%s)', name));
                    variables{k}.type    = evalin('base', sprintf('getstructure(%s)', name));
                    if options.RankProperty
                        variables{k}.rank = evalin('base', sprintf('size(%s{1},2)', name));
                    end 
                else  
                    variables{k} = [];
                end 
            catch
                variables{k} = [];
            end 
        end 
        try 
            evalin('base', sprintf('clear %s', selfunname));
        end 
        variables = variables(~cellfun(@isempty, variables));
    end 
    
    function variables = populateFromFile(path)
        data = load(path);
        variables = fieldnames(data);
        for k = 1:numel(variables)
            try 
                name = variables{k};
                if options.SelectionFcn(data.(name))
                    variables{k} = struct();
                    variables{k}.varname = name;
                    variables{k}.order   = getorder(data.(name));
                    variables{k}.size    = getsize(data.(name)); 
                    variables{k}.type    = getstructure(data.(name));
                    if options.RankProperty
                        tmp = data.(name);
                        variables{k}.rank = size(tmp, 2);
                    end 
                else  
                    variables{k} = [];
                end 
            catch 
                variables{k} = [];
            end 
        end 
        variables = variables(~cellfun(@isempty, variables));
    end 
    
    function fillList(variables)
        lst = findall(fh, 'Tag', 'lstvariables');
        set(lst, 'String', cellfun(@(s) s.varname, variables, 'UniformOutput', false));
        set(lst, 'Value', 1);
        cbSelectionChanged(lst);
    end 
    
    function cbResize(~,~)
        margin = 5; 
        pos = get(fh, 'Position');
        %% Compute heights for panels
        panels = get(fh, 'Children');
        if isempty(panels), return; end 
        heights = cell2mat(get(panels, 'Position'));
        heights = heights(:,4);
        % Top panel (= last of children) has a dynamic height
        heights(end) = pos(end) - sum(heights(1:end-1)) - (2+numel(panels))*margin;
        %% Place panels at right positions
        for k = 1:numel(panels)
            tmp = [2*margin, sum(heights(1:k-1))+(1+k)*margin, pos(3)-3*margin, heights(k)];
            set(panels(k), 'Units', 'Pixels', 'Position', tmp);
        end 
        %% Position elements in selection panel
        tmp = get(panels(3), 'Position');
        panelwidth = tmp(3);
        panelheight = tmp(4);
        h = findall(panels(3), 'Tag', 'lblcomment');
        tmp = get(h, 'Position');
        set(h, 'Position', [2*margin, panelheight-20-tmp(4), panelwidth-4*margin, tmp(4)]);
        h = findall(panels(3), 'Tag', 'lstvariables');
        set(h, 'Position', [2*margin, 2*margin, panelwidth-4*margin, panelheight-tmp(4)-20-2*margin]);
        
        %% Position buttons in finalize panel
        panelwidth = get(panels(1), 'Position');
        panelwidth = panelwidth(3);
        buttons = get(panels(1), 'Children');
        widths = cell2mat(get(buttons, 'Position'));
        widths = widths(:,3);
        totalwidth = sum(widths) + (numel(buttons)-1)*2*margin;
        xoffset = (panelwidth - totalwidth)/2;
        for k = numel(buttons):-1:1
            tmp = get(buttons(k), 'Position');
            tmp(1) = xoffset;
            set(buttons(k), 'Position', tmp);
            xoffset = xoffset + widths(k) + 2*margin;
        end 
    end 
    
    function cbSelectionChanged(src, evt)
        import tensorlab.auxiliary.size2str;
        btn = findall(fh, 'Tag', 'btnselect');
        if numel(variables) > 0 &&  src.Value > 0, state = 'on';
        else state = 'off'; end
        set(btn, 'Enable', state);
        
        if numel(variables) == 0 || src.Value <= 0, return; end
        
        var = variables{src.Value};
        set(findall(fh, 'Tag', 'txtvariable'), 'String', var.varname);
        set(findall(fh, 'Tag', 'txtorder'), 'String', var.order);
        set(findall(fh, 'Tag', 'txtsize'), 'String', size2str(var.size, 'x', false));
        if options.RankProperty
            set(findall(fh, 'Tag', 'txtrank'), 'String', var.rank);
        else 
            set(findall(fh, 'Tag', 'txttype'), 'String', var.type);
        end 
    end 
    
    function cbSelect(~,~)
        lst = findall(fh, 'Tag', 'lstvariables');
        idx = get(lst, 'Value');
        varargout = {variables{idx}.varname};
        uiresume(fh);
    end 
    
    function cbCancel(~,~)
        uiresume(fh);
    end 
    
end
