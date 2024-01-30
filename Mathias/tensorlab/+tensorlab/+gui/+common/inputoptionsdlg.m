function varargout = inputoptionsdlg(name, options)
%INPUTOPTIONSDLG Input options dialog.
%   VALUES = INPUTOPTIONSDLG(NAME, FIELDS) creates an input dialog with a given NAME. VALUES is a
%   cell of values in the order defined by OPTIONS. If a user hits cancel, VALUES is empty. OPTIONS
%   is a struct array with the following fields:
%   
%       - name            unique name of option
%       - string          text to display
%       - tooltip         text shown as tooltip
%       - value           value to display
%       - default         default value to use when user hits reset
%       - validation      function handle accepting a single number of string and returns true if
%                         valid.
%       - error           error message if a value is not valid

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/08/19   NV      Initial version
    
    % Global field used to focus ok button
    btnok = []; 
    % Timer used to prevent pressing ok if pressing enter on a text field
    btntimer = tic;
    
    % Create figure
    fh = figure();
    create();
    cbResize();
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
                'Name', 'Set algorithm options', ...
                'Visible', 'off');
        
        styles = getDefaultStyles();
        
        ph1 = uipanel(fh, styles.panel, ...
                     'Title', sprintf('%s options', name));
        uicontrol(ph1, styles.comment, ...
                  'Tag', 'lblcomment', ...
                  'String', sprintf('Set options for the %s algorithm', name), ...
                  'Position', [0 0 0 40]);

        props = [];
        for k = 1:numel(options)
            props(k,1) = uicontrol(ph1, styles.label, 'String', options(k).string, ...
                                   'Position', [0 0 0 15]);
            props(k,2) = uicontrol(ph1, styles.edit, 'Tag', options(k).name, ...
                                   'String', options(k).value, 'UserData', k, ...
                                   'Callback', @cbValidate, 'Position', [0 0 0 25], ...
                                   'HorizontalAlignment', 'center');
            options(k).handle = props(k,2);
        end 
        layoutGrid(props, [10 10], [150 175], 5, true);
        
        ph2 = uipanel(fh, styles.panel, ...
                     'Position', [0 0 0 55]);
        btnok = uicontrol(ph2, styles.accentbutton, ...
                          'String', 'Ok', ...
                          'Tag', 'btnok', ...
                          'Callback', @cbOk, ...
                          'Position', [10 10 80 35], ...
                          'KeyReleaseFcn', @cbReturn);
        uicontrol(ph2, styles.button, ...
                  'String', 'Cancel', ...
                  'Callback', @cbCancel, ...
                  'Position', [10 10 80 35]);
        uicontrol(ph2, styles.button, ...
                  'String', 'Reset to default', ...
                  'Callback', @cbReset, ...
                  'Position', [10 10 140 35]);

        % Position figure elements
        pos = get(fh, 'Position');
        %% Set minimum size
        pos(3:4) = [350 30*numel(options)+70];
        set(fh, 'Position', pos);
        movegui(fh, 'center');
        set(fh, 'Visible', 'On');
        uicontrol(btnok);
        setminsize(fh);
        set(fh, 'SizeChangedFcn', @cbResize);
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
    
    function cbOk(~,~)
        vals = {options.value};
        varargout = {vals};
        uiresume(fh);
    end 
    
    function cbCancel(~,~)
        uiresume(fh);
    end 
    
    function cbReset(~,~)
        for k = 1:numel(options)
            set(options(k).handle, 'String', options(k).default);
        end
    end 
    
    function cbReturn(src,evt)
        if toc(btntimer) < 0.15, return; end
        if strcmpi(evt.Key, 'return'), cbOk(); end
    end 
    
    function cbValidate(src,evt)
        str = get(src, 'String');
        k = get(src,'Userdata');
        if isnumeric(options(k).default)
            num = str2num(str);
            if numel(num) ~= 1 || (isa(options(k).validation,'function_handle') && ...
                                  ~options(k).validation(num))
                errordlg(options(k).error);
                set(src, 'String', options(k).value);
            else 
                options(k).value = num;
            end 
        elseif ischar(options(k).default)
            if isa(options(k).validation,'function_handle') && ~options(k).validation(str)
                errordlg(options(k).error);
                set(src, 'String', options(k).value);
            else                 
                options(k).value = str;
            end 
        end 
        btntimer = tic;
        uicontrol(btnok);
    end 
end
