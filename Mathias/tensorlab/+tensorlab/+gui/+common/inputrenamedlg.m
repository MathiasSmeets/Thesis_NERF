function varargout = inputrenamedlg(field, value, title)
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
                'Name', title, ...
                'Visible', 'off');
        
        styles = getDefaultStyles();
        
        ph1 = uipanel(fh, styles.panel, ...
                     'Title', title);
        uicontrol(ph1, styles.comment, ...
                  'Tag', 'lblcomment', ...
                  'String', sprintf('Set the new value for %s.', field), ...
                  'Position', [0 0 0 40]);

        props = [];
        props(1,1) = uicontrol(ph1, styles.label, 'String', 'Old value', ...
                               'Position', [0 0 0 15]);
        props(1,2) = uicontrol(ph1, styles.label, 'FontWeight', 'normal', ...
                               'String', value, 'Position', [0 0 0 15]);
        props(2,1) = uicontrol(ph1, styles.label, 'String', 'New value', ...
                               'Position', [0 0 0 15]);
        props(2,2) = uicontrol(ph1, styles.edit, 'Tag', 'txtnewvalue', ...
                               'String', value, 'Callback', @cbValidate, ...
                               'Position', [0 0 0 25]);
        
        layoutGrid(props, [10 10], [150 175], 5, true);
        
        ph2 = uipanel(fh, styles.panel, ...
                     'Title', 'Finalize', ...
                     'Position', [0 0 0 70]);
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

        % Position figure elements
        pos = get(fh, 'Position');
        %% Set minimum size
        pos(3:4) = [350 160];
        set(fh, 'Position', pos);
        movegui(fh, 'center');
        set(fh, 'Visible', 'On');
        uicontrol(btnok);
        setminsize(fh);
        set(fh, 'SizeChangedFcn', @cbResize);
    end
    
    function cbResize(~,~)
        import tensorlab.gui.common.*;

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
        %% Place fields in rename panel
        chd = get(panels(end), 'Children');
        tmp = get(panels(end), 'Position');
        % comment 
        set(chd(5), 'Position', [2*margin, tmp(4)-4*margin-20, tmp(3)-4*margin, 20]);
        % fields 
        layoutGrid(chd([4 3; 2 1]), [10 tmp(4)-3*margin-20-55], [150 tmp(3)-150-4*margin], 5, true);

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
        vals = get(findall(fh, 'Tag', 'txtnewvalue'), 'String');
        varargout = {vals};
        uiresume(fh);
    end 
    
    function cbCancel(~,~)
        uiresume(fh);
    end 

    function cbReturn(src,evt)
        if toc(btntimer) < 0.15, return; end
        if strcmpi(evt.Key, 'return'), cbOk(); end
    end 
    
    function cbValidate(src,evt)
        str = get(src, 'String');
        k = get(src,'Userdata');
        if isempty(get(src, 'String')), 
            errordlg(sprintf('The new value for %s cannot be empty', field));
        end
        btntimer = tic;
        uicontrol(btnok);
    end 
end
