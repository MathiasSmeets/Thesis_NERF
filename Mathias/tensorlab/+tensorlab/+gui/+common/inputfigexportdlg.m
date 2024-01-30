function varargout = inputfigexportdlg()
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
    
    % Create figure
    fh = figure();
    create();
    cbResize();
    % Wait for user
    varargout = {[],[]};
    uiwait(fh);
    try 
        close(fh);
    end
    
    function create()
        import tensorlab.gui.common.*;
        set(fh, 'WindowStyle', 'Modal', ...
                'NumberTitle', 'off', ...
                'Units', 'pixels', ...
                'Name', 'Export figure', ...
                'Visible', 'off');
        
        styles = getDefaultStyles();
        
        ph1 = uipanel(fh, styles.panel, 'Title', 'Export figure');
        uicontrol(ph1, styles.comment, ...
                  'Tag', 'lblcomment', ...
                  'String', 'Select file formats for this figure.', ...
                  'Position', [0 0 0 40]);

        ph1a = uipanel(ph1, styles.place, 'Tag', 'taboptions');

        extensions = {'png', 'fig', 'eps', 'pdf'};
        props(1,1) = uicontrol(ph1a, styles.label, 'String', 'File formats', ...
                               'Position', [0 0 0 15]);
        props(1,2) = uipanel(ph1a, styles.place, 'Tag', 'tabextensions', 'Position', [0 0 0 25]);
        for k = 1:numel(extensions)
            uicontrol(props(1,2), styles.checkbox, 'Tag', extensions{k}, ...
                      'String', extensions{k}, 'UserData', extensions{k}, ...
                      'Position', [0 0 0 25]);
        end 
        props(2,1) = uicontrol(ph1a, styles.label, 'String', 'Compress', ...
                               'Position', [0 0 0 15]);
        props(2,2) = uicontrol(ph1a, styles.checkbox, 'Tag', 'zip', ...
                               'String', 'collect in zip', 'UserData', 'zip', ...
                               'Position', [0 0 0 25]);
        layoutGrid(props, [10 10], [150 175], 5, true);
        
        ph2 = uipanel(fh, styles.panel, ...
                     'Title', 'Export', ...
                     'Position', [0 0 0 70]);
        btnok = uicontrol(ph2, styles.accentbutton, ...
                          'String', 'Export', ...
                          'Tag', 'btnexport', ...
                          'Callback', @cbExport, ...
                          'Position', [10 10 80 35], ...
                          'KeyReleaseFcn', @cbReturn);
        uicontrol(ph2, styles.button, ...
                  'String', 'Cancel', ...
                  'Callback', @cbCancel, ...
                  'Position', [10 10 80 35]);

        % Position figure elements
        pos = get(fh, 'Position');
        %% Set minimum size
        pos(3:4) = [350 170];
        set(fh, 'Position', pos);
        movegui(fh, 'center');
        set(fh, 'Visible', 'On');
        pause(0.1);
        drawnow();
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
        
        %% Position panel 1
        h = findall(fh, 'Tag', 'lblcomment');
        tmp = get(get(h, 'Parent'), 'Position');
        set(h, 'Position', [9 tmp(4)-45 tmp(3)-20, 20]);
        
        h = findall(fh, 'Tag', 'taboptions');
        tmp = get(get(h, 'Parent'), 'Position');
        set(h, 'Position', [10 tmp(4)-100 tmp(3)-20, 60]);
        chd = get(h, 'Children');
        chd = reshape(chd(end:-1:1), [2 2]).';
        layoutGrid(chd, [0 0], [110 tmp(3)-20-110-5], 5, true);

        h = findall(fh, 'Tag', 'tabextensions');
        tmp = get(h, 'Position');
        chd = get(h, 'Children');
        layoutGrid(chd(:).', [0 0], ones(1,numel(chd))*tmp(3)/numel(chd), 5, true);

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
    
    function cbExport(~,~)
        ext = get(findall(fh, 'Tag', 'tabextensions'), 'Children');
        val = cell2mat(get(ext, 'Value')) == 1;
        sel = get(ext(val), 'Userdata');
        zip = get(findall(fh, 'Tag', 'zip'), 'Value') == 1;

        varargout = {sel, zip};
        uiresume(fh);
    end 
    
    function cbCancel(~,~)
        uiresume(fh);
    end 
    
    function cbReturn(src,evt)
        if strcmpi(evt.Key, 'return'), cbExport(); end
    end 
    
end
