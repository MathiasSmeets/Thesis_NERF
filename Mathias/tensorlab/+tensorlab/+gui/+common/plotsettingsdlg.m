function varargout = plotsettingsdlg(modesettings, termsettings, plottypes)
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
% - 2018/08/23   NV      Initial version
    
    % Global field used to focus ok button
    btnok = []; 
    % Timer used to prevent pressing ok if pressing enter on a text field
    btntimer = tic;
        
    if nargin < 3, plottypes = {'Line', 'Bar', 'Dot'}; end
    
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
                'Name', 'Set plot options', ...
                'Visible', 'off');
        
        styles = getDefaultStyles();
        
        %% Mode settings
        ph1 = uipanel(fh, styles.panel, 'Title', 'Settings for modes', ...
                      'Position', [0 0 0 numel(modesettings)*30+50]);

        props = [];
        props(1,1) = uicontrol(ph1, styles.label, 'Position', [0 0 0 15]);
        props(1,2) = uicontrol(ph1, styles.label, 'String', 'Name', 'Position', [0 0 0 15]);
        props(1,3) = uicontrol(ph1, styles.label, 'String', 'Plot type', 'Position', [0 0 0 15]);
        
        for k = 1:numel(modesettings)
            props(k+1,1) = uicontrol(ph1, styles.label, 'String', sprintf('Mode %d', k), ...
                                   'Position', [0 0 0 15]);
            props(k+1,2) = uicontrol(ph1, styles.edit, ...
                                     'String', modesettings(k).name, 'UserData', k, ...
                                     'Position', [0 0 0 25]);
            props(k+1,3) = uicontrol(ph1, styles.dropdown, ...
                                     'String', plottypes, 'UserData', k, ...
                                     'Position', [0 0 0 25]);
            set(props(k+1,3), 'Value', find(strcmpi(plottypes, modesettings(k).type)));
            modesettings(k).txthandle = props(k+1,2);
            modesettings(k).plothandle = props(k+1,3);
        end 
        layoutGrid(props, [10 10], [100 150 150], 5, true, [0 5 0]);
        
        %% Terms settings
        ph2 = uipanel(fh, styles.panel, 'Title', 'Settings for rank-1 terms', 'Position', [0 0 0 150]);

        uicontrol(ph2, styles.comment, 'Position', [0 0 0 15], 'Tag', 'lblr1tcomment', ...
                  'String', ['Set rank-1 term names according to template. %d is replaced by the ' ...
                            'number of the rank-1 term.']);

        p = uipanel(ph2, styles.place, 'Position', [0 0 400 30], 'Tag', 'tabtemplate');
        uicontrol(p, styles.label, 'String', 'Template', 'Position', [0 5 95 15]);
        uicontrol(p, styles.edit, 'String', 'Term %d', 'Tag', 'txttemplate', ...
                  'Position', [100 1 145 25]);
        uicontrol(p, styles.button, 'String', 'Apply', 'Callback', @cbApply, ...
                  'Position', [250 1 150 25]);
        
        p = uipanel(ph2, styles.place, 'Position', [0 0 400 30], 'Tag', 'tabscrollouter');
        p = uipanel(p, styles.place, 'Position', [0 0 300 numel(termsettings)*30], 'Tag', 'tabscrollinner');
        
        props = [];
        for k = 1:numel(termsettings)
            props(k,1) = uicontrol(p, styles.label, 'String', sprintf('Term %d', k), ...
                                   'Position', [0 0 0 15]);
            props(k,2) = uicontrol(p, styles.edit, 'String', termsettings(k).name, 'UserData', k, ...
                                   'Position', [0 0 0 25]);
            termsettings(k).handle = props(k,2);
        end 
        layoutGrid(props, [0 2], [100 150], 5, true);
        
        uicontrol(ph2, styles.slider, 'Tag', 'sldterms', 'Callback', @cbScroll);
                
        %% Finalize 
        ph3 = uipanel(fh, styles.panel, ...
                     'Position', [0 0 0 55]);
        btnok = uicontrol(ph3, styles.accentbutton, ...
                          'String', 'Ok', ...
                          'Tag', 'btnok', ...
                          'Callback', @cbOk, ...
                          'Position', [10 10 80 35], ...
                          'KeyReleaseFcn', @cbReturn);
        uicontrol(ph3, styles.button, ...
                  'String', 'Cancel', ...
                  'Callback', @cbCancel, ...
                  'Position', [10 10 80 35]);
        uicontrol(ph3, styles.button, ...
                  'String', 'Reset to default', ...
                  'Callback', @cbReset, ...
                  'Position', [10 10 140 35]);

        % Position figure elements
        pos = get(fh, 'Position');
        %% Set minimum size
        drawnow();
        pause(0.1);
        pos(3:4) = [425 30*(numel(modesettings)+1)+min(5,numel(termsettings))*30+170];
        set(fh, 'Position', pos);
        movegui(fh, 'center');
        set(fh, 'Visible', 'On');
        uicontrol(btnok);
        setminsize(fh);
        set(fh, 'SizeChangedFcn', @cbResize);
    end
    
    function cbScroll(~,~)
        panels = get(fh, 'Children');
        h = findall(panels(2), 'Tag', 'sldterms');
        val = get(h, 'Value');
        if round(val) ~= val, 
            set(h, 'Value', round(val)); 
            return;
        end
       
        h = findall(panels(2), 'Tag', 'tabscrollouter');
        tmp = get(h, 'Position');
        outerheight = tmp(4);
        h = findall(panels(2), 'Tag', 'tabscrollinner');

        nbrows = floor(outerheight / 30);
        nbterms = numel(termsettings);
        if nbrows >= nbterms, val = 1; 
        else, val = nbterms - nbrows - val + 1; end
        tmp = get(h, 'Position');
        offset = outerheight - tmp(4) - (1-val)*30 + 5;
        set(h, 'Position', [0, offset, tmp(3:4)]);
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
        heights(2) = pos(end) - sum(heights([1 3])) - (2+numel(panels))*margin;
        % Place panels at right positions
        for k = 1:numel(panels)
            tmp = [2*margin, sum(heights(1:k-1))+(1+k)*margin, pos(3)-3*margin, heights(k)];
            set(panels(k), 'Units', 'Pixels', 'Position', tmp);
        end 

        %% Do middle positions
        pos = get(panels(2), 'Position');
        h = findall(panels(2), 'Tag', 'lblr1tcomment');
        set(h, 'Position', [10, pos(4)-50, pos(3)-20, 30]);
        h = findall(panels(2), 'Tag', 'tabtemplate');
        set(h, 'Position', [10, pos(4)-80, pos(3)-20, 30]);
        delete(findall(panels(2), 'Type', 'Axes'));
        drawline(panels(2), [10 pos(4)-90], pos(3)-22);
        h = findall(panels(2), 'Tag', 'tabscrollouter');
        outerheight = floor((pos(4)-100)/30)*30;
        set(h, 'Position', [10, 5+mod(pos(4)-100,30), pos(3)-30, outerheight]);
        
        h = findall(panels(2), 'Tag', 'sldterms');
        set(h, 'Position', [pos(3)-30, 10, 20, pos(4)-105]);
        nbrows = floor(outerheight / 30);
        nbterms = numel(termsettings);
        if nbterms <= nbrows
            set(h, 'Visible', 'off');
        else 
            set(h, 'Visible', 'on');
            set(h, 'Min', 0, 'Max', nbterms-nbrows, 'Value', nbterms-nbrows, ...
                   'SliderStep', [1 1]./(nbterms-nbrows));
        end 
        
        cbScroll();
        
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
    
    function cbApply(~,~)
        str = get(findall(fh, 'Tag', 'txttemplate'), 'String');
        for k = 1:numel(termsettings)
            termsettings(k).name = sprintf(str, k*ones(1, sum(str == '%')));
            set(termsettings(k).handle, 'String', termsettings(k).name);
        end 
    end 
    
    function cbOk(~,~)
        for k = 1:numel(modesettings)
            modesettings(k).name = get(modesettings(k).txthandle, 'String');
            modesettings(k).type = plottypes{get(modesettings(k).plothandle, 'Value')};
        end
        for k = 1:numel(termsettings)
            termsettings(k).name = get(termsettings(k).handle, 'String');
        end
        modesettings = rmfield(modesettings, 'txthandle');
        modesettings = rmfield(modesettings, 'plothandle');
        termsettings = rmfield(termsettings, 'handle');
        varargout = {modesettings, termsettings};
        uiresume(fh);
    end 
    
    function cbCancel(~,~)
        uiresume(fh);
    end 
    
    function cbReset(~,~)
        for k = 1:numel(modesettings)
            set(modesettings(k).txthandle, 'String', sprintf('Mode %d', k));
            set(modesettings(k).plothandle, 'Value', 1);
        end
        set(findall(fh, 'Tag', 'txttemplate'), 'String', 'Term %d');
        cbApply();
    end 
    
    function cbReturn(src,evt)
        if toc(btntimer) < 0.15, return; end
        if strcmpi(evt.Key, 'return'), cbOk(); end
    end 

end
