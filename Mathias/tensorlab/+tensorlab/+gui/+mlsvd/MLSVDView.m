classdef MLSVDView < tensorlab.gui.common.GuiView

    properties
        
    end
    
    methods
        
        function this = MLSVDView()
            this = this@tensorlab.gui.common.GuiView();
        end 
        
        function createLayout(this)
            import tensorlab.gui.common.*;
            
            styles = this.styles;
            fh = this.fh;
            set(fh, 'MenuBar', 'none');
            set(fh, 'ToolBar', 'none');
            set(fh, 'DockControls', 'off');
            set(fh, 'WindowStyle', 'normal');
            set(fh, 'Units', 'pixels');
            set(fh, 'NumberTitle', 'off', 'Name', 'MLSVD GUI');
            set(fh, 'Resize', 'off');
            set(fh, 'Position', [0 0 855 432]);
            
            %% Create menu bar
            mnu = uimenu(fh, 'Label', 'Menu', 'Tag', 'mnumain');
            sub = uimenu(mnu, 'Label', 'Load data', 'Tag', 'mnuhelp');
            uimenu(sub, 'Label', 'From file', 'Tag', 'mnuloadfile', 'Accelerator', 'F');
            uimenu(sub, 'Label', 'From workspace', 'Tag', 'mnuloadworkspace', 'Accelerator', 'W');
            uimenu(mnu, 'Label', 'Reset data', 'Tag', 'mnureset');
            uimenu(mnu, 'Label', 'Generate code', 'Tag', 'mnugeneratecode', 'Accelerator', 'O');
            uimenu(mnu, 'Label', 'Exit', 'Tag', 'mnuexit');
            mnu = uimenu(fh, 'Label', 'Help');
            uimenu(mnu, 'Label', 'About', 'Tag', 'mnuabout');
            uimenu(mnu, 'Label', 'Website', 'Tag', 'mnuwebsite');
            uimenu(mnu, 'Label', 'Citing', 'Tag', 'mnucite');
            uimenu(mnu, 'Label', 'Userguide', 'Tag', 'mnuserguide');

            %% Create main colums
            col1 = uipanel(fh, 'BorderType', 'none', 'Units', 'Pixels', ...
                           'Position', [0 0 450 0]);
            col2 = uipanel(fh, 'BorderType', 'none', 'Units', 'Pixels', ...
                           'Position', [0 0 400 0]);
            
            %% Data
            tabdata = uipanel(col1, styles.panel, 'Title', 'Data', ...
                              'Tag', 'tabdata', ...
                              'Position', [0 0 0 170]);
            props = [];
            props(1,:) = this.addRow(tabdata, 'load', 'Load data', 14);
            p = uipanel(tabdata, 'BorderType', 'none', 'Unit', 'pixels', ...
                        'Position', [0 0 300 30]);
            uicontrol(p, styles.button, 'String', 'From file...', ...
                      'Tag', 'btnfromfile', ...
                      'Position', [0 5 100 25]);
            uicontrol(p, styles.button, 'String', 'From workspace...', ...
                      'Tag', 'btnfromworkspace', ...
                      'Position', [105 5 150 25]);
            uicontrol(p, styles.button, 'String', 'Reset', ...
                      'Tag', 'btnreset', ...
                      'Position', [260 5 60 25]);
            props(1,2) = p;
            props(2,:) = this.addRow(tabdata, 'path', 'Path', 40);
            props(3,:) = this.addRow(tabdata, 'variable', 'Variable', 20);
            props(4,:) = this.addRow(tabdata, 'order', 'Order', 20);
            props(5,:) = this.addRow(tabdata, 'size', 'Size', 20);
            props(6,:) = this.addRow(tabdata, 'type', 'Type', 20);
            layoutGrid(props, [10 5], [100 330], 0, true);
            
            %% Inspect
            tabinspect = uipanel(col1, styles.panel, ...
                                 'Title', 'Exploration', ...
                                 'Position', [0 0 0 50]);
            uicontrol(tabinspect, styles.comment, ...
                      'String', 'Inspect a range of (approximate) multilinear singular values', ...
                      'Position', [10 5 330 20]);
                      
            h = uicontrol(tabinspect, styles.button, 'String', 'Inspect', ...
                          'Tag', 'btninspect', ...
                          'Position', [350 6 80 25]);
            
            %% Settings
            tabsettings = uipanel(col1, styles.panel, 'Title', 'Settings', ...
                                  'Tag', 'tabsettings', ...
                                  'Position', [0 0 0 190]);
            
            p = uipanel(tabsettings, 'BorderType', 'none', 'Units', 'Pixels', ...
                        'Position', [10 140 450 50]);
            uicontrol(p, styles.label, 'Tag', 'lblsizecore', 'String', 'Multilinear rank', ...
                      'Position', [0 0 200 20]);
            uicontrol(p, styles.edit, 'Tag', 'txtsizecore', 'String', '', ...
                      'HorizontalAlignment', 'center', 'Position', [215 0 175 25]);
            uicontrol(p, styles.button, 'Tag', 'btnhelprank', 'String', '?', ...
                      'Position', [395 1 25 25]);
            
            props = [];
            props(1,1) = uicontrol(tabsettings, styles.label, 'String', 'Algorithm');
            props(1,2) = uibuttongroup(tabsettings, styles.buttongroup, ...
                                       'Tag', 'bgalgorithm', ...
                                       'Position', [0 0 130 61]);
            props(1,3) = uicontrol(tabsettings, styles.label, 'String', 'Refinement');
            props(1,4) = uibuttongroup(tabsettings, styles.buttongroup, ...
                                       'Tag', 'bgrefinement', ...
                                       'Position', [0 0 130 61]);
            layoutGrid(props, [10 55], [100 115 100 100], 10);
            
            h = props(1,2);
            btns = [];
            btns(1,1) = uicontrol(h, styles.radio, 'Tag', 'rbmlsvdrsi', 'String', 'mlsvd_rsi');
            btns(2,1) = uicontrol(h, styles.radio, 'Tag', 'rblmlraaca', 'String', 'lmlra_aca');
            btns(3,1) = uicontrol(h, styles.radio, 'Tag', 'rbmlsvd', 'String', 'mlsvd');
            layoutGrid(btns, [0 5], 90);
            
            h = props(1,4);
            btns = [];
            btns(1,1) = uicontrol(h, styles.radio, 'Tag', 'rbnone', ...
                                  'String', 'none', 'UserData', 'none');
            btns(2,1) = uicontrol(h, styles.radio, 'Tag', 'rblmlraminf', ...
                                  'String', 'lmlra_minf', 'UserData', 'minf');
            btns(3,1) = uicontrol(h, styles.radio, 'Tag', 'rblmlranls', ...
                                  'String', 'lmlra_nls', 'UserData', 'nls');
            layoutGrid(btns, [0 5], 90);
            
            uicontrol(tabsettings, styles.accentbutton, 'Tag', 'btncompute', 'String', 'Run', ...
                      'Position', [0 10 150 35]);
            
            %% Results
            tabresults = uipanel(col2, styles.panel, 'Title', 'Results', ...
                                 'Tag', 'tabresults', ...
                                 'Position', [0 0 0 290]);
            
            uicontrol(tabresults, styles.comment, 'Tag', 'lblallentries', ...
                      'String', 'All entries', 'HorizontalAlignment', 'center', ...
                      'Position', [0 255 100 15]);

            p = uipanel(tabresults, 'BorderType', 'none', 'Units', 'Pixels', ...
                        'Position', [0 165 300 80]);
            props = [];
            props(1,:) = this.addRow(p, 'relerr', 'Relative error', 20);
            props(2,:) = this.addRow(p, 'relfit', 'Relative fit', 20);
            props(3,:) = this.addRow(p, 'abserr', 'Absolute error', 20);
            props(4,:) = this.addRow(p, 'ratio',  'Compression ratio', 20);
            layoutGrid(props, [10 0], [150 290]);

            this.drawLine(tabresults, [10 155], 365);
            uicontrol(tabresults, styles.comment, 'Tag', 'lblsingleentries', ...
                      'String', 'Single entries','HorizontalAlignment', 'center', ...
                      'Position', [0 135 100 15]);
                        
            p = uipanel(tabresults, 'BorderType', 'none', 'Units', 'Pixels', ...
                        'Position', [0 45 300 80]);
            props = [];
            props(1,:) = this.addRow(p, 'maxabserr', 'Max absolute error', 20);
            props(2,:) = this.addRow(p, 'maxabserrind', 'Index', 20);
            props(3,:) = this.addRow(p, 'maxrelerr', 'Max relative error', 20);
            props(4,:) = this.addRow(p, 'maxrelerrind',  'Index', 20);
            layoutGrid(props, [10 0], [150 290]);
            
            this.drawLine(tabresults, [10 40], 365);

            p = uipanel(tabresults, 'BorderType', 'none', 'Units', 'Pixels', ...
                        'Tag', 'tabresultbtns', 'Position', [0 10 285 25]);
            uicontrol(p, styles.button, 'Tag', 'btnvisualize', ...
                      'String', 'Visualize residual', ...
                      'Position', [0 0 140 25]);
            uicontrol(p, styles.button, 'Tag', 'btnerrslice', ...
                      'String', 'Plot slice errors', ...
                      'Position', [145 0 140 25]);
            
            %% Export
            tabexport = uipanel(col2, styles.panel, 'Title', 'Export to workspace', ...
                                'Tag', 'tabexport', ...
                                'Position', [0 0 0 125]);

            p = uipanel(tabexport, 'BorderType', 'none', 'Units', 'Pixels', ...
                        'Tag', 'tabexportnames', 'Position', [5 55 380 50]);
            w = 120; h = 20; m = 5;
            uicontrol(p, styles.label, 'Tag', 'lblexportU', 'String', 'Factors', ...
                      'Position', [m, h+m, w, h]);
            uicontrol(p, styles.edit, 'Tag', 'txtexportU', 'String', '', ...
                      'Position', [m, m, w, h]);
            uicontrol(p, styles.label, 'Tag', 'lblexportS', 'String', 'Core', ...
                      'Position', [w+2*m, h+m, w, h]);
            uicontrol(p, styles.edit, 'Tag', 'txtexportS', 'String', '', ...
                      'Position', [w+2*m, m, w, h]);
            uicontrol(p, styles.label, 'Tag', 'lblexportT', 'String', 'Tensor', ...
                      'Position', [2*w+3*m, h+m, w, h]);
            uicontrol(p, styles.edit, 'Tag', 'txtexportT', 'String', '', ...
                      'Position', [2*w+3*m, m, w, h]);
                        
            h = uicontrol(tabexport, styles.accentbutton, 'String', 'Export', ...
                          'Tag', 'btnexport', ...
                          'Position', [0 10 150 35]);
        end 
        
        function resize(this, ~, ~)
            pos = get(this.fh, 'Position');
            % layout columns
            cols = findall(this.fh, 'Type', 'uipanel', 'Parent', this.fh);
            cols = cols(end:-1:1);
            set(cols, 'Units', 'Pixels');
            widths = get(cols, 'Position');
            widths = cellfun(@(w) w(3), widths);
            margin = 5;
            offset = 5;
            for k = 1:numel(cols)
                tmp = [offset, margin, widths(k), pos(4)-2*margin];
                set(cols(k), 'Position', tmp);
                offset = offset + widths(k);
            end 
            % layout elements in column
            for k = 1:numel(cols)
                rows = get(cols(k),'Children');
                pos = get(cols(k), 'Position');
                if isempty(rows), continue; end
                set(rows, 'Units', 'Pixels');
                heights = get(rows, 'Position');
                if iscell(heights)
                    heights = cellfun(@(h) h(4), heights);
                else 
                    heights = heights(4);
                end 
                margin = 5;
                offset = margin;
                for l = 1:numel(rows)
                    tmp = [margin, offset, pos(3)-2*margin, heights(l)];
                    set(rows(l), 'Position', tmp);
                    offset = offset + heights(l) + margin;
                end 
            end 
            % speficics
            this.setPosition(findall(this.fh, 'Tag', 'lblallentries'), 'hcenter');
            this.setPosition(findall(this.fh, 'Tag', 'lblsingleentries'), 'hcenter');
            this.setPosition(findall(this.fh, 'Tag', 'tabresultbtns'), 'hcenter');
            
            this.setPosition(findall(this.fh, 'Tag', 'btncompute'), 'hcenter');
            %this.setPosition(findall(this.fh, 'Tag', 'btninspect'), 'hcenter');
            this.setPosition(findall(this.fh, 'Tag', 'btnexport'), 'hcenter');
        end 
        
        % function reset(this, ~, ~)
            
        % end 
        
    end 
    
end 
