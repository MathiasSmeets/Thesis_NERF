classdef CPDView < tensorlab.gui.common.GuiView

    properties
        
    end
    
    methods
        
        function this = CPDView()
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
            set(fh, 'NumberTitle', 'off', 'Name', 'CPD GUI');
            set(fh, 'Position', [0 0 1110 585]);
            set(fh, 'Resize', 'off');
            
            %% Create menu bar
            mnu = uimenu(fh, 'Label', 'Menu', 'Tag', 'mnumain');
            sub = uimenu(mnu, 'Label', 'Load data', 'Tag', 'mnuhelp');
            uimenu(sub, 'Label', 'From file', 'Tag', 'mnuloadfile', 'Accelerator', 'F');
            uimenu(sub, 'Label', 'From workspace', 'Tag', 'mnuloadworkspace', 'Accelerator', 'W');
            uimenu(mnu, 'Label', 'Reset data', 'Tag', 'mnureset');
            uimenu(mnu, 'Label', 'Load session', 'Tag', 'mnuloadsession', 'Accelerator', 'L');
            uimenu(mnu, 'Label', 'Save session', 'Tag', 'mnusavesession', 'Accelerator', 'S');
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
                           'Position', [0 0 450 0]);
            col3 = uipanel(fh, 'BorderType', 'none', 'Units', 'Pixels', ...
                           'Position', [0 0 200 0]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Column 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            margin = 10;
            width  = 450 - 2*margin;
            
            %% Data
            tabdata = uipanel(col1, styles.panel, 'Title', 'Data', ...
                              'Tag', 'tabdata', ...
                              'Position', [0 0 0 200]);
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
            layoutGrid(props, [10 35], [100 325], 0, true);
            
            h = uicontrol(tabdata, styles.button, ...
                          'String', 'Inspect multilinear singular values', ...
                          'Tag', 'btninspect', ...
                          'Position', [0 5 300 25]);

            %% Initialization
            tabinit = uipanel(col1, styles.panel, 'Title', 'Initialization', ...
                              'Tag', 'tabinit', ...
                              'Position', [0 0 0 180]);
            h = uibuttongroup(tabinit, styles.buttongroup, ...
                              'Tag', 'bginitmethod', ...
                              'Position', [margin margin width 180]);
            uicontrol(h, styles.radio, 'Tag', 'rbmanual', 'String', 'Manual initialization', ...
                      'Position', [0 130 width 20]);
            uicontrol(h, styles.radio, 'Tag', 'rbauto', 'String', 'Automatic initialization', ...
                      'Position', [0 60 width 20]);
            
            p = uipanel(tabinit, 'BorderType', 'none', 'Units', 'Pixels', ...
                        'Tag', 'tabinitmanual', ...
                        'Position', [margin+30 margin+80 width-30 45]);
            props = [];
            props(1,1) = uicontrol(p, styles.label, 'String', 'Variable');
            props(1,2) = uicontrol(p, styles.value, 'String', 'T', ...
                                   'Tag', 'txtinittensor');
            props(1,3) = uicontrol(p, styles.button, ...
                          'String', 'Select', ...
                          'Tag', 'btninitselect', ...
                          'Position', [0 5 300 25]);
            props(2,1) = uicontrol(p, styles.label, 'String', 'Rank');
            props(2,2) = uicontrol(p, styles.value, 'String', '', ...
                                   'Tag', 'txtinitmanualrank');
            props(2,3) = uicontrol(p, styles.label);
            layoutGrid(props, [0 0], [120 width-250-margin 100], 0, true);
            
            
            p = uipanel(tabinit, 'BorderType', 'none', 'Units', 'Pixels', ...
                        'Tag', 'tabinitauto', ...
                        'Position', [margin+30 margin width-30 55]);
            uicontrol(p, styles.label, 'Tag', 'lblrank', 'String', 'Rank', ...
                      'Position', [0 30 200 20]);
            uicontrol(p, styles.edit, 'Tag', 'txtrank', 'String', '', ...
                      'HorizontalAlignment', 'center', 'Position', [210 30 150 25]);
            uicontrol(p, styles.button, 'Tag', 'btnhelprank', 'String', '?', ...
                      'Position', [365 31 25 25]);
            
            props = [];
            props(1,1) = uicontrol(p, styles.label, 'String', 'Algorithm');
            props(1,2) = uibuttongroup(p, styles.buttongroup, ...
                                       'Tag', 'bginitalgorithm', ...
                                       'Position', [0 0 width 10]);
            layoutGrid(props, [0 0], [120 width-120], 10);
            
            h = props(1,2);
            btns = [];
            btns(1,1) = uicontrol(h, styles.radio, 'Tag', 'rbgevd', 'String', 'GEVD');
            btns(1,2) = uicontrol(h, styles.radio, 'Tag', 'rbrandn', 'String', 'randn');
            btns(1,3) = uicontrol(h, styles.radio, 'Tag', 'rbrand', 'String', 'rand');
            layoutGrid(btns, [0 5], [90 90 90]);
            
            %% Settings
            tabsettings = uipanel(col1, styles.panel, 'Title', 'Settings and constraints', ...
                                  'Tag', 'tabsettings', ...
                                  'Position', [0 0 0 180]);

            p = uipanel(tabsettings, 'BorderType', 'none', 'Units', 'Pixels', ...
                        'Tag', 'tabalgorithm', ...
                        'Position', [margin margin+130 width 55]);
            
            uicontrol(p, styles.label, 'String', 'Algorithm',...
                      'Position', [0 4 100 15]);
            bg = uibuttongroup(p, styles.buttongroup, ...
                               'Tag', 'bgalgorithm', ...
                               'Position', [150 4 200 25]);
            uicontrol(p, styles.button, ...
                      'String', 'Options', ...
                      'Tag', 'btnoptions', ...
                      'Position', [width-100-margin 0 100 25]);
            
            btns = [];
            btns(1,1) = uicontrol(bg, styles.radio, 'Tag', 'rbnls', 'String', 'NLS');
            btns(1,2) = uicontrol(bg, styles.radio, 'Tag', 'rbals', 'String', 'ALS');
            layoutGrid(btns, [0 0], [90 90]);

            props = [];
            props(1,:) = this.addRow(tabsettings, 'symmetry', 'Symmetry in modes', [15 25], true);
            props(2,:) = this.addRow(tabsettings, 'nonneg', 'Nonnegative factors', [15 25], true);
            layoutGrid(props, [margin 80], [150 width-160], 5, true);

            uicontrol(tabsettings, styles.checkbox, 'Tag', 'chkplotcvg', ...
                      'String', 'Plot convergence', ...
                      'Position', [margin 55 width 20]);

            uicontrol(tabsettings, styles.accentbutton, 'Tag', 'btncompute', 'String', 'Run', ...
                      'Position', [0 10 150 35]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Column 2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            margin = 10;
            width  = 450 - 2*margin;
            
            %% Model
            tabmodel = uipanel(col2, styles.panel, 'Title', 'Model', ...
                               'Tag', 'tabmodel', ...
                               'Position', [0 0 0 180]);

            uipanel(tabmodel, styles.place, 'Tag', 'tabfactors', ...
                    'Position', [margin margin width-margin 150]);
            
            %% Results
            tabresults = uipanel(col2, styles.panel, 'Title', 'Results', ...
                                 'Tag', 'tabresults', ...
                                 'Position', [0 0 0 325]);
            
            height = 27;

            offset = 315;
            offset = offset - 55;
            p = uipanel(tabresults, styles.place, 'Position', [margin offset width 55]);
            uicontrol(p, styles.label, 'String', 'Stop criterion', ...
                      'Position', [0  height 150 15]);
            uicontrol(p, styles.value, 'Tag', 'txtstopcriterion', ...
                      'String', 'Max iterations', ...
                      'Position', [155 height 260 15]);
            uicontrol(p, styles.button, 'Tag', 'btnshowcurves', ...
                      'String', 'Show curves', ...
                      'Position', [155  0 265 25]);
            
            offset = offset - height;
            p = uipanel(tabresults, styles.place, 'Position', [margin offset width height]);
            uicontrol(p, styles.label, 'String', 'Residual', ...
                      'Position', [0   4 150 15]);
            uicontrol(p, styles.button, 'Tag', 'btnvisualize', ...
                      'String', 'Visualize', ...
                      'Position', [155 0 130 25]);
            uicontrol(p, styles.button, 'Tag', 'btnerrslice', ...
                      'String', 'Plot slice errors', ...
                      'Position', [290 0 130 25]);

            offset = offset - height;
            p = uipanel(tabresults, styles.place, 'Position', [margin offset width height]);
            uicontrol(p, styles.label, 'String', 'Relative error', ...
                      'Position', [0   4 150 15]);
            uicontrol(p, styles.value, 'Tag', 'txtrelerr', ...
                      'Position', [155 4 130 15]);

            offset = offset - height;
            p = uipanel(tabresults, styles.place, 'Position', [margin offset width height]);
            uicontrol(p, styles.label, 'String', 'Corcondia', ...
                      'Position', [0   4 150 15]);
            uicontrol(p, styles.value, 'Tag', 'txtcorcondia', ...
                      'Position', [155 4 130 15]);
            
            offset = offset - height;
            p = uipanel(tabresults, styles.place, 'Position', [margin offset width height]);
            uicontrol(p, styles.label, 'String', 'Condition number', ...
                      'Position', [0   4 150 15]);
            uicontrol(p, styles.value, 'Tag', 'txtcondition', ...
                      'Position', [155 4 130 15]);
            uicontrol(p, styles.button, 'Tag', 'btncomputecondition', ...
                      'String', 'Compute', 'Visible', 'off', ...
                      'Position', [290 0 130 25]);

            offset = offset - 45;
            p = uipanel(tabresults, styles.place, 'Position', [margin offset width 40]);
            uicontrol(p, styles.label, 'String', 'Ratio max/min singular value', ...
                      'Position', [0   0 140 35]);
            uicontrol(p, styles.value, 'Tag', 'txtratiosv', ...
                      'Position', [155 1 130 35]);
            
            offset = offset - height;
            p = uipanel(tabresults, styles.place, 'Position', [margin offset width height]);
            uicontrol(p, styles.label, 'String', 'Norm rank-1 terms', ...
                      'Position', [0   4 150 15]);
            uicontrol(p, styles.value, 'Tag', 'txtnorms', ...
                      'Position', [155 4 210 15]);
            uicontrol(p, styles.button, 'Tag', 'btnshownorms', ...
                      'String', 'Plot', ...
                      'Position', [370 0 50 25]);

            offset = offset - height;
            p = uipanel(tabresults, styles.place, 'Position', [margin offset width height]);
            uicontrol(p, styles.label, 'String', 'Angles rank-1 terms', ...
                      'Position', [0   4 150 15]);
            uicontrol(p, styles.value, 'Tag', 'txtangles', ...
                      'Position', [155 4 210 15]);
            uicontrol(p, styles.button, 'Tag', 'btnshowangles', ...
                      'String', 'Plot', ...
                      'Position', [370 0 50 25]);
            
            h = uicontrol(tabresults, styles.accentbutton, 'String', 'Plot factors', ...
                          'Tag', 'btnfactors', ...
                          'Position', [0 10 150 35]);
            
            %% Export
            tabexport = uipanel(col2, styles.panel, 'Title', 'Export to workspace', ...
                                'Tag', 'tabexport', ...
                                'Position', [0 0 0 55]);
            
                        
            p = uipanel(tabexport, styles.place, 'Position', [margin margin width 30]);
            uicontrol(p, styles.label, 'String', 'Factors', ...
                      'Position', [0   5 150 15]);
            uicontrol(p, styles.edit, 'Tag', 'txtexportU', ...
                      'String', '', ...
                      'Position', [155 1 130 24]);
            uicontrol(p, styles.button, 'Tag', 'btnexport', ...
                      'String', 'Export', ...
                      'Position', [290 1 130 25]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Column 3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            margin = 10;
            width  = 200 - 3*margin;

            tabhistory = uipanel(col3, styles.panel, 'Title', 'History', ...
                                 'Tag', 'tabhistory', ...
                                 'Position', [0 0 0 570]);

            
            p = uipanel(tabhistory, styles.place, 'Position', [margin margin width 540]);
            uicontrol(p, styles.list, 'Tag', 'lsthistory', 'Position', [1 60 width-1 480]);
            uicontrol(p, styles.button, 'Tag', 'btnhistrename', ...
                      'String', 'Rename', ...
                      'Position', [0 30 width 25]);
            uicontrol(p, styles.button, 'Tag', 'btnhistdelete', ...
                      'String', 'Delete', ...
                      'Position', [0 0 width 25]);

            
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
            % this.setPosition(findall(this.fh, 'Tag', 'tabresultbtns'), 'hcenter');           
            this.setPosition(findall(this.fh, 'Tag', 'btncompute'), 'hcenter');
            this.setPosition(findall(this.fh, 'Tag', 'btninspect'), 'hcenter');
            this.setPosition(findall(this.fh, 'Tag', 'btnfactors'), 'hcenter');
            % this.setPosition(findall(this.fh, 'Tag', 'btnexport'), 'hcenter');
        end 
        
        function reset(this, ~, ~)
            
        end 
        
    end 
    
end 
