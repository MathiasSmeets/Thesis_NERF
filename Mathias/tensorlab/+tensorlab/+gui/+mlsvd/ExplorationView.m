classdef ExplorationView < tensorlab.gui.common.GuiView
    
    properties
        %fh 
        tableft
        tabright
        tabaxes
        width = 350;
        
        %styles;
        
        % Right tab element
        ax 
        maxAxes = 4;
        useSlider = false;
        
        cbEndDrag;
        linepos;
    end 
    
    methods
        
        function this = ExplorationView()
            this = this@tensorlab.gui.common.GuiView('HasMinimumSize', true);
        end

        function createDefaultStyles(this)
            this.styles = tensorlab.gui.common.getDefaultStyles();
        end 
        
        function createLayout(this)
            import tensorlab.gui.common.layoutGrid;
            
            fh = this.fh;
            styles = this.styles;
            % Create empty figure
            set(fh, 'MenuBar', 'none');
            set(fh, 'ToolBar', 'none');
            set(fh, 'Visible', 'off');
            set(fh, 'WindowStyle', 'normal');
            set(fh, 'Units', 'pixels');
            set(fh, 'NumberTitle', 'off', 'Name', 'Inspect multilinear singular value estimates');
                        
            % Create left and right panel. The sizes are set by resize.
            this.tableft  = uipanel(fh, 'BorderType', 'none');
            this.tabright = uipanel(fh, 'BorderType', 'none');
            this.tabaxes  = uipanel(this.tabright, 'BorderType', 'none', 'Tag', 'tabaxes');            
            h = uipanel(this.tabright, 'BorderType', 'none', 'Tag', 'tabslider', ...
                        'Visible', 'off');
            uicontrol(h, styles.slider, 'Tag', 'sldaxes', 'Visible', 'off');
            
            %% Create left panels
            tabproperties = uipanel(this.tableft, styles.panel, ...
                                    'Title', 'Tensor properties', ...
                                    'Position', [0 0 this.width-10 110]);
            props(1,1) = uicontrol(tabproperties, styles.label, 'String', 'Tensor');
            props(1,2) = uicontrol(tabproperties, styles.value, 'Tag', 'txttensor');
            props(2,1) = uicontrol(tabproperties, styles.label, 'String', 'Order');
            props(2,2) = uicontrol(tabproperties, styles.value, 'Tag', 'txtorder');
            props(3,1) = uicontrol(tabproperties, styles.label, 'String', 'Size');
            props(3,2) = uicontrol(tabproperties, styles.value, 'Tag', 'txtsize');
            props(4,1) = uicontrol(tabproperties, styles.label, 'String', 'Type');
            props(4,2) = uicontrol(tabproperties, styles.value, 'Tag', 'txttype');
            layoutGrid(props, [10 5], [150 130]);
            
            %% Compute panel
            tabcompute = uipanel(this.tableft, styles.panel, ...
                                 'Tag', 'tabcompute', ...
                                 'Title', 'Estimate', ...
                                 'Position', [0 0 this.width-10 195]);
            
            uicontrol(tabcompute, styles.comment, 'String', ['Indicate how many singular values ' ...
                                'you want to assess and an algorithm for their estimation.'], ...
                      'Position', [10 130 this.width-30, 45]);
            
            props = [];
            props(1,1) = uicontrol(tabcompute, styles.label, 'String', 'Range of inspection', ...
                                   'ToolTipString', 'One number per mode');
            props(1,2) = uicontrol(tabcompute, styles.edit,  'Tag', 'txtsizeupper', ...
                                   'HorizontalAlignment', 'center');
            props(2,1) = uicontrol(tabcompute, styles.label, 'String', 'Algorithm');
            props(2,2) = uibuttongroup(tabcompute, styles.buttongroup, ...
                                       'Tag', 'bgalgorithm', ...
                                       'Position', [0 0 130 61]);            
            layoutGrid(props, [10 45], [160 150], 10);
            h = props(2,2);
            props = [];
            props(1,1) = uicontrol(h, styles.radio, 'Tag', 'rbmlsvdrsi', 'String', 'mlsvd_rsi');
            props(2,1) = uicontrol(h, styles.radio, 'Tag', 'rblmlraaca', 'String', 'lmlra_aca');
            props(3,1) = uicontrol(h, styles.radio, 'Tag', 'rbmlsvd', 'String', 'mlsvd');
            layoutGrid(props, [0 5], 160);
            uicontrol(tabcompute, styles.button, ...
                      'String', 'Estimate', ...
                      'Tag', 'btncompute', ...
                      'Position', [75 10 this.width-150 30]);
            
            %% Plot options panel
            taboptions = uipanel(this.tableft, styles.panel, ...
                                 'Tag',   'taboptions', ...
                                 'Title', 'Plot options', ...
                                 'Position', [0 0 this.width-10 80]);
            props = [];
            props(1,1) = uicontrol(taboptions, styles.label, 'String', 'YScale');
            props(1,2) = uibuttongroup(taboptions, styles.buttongroup, ...
                                       'Tag', 'bgyscale', ...
                                       'Position', [0 0 130 23]);
            props(2,1) = uicontrol(taboptions, styles.label, 'String', 'Plot');
            props(2,2) = uibuttongroup(taboptions, styles.buttongroup, ...
                                       'Tag', 'bgplotfun', ...
                                       'Position', [0 0 130 23]);
            layoutGrid(props, [10 5], [100 210], 1);

            bg1 = props(1,2);
            bg2 = props(2,2);
            props = [];
            props(1,1) = uicontrol(bg1, styles.radio, 'Tag', 'rblin', 'String', 'Linear');
            props(1,2) = uicontrol(bg1, styles.radio, 'Tag', 'rblog', 'String', 'Log');
            layoutGrid(props, [0 5], [90 130]);
            props = [];
            props(1,1) = uicontrol(bg2, styles.radio, 'Tag', 'rbmlsv', 'String', 'Values', ...
                                   'ToolTipString', 'Plot multilinear singular values sv{n}.');
            props(1,2) = uicontrol(bg2, styles.radio, 'Tag', 'rbenergy', 'String', 'Sum of squares', ...
                                   'ToolTipString', 'Plot sum of squares, computed as cumsum(sv{n}.^2)/sum(sv{n}.^2)');
            layoutGrid(props, [0 5], [90 130]);

            %% Choice panel
            height = 170;
            tabchoose = uipanel(this.tableft, styles.panel, ...
                                'Tag', 'tabchoose', ...
                                'Title', 'Evaluate', ...
                                'Position', [0 0 this.width-10 height]);
            uicontrol(tabchoose, styles.comment, ...
                      'String', ['Select a multilinear rank below or drag the lines on the graphs ' ...
                                'on the right.'], ...                      
                      'Position', [10 height-70 this.width-25, 50]);
            props = [];
            props(1,1) = uicontrol(tabchoose, styles.label, 'String', 'Multilinear rank');
            props(1,2) = uicontrol(tabchoose, styles.edit,  'Tag', 'txtsizecore', ...
                                   'HorizontalAlignment', 'center');
            props(2,1) = uicontrol(tabchoose, styles.label, 'String', 'Absolute error',...
                                   'ToolTipString', 'Computed as froblmlrares(T, factors, core).');
            props(2,2) = uicontrol(tabchoose, styles.value, 'Tag', 'txtabserr', ...
                                   'HorizontalAlignment', 'left');
            props(3,1) = uicontrol(tabchoose, styles.label, 'String', 'Relative error', ...
                                   'ToolTipString', 'Computed as froblmlrares(T, factors, core)/frob(T).');
            props(3,2) = uicontrol(tabchoose, styles.value, 'Tag', 'txtrelerr', ...
                                   'HorizontalAlignment', 'left');
            props(4,1) = uicontrol(tabchoose, styles.label, 'String', 'Relative fit', ...
                                   'ToolTipString', 'Computed as 1 - relative fit.');
            props(4,2) = uicontrol(tabchoose, styles.value, 'Tag', 'txtrelfit', ...
                                   'HorizontalAlignment', 'left');
            props(5,1) = uicontrol(tabchoose, styles.label, 'String', 'Compression ratio', ...
                                   'ToolTipString', ['Ratio of number of variables (numel(core) ' ...
                                   '+ sum(cellfun(@numel,factors))), and number of tensor ' ...
                                   'entries.']);
            props(5,2) = uicontrol(tabchoose, styles.value, 'Tag', 'txtratio', ...
                                   'HorizontalAlignment', 'left');
            layoutGrid(props, [10 5], [160 150], 2);

            
            %% Finalize panel
            tabfinalize   = uipanel(this.tableft, styles.panel, ...
                                    'Tag', 'tabfinalize', ...
                                    'Title', 'Export', ...
                                    'Position', [0 0 this.width-10 60]);
            uicontrol(tabfinalize, styles.button, ...
                      'String', 'Select and close', ...
                      'Tag', 'btnselect', ...
                      'Position', [75 10 this.width-150 30]);                   
        end 
        
        function createAxes(this, N)
            h = findobj(this.tabright, 'Tag', 'tabaxes');
            Naxes = min(N, this.maxAxes);
            this.ax = cell(1, Naxes);
            for n = 1:Naxes
                this.ax{n} = subplot(Naxes, 1, n, 'parent', h);
                title(sprintf('Mode-%d singular values', n));
            end 
            this.useSlider = Naxes < N;
            if this.useSlider
                h = findobj(this.fh, 'Tag', 'sldaxes');
                set(h, 'Min', 0, 'Max', N-Naxes, 'SliderStep', [1 1]./(N-Naxes));
                set(h, 'Value', N-Naxes);
                this.resize(); % Call to display axes
            end 
            this.linepos = nan(1,N);
        end         
        
        function resize(this, h, evt)
            width = this.width;
            margin = 4;
            origpos = get(this.fh, 'Position');
            pos = origpos;
            pos(3) = max(pos(3), 2*width); 
            % Position left tab tabs
            tabs = get(this.tableft, 'Children');
            heights = cell2mat(get(tabs, 'Position'));
            heights = heights(:,4)+margin;
            heights([1 end]) = heights([1 end]) + 1*margin;
            for k = 1:numel(tabs)
                wh  = get(tabs(k), 'Position');
                tmp = [5, sum(heights(1:k-1))+2*margin, wh(3:4)];
                set(tabs(k), 'Units', 'pixels', 'Position', tmp);
            end
            pos(4) = sum(heights)+margin;
            if nargin < 2 || isempty(h)
                tmp = pos;
                tmp(4) = tmp(4) - 7*margin;
                set(this.fh, 'Position', tmp);
            end     
            set(this.tableft,  ...
                'Units', 'pixels', ...
                'Position', [5, max(origpos(4)-pos(4),0), width, pos(4)]);
            tmpright = [width+5, 0, pos(3)-width+margin, max(origpos(4),pos(4))];
            set(this.tabright, ...
                'Units', 'pixels', ...
                'Position', tmpright);
            tmp = [0, margin, tmpright(3), tmpright(4)];
            if this.useSlider, tmp(3) = tmp(3) - 25; end 
            this.setProperty('tabaxes', 'Units', 'pixels');
            this.setProperty('tabaxes', 'Position', tmp);
                        
            if this.useSlider
                tmp = [tmpright(3)-25-2*margin, margin, 25, tmpright(4)];
                this.setProperty('tabslider', 'Units', 'pixels');
                this.setProperty('tabslider', 'Position', tmp);
                this.setProperty('tabslider', 'Visible', 'on');
                tmp = [margin, margin, 20, tmpright(4)-2*margin];
                this.setProperty('sldaxes', 'Position', tmp);
                this.setProperty('sldaxes', 'Visible', 'on');
            end 
        end 
        
        function setText(this, name, str)
            h = findall(this.fh, 'tag', name);
            if ~isempty(h), set(h, 'String', str); end
        end 
        
        function text = getText(this, name, str)
            h = findall(this.fh, 'tag', name);
            if ~isempty(h), text = get(h, 'String'); 
            else text = ''; end
        end

        function setProperty(this, name, prop, val)
            h = findall(this.fh, 'tag', name);
            if ~isempty(h), set(h, prop, val); 
            else error('Object %s not found!', name); end
        end 
        
        function setCallback(this, name, cb)
            h = findall(this.fh, 'tag', name);
            if ~isempty(h), set(h, 'Callback', cb); end            
        end

        function setPanelState(this, panel, isenabled)
            if ~ischar(isenabled)
                if isenabled, state = 'on';
                else state = 'off'; end
            else 
                state = isenabled;
            end 
            h = findall(this.fh, 'tag', panel);
            if ~isempty(h)
                set(findall(h, '-property', 'Enable'), 'Enable', state);
            end 
        end 
        
        function updateGraphs(this, data, size_core, varargin)
            p = inputParser();
            p.addOptional('YScale', 'linear');
            p.addOptional('Title', 'Mode-%d singular values');
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            options = p.Results;
            
            offset = this.getoffset();
            for n = offset+(1:numel(this.ax))
                h = this.getaxis(n);
                if ~isempty(data)
                    set(h, 'NextPlot', 'replace');
                    stem(h, 1:numel(data{n}), data{n}, 'filled', ...
                         'ButtonDownFcn', @this.cbSelect, ...
                         'UserData', n);
                    set(h, 'YScale', options.YScale);
                    set(h, 'NextPlot', 'add');
                    line(h, 'XData', size_core([n n]), 'YData', ylim(h), ...
                         'linewidth', 1.5, 'Color', [1 0 0], 'Tag', 'select', ...
                         'ButtonDownFcn', @this.cbStartDrag);
                    xlim(h, [1, numel(data{n})]);
                end 
                title(h, sprintf(options.Title, n));
            end 
            this.linepos = size_core;
        end
        
        function offset = getoffset(this)
            if this.useSlider
                h = findall(this.fh, 'Tag', 'sldaxes');
                offset = get(h, 'Max') - round(get(h, 'Value'));
            else 
                offset = 0;
            end 
        end 
        
        function ax = getaxis(this, mode)
            ax = this.ax{mode-this.getoffset()};
        end 
        
        function cbSelect(this, src, ev)
            xval = max(round(ev.IntersectionPoint(1)), 1);
            offset = this.getoffset();
            for n = offset+(1:numel(this.ax))
                try 
                    x = get(findall(this.getaxis(n), 'Tag', 'select'), 'XData');
                    this.linepos(n) = x(1);
                end 
            end 
            this.linepos(src.UserData) = xval;
            this.cbEndDrag();
        end 
        
        function clearGraphs(this)
            for n = 1:numel(this.ax)
                cla(this.ax{n});
            end 
        end
        
        function cbStartDrag(this, ev, src)
            set(this.fh, 'Pointer', 'fleur');
            set(this.fh, 'windowbuttonmotionfcn', {@this.cbDrag, src})
            set(this.fh, 'windowbuttonupfcn',     @this.cbStopDrag)
        end 
        
        function cbDrag(this, ev, ~, src)
            coords = get(gca, 'CurrentPoint');
            x = median([xlim(gca) round(coords(1))]);
            set(src.Source, 'XData', [x x]);
        end 
        
        function cbStopDrag(this, ev, ~)
            set(this.fh, 'Pointer', 'arrow');
            set(this.fh, 'windowbuttonmotionfcn', '');
            set(this.fh, 'windowbuttonupfcn', '');
            if isa(this.cbEndDrag, 'function_handle')
                offset = this.getoffset();
                for n = offset+(1:numel(this.ax))
                    try 
                        x = get(findall(this.getaxis(n), 'Tag', 'select'), 'XData');
                        this.linepos(n) = x(1);
                    end 
                end 
                this.cbEndDrag();
            end 
        end 
    end 
end 