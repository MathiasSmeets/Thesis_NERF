classdef GuiView < handle
    
    properties
        fh
        styles 
        handles
    end 
    
    methods
        
        function this = GuiView(varargin)
            p = inputParser;
            p.addOptional('HasMinimumSize', false);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            options = p.Results;

            this.fh = figure('Visible', 'off');
            drawnow();
            pause(0.1);
            this.styles = tensorlab.gui.common.getDefaultStyles();
            this.createLayout();
            % movegui(this.fh, 'onscreen');
            % movegui(this.fh, 'center');
            this.resize();
            set(this.fh, 'Visible', 'on');
            if strcmpi(get(this.fh, 'Resize'), 'on')
                if options.HasMinimumSize
                    tensorlab.gui.common.setminsize(this.fh);
                end 
                set(this.fh, 'SizeChangedFcn', @this.resize);
            end 
        end
        
        function row = addRow(this, parent, name, txt, height, allowedit)
            if numel(height) < 2, height = ones(1,2)*height; end
            row(1,1) = uicontrol(parent, this.styles.label, ...
                                 'Tag', ['lbl' name], ...
                                 'String', txt, ...
                                 'Position', [0 0 0 height(1)]);
            if nargin >= 6 && allowedit
                row(1,2) = uicontrol(parent, this.styles.edit, ...
                                     'Tag', ['txt' name], ...
                                     'Position', [0 0 0 height(2)]);
            else 
                row(1,2) = uicontrol(parent, this.styles.value, ...
                                     'Tag', ['txt' name], ...
                                     'Position', [0 0 0 height(2)]);
            end 
        end 
        
        function setPosition(this, h, pos, h2, margin)
            parent = get(h, 'Parent');
            hdims = get(h, 'Position');
            if strcmpi(pos, 'hcenter')
                pdims = get(parent, 'Position');
                hdims(1) = 0.5*(pdims(3)-hdims(3));
            elseif strcmpi(pos, 'above')
                if nargin < 5, margin = 5; end
                h2dims = get(h2, 'Position');
                hdims(2) = h2dims(2) + h2dims(4) + margin;
            end 
            set(h, 'Position', hdims);
        end
                
        function setText(this, name, str)
            if isempty(str), str = ''; end 
            h = this.findall(this.fh, 'tag', name);
            if ~isempty(h), set(h, 'String', str); end
        end 
        
        function text = getText(this, name, str)
            h = this.findall(this.fh, 'tag', name);
            if ~isempty(h), text = get(h, 'String'); 
            else text = ''; end
        end

        function setProperty(this, name, prop, val)
            h = this.findall(this.fh, 'tag', name);
            if ~isempty(h), set(h, prop, val); 
            else error('Object %s not found!', name); end
        end 
        
        function val = getProperty(this, name, prop)
            h = this.findall(this.fh, 'tag', name);
            if ~isempty(h), val = get(h, prop); 
            else error('Object %s not found!', name); end
        end 
        
        function setCallback(this, name, cb)
            h = this.findall(this.fh, 'tag', name);
            if ~isempty(h), set(h, 'Callback', cb); end            
        end

        function setEnabled(this, h, isenabled)
            if ~ischar(isenabled)
                if isenabled, state = 'on';
                else state = 'off'; end
            else 
                state = isenabled;
            end 
            h = this.findall(this.fh, 'Tag', h);
            if ~isempty(h), set(h, 'Enable', state); end 
        end 

        
        function setPanelState(this, panel, isenabled)
            if ~ischar(isenabled)
                if isenabled, state = 'on';
                else state = 'off'; end
            else 
                state = isenabled;
            end 
            h = this.findall(this.fh, 'Tag', panel);
            if ~isempty(h)
                set(findall(h, '-property', 'Enable'), 'Enable', state);
            end 
        end 
        
        function drawLine(this, h, pos, width)
            tensorlab.gui.common.drawline(h, pos, width);
        end 
        
        function match = findall(this, h, varargin)
            if isempty(this.handles)
                handles = findall(ancestor(h, 'Figure'), '-property', 'Tag');
                tags    = get(handles, 'Tag');
                handles = handles(~cellfun(@isempty, tags));
                tags    = tags(~cellfun(@isempty, tags));
                % remove nonunique tags
                [~,~,i] = unique(tags);
                ind     = find(histc(i,1:max(i)) > 1);
                for k = 1:numel(handles)
                    if any(ind-i(k)==0), continue; end
                    this.handles.(tags{k}) = handles(k);
                end 
            end 
            if nargin >= 4 && strcmpi(varargin{1},'Tag') && isfield(this.handles, varargin{2})
                match = this.handles.(varargin{2});
            else 
                match = findall(h, varargin{:});
            end 
        end 
        
    end 
    
    methods 
        
        function createLayout(this)
        % do something
        end 
        
        function resize(this, ~, ~)
            
        end 
        
        function reset(this, ~, ~)
            
        end 
        
    end 
            
end