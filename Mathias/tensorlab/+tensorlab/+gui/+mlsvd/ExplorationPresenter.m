classdef ExplorationPresenter < handle 
     
    properties
        model 
        view
        
        listeners
        cbResult
    end 
    
    methods
        
        function this = ExplorationPresenter(model, callback)
            this.model = model;
            if nargin >= 2, this.cbResult = callback; end
            this.view  = tensorlab.gui.mlsvd.ExplorationView();
            this.create();
            this.reset();
            % Ugly hack for matlabs figure scaling behavior
            for k = 1:5
                drawnow();
                pause(0.1);
                pos = get(this.view.fh, 'Position');
                if all(round(pos(3:4) == [700 649])), break; end 
                pos(3:4) = [700 649];
                set(this.view.fh, 'Position', pos);
            end 
            movegui(this.view.fh, 'center');
            set(this.view.fh, 'Visible', 'on');
        end 
        
        function create(this)
            import tensorlab.gui.common.*;
            view = this.view;
            
            %% Add callbacks
            view.setCallback('btncompute',   @this.cbCompute);
            view.setCallback('btnselect',    @this.cbFinalize);
            view.setCallback('txtsizecore',  @this.cbSizeCore);
            view.setCallback('txtsizeupper', @this.cbSizeUpper);
            view.setCallback('sldaxes',      @this.cbModesChanged);
            view.setProperty('bgalgorithm', 'SelectionChangedFcn', @this.cbAlgorithm);
            view.setProperty('bgyscale',    'SelectionChangedFcn', @this.cbYScale);
            view.setProperty('bgplotfun',   'SelectionChangedFcn', @this.cbPlotFun);
            
            view.cbEndDrag = @this.cbEndDrag;
            
            %% Change finalization if no export
            if ~isa(this.cbResult, 'function_handle')
                view.setProperty('tabfinalize', 'Title', 'Close');
                view.setText('btnselect', 'Close');
            end 
            
            %% Add input validation
            try 
                jtxtsizeupper = findjobj(view.fh, 'txtsizeupper');
                set(jtxtsizeupper, 'KeyPressedCallback', @this.cbSizeValidation);
                jtxtsizecore = findjobj(view.fh, 'txtsizecore');
                set(jtxtsizecore, 'KeyPressedCallback', @this.cbSizeValidation);
            end 
            
            %% Remove old listeners
            delete(this.listeners);
            fields = {'T', 'varname', 'N', 'type', 'size_tens', 'size_upper', 'size_core', ...
                      'algorithm', 'ratio', 'abserr', 'relerr', 'yscale', 'plotfun'}; 
            l = event.proplistener.empty();
            for f = fields
                l(end+1) = addlistener(this.model, f, 'PostSet', @this.handleModelEvents);
            end 
            this.listeners = l;
            set(view.fh, 'CloseRequestFcn', @this.cbClose);
            
            %% Trigger event for algorithm
            this.model.algorithm = this.model.algorithm;

            view.createAxes(this.model.N);
        end 
                
        function reset(this)
            import tensorlab.auxiliary.*;
            view  = this.view;
            model = this.model;
            
            model.plotfun   = 'mlsv';
            algos = model.getSupportedAlgorithms();
            model.algorithm = algos{1};
            model.yscale    = 'linear';
            
            if strcmpi(this.model.type, 'sparse')
                view.setProperty('rbmlsvd', 'String', 'mlsvds');
            else 
                view.setProperty('rbmlsvd', 'String', 'mlsvd');
            end 
            
            view.setText('txttensor', model.varname);
            view.setText('txtorder', model.N);
            view.setText('txtsize', size2str(model.size_tens, 'x', false));
            view.setText('txttype', model.type);
            view.setText('txtsizeupper', '');
            view.setText('txtsizecore', '');
            switch model.algorithm
              case 'mlsvd_rsi', view.setProperty('rbmlsvdrsi', 'Value', 1);
              case 'mlsvd',     view.setProperty('rbmlsvd', 'Value', 1);
              case 'lmlra_aca', view.setProperty('rblmlraaca', 'Value', 1);
            end 
            this.setPanelState('tabchoose', false);
            this.setPanelState('taboptions', false);
            this.setPanelState('tabfinalize', false);

            model.size_upper = min(model.size_tens, 50);
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Callbacks
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cbClose(this, h, evt)
            try 
                delete(this.listeners);
            end 
            delete(this.view.fh);
        end 
        
        function cbCompute(this, ~, ~, ~)
            import tensorlab.auxiliary.*;
            
            switch lower(this.model.algorithm)
              case 'mlsvd',     fun = @mlsvd; 
              case 'mlsvds',    fun = @mlsvds; 
              case 'lmlra_aca', fun = @lmlra_aca; 
              case 'mlsvd_rsi', fun = @mlsvd_rsi; 
            end 
            set(this.view.fh,'pointer','watch');
            drawnow();
            fprintf('Performing %s on %s... ', this.model.algorithm, this.model.varname);
            try 
                if strcmpi(this.model.algorithm, 'lmlra_aca')
                    [U, S, sv] = fun(this.model.T, this.model.size_upper, 'FillCore', true);
                    sv = sv.sv;
                else 
                    [U, S, sv] = fun(this.model.T, this.model.size_upper);
                end 
                fprintf('finished\r\n');
            catch 
                fprintf('error\r\n');
            end 
            set(this.view.fh,'pointer','arrow')
            drawnow();
            for n = 1:numel(sv)
                sv{n} = sv{n}(1:min(this.model.size_upper(n),numel(sv{n})));
            end
            this.model.U  = U;
            this.model.S  = S;
            this.model.sv = sv;

            this.model.size_computed = this.model.size_upper;
            this.model.size_core     = this.model.size_upper;
            this.setPanelState('tabchoose', true);
            this.setPanelState('taboptions', true);
        end 

        function cbComputeErrors(this, ~, ~, ~)
            model = this.model;
            idx = arrayfun(@(n,m) 1:min(n,m), model.size_core, cellfun(@numel,model.sv), ...
                           'UniformOutput', false);
            U = cellfun(@(u,ind) u(:,ind), model.U, idx, 'UniformOutput', false);
            S = model.S(idx{:});
            residual = lmlrares(model.T, U, S); 
            model.abserr = frob(residual);
            model.relerr = model.abserr/frob(model.T);
            model.ratio  = prod(model.size_tens) / (sum(cellfun(@numel, U)) + numel(S));
        end 
        
        function updateGraphs(this)
            options = struct;
            options.YScale = this.model.yscale;
            switch this.model.plotfun
              case 'mlsv'
                options.Title = 'Mode-%d singular values';
                data = this.model.sv;
              case 'energy'
                options.Title = 'Relative sum of squares in mode %d';
                data = cellfun(@(sv) cumsum(sv.^2)/sum(sv.^2), this.model.sv, 'UniformOutput', ...
                               false);
              otherwise
                return;
            end 
            this.view.updateGraphs(data, this.model.size_core, options);     
        end 
        
        function clearGraphs(this)
            this.view.clearGraphs();     
        end         
        
        function cbFinalize(this, ~, ~, ~)
            try 
                if isa(this.cbResult, 'function_handle')
                    this.cbResult(this.model.size_core);
                end 
            end 
            close(this.view.fh);
        end 
        
        function cbSizeCore(this, ~, ~, ~)
            import tensorlab.auxiliary.*;
            
            size_core = this.view.getText('txtsizecore');
            if isempty(size_core), this.model.size_core = []; return; end
            this.view.setProperty('txtsizecore', 'ForegroundColor', 'red');
            this.setPanelState('tabfinalize', false);
            if isempty(regexp(size_core, '^[0-9,; ]+$'))
                errordlg(sprintf(['The multilinear rank should consist of %d nonnegative numbers, ' ...
                                  'separated by commas, semicolons or spaces.'], this.model.N));
                this.model.size_core = this.model.size_core;
                return;
            end 
            size_core = regexp(size_core, '[ ,;]+', 'split');
            size_core(cellfun(@isempty, size_core)) = [];
            size_core = cellfun(@str2num, size_core);
            if numel(size_core) ~= this.model.N || any(~isnonnegint(size_core,true)) 
                errordlg(sprintf(['The multilinear rank should consist of %d nonnegative numbers, ' ...
                                  'separated by commas, semicolons or spaces.'], this.model.N));
                this.model.size_core = this.model.size_core;
                return;
            end 
            if any(size_core > this.model.size_upper)
                errordlg('Each mode-n rank is limited by the chosen upper limits.');
                this.model.size_core = this.model.size_core;
                return;
            end 
            [m,i] = max(size_core); 
            if m > prod(size_core([1:i-1 i+1:end]))
                errordlg(['The maximal mode-n rank in a given mode is (theoretically) limited by the ' ...
                          'product of the mode-n ranks of the other modes.']);
                this.model.size_core = this.model.size_core;
                return;
            end
            this.view.setProperty('txtsizecore', 'ForegroundColor', 'black');
            this.setPanelState('tabfinalize', true);
            
            this.model.size_core = size_core;
        end 
        
        function cbSizeUpper(this, ~, ~, ~)
            import tensorlab.auxiliary.*;
            
            isvalid = true;
            size_core = this.view.getText('txtsizeupper');
            if isempty(size_core), this.model.size_core = []; 
                isvalid = false; 
            else 
                this.view.setProperty('txtsizeupper', 'ForegroundColor', 'red');
                this.view.setProperty('btncompute', 'Enable', 'off')
                this.setPanelState('tabchoose', false);
                this.setPanelState('taboptions', false);
                this.setPanelState('tabfinalize', false);
                if isempty(regexp(size_core, '^[0-9,; ]+$'))
                    errordlg(sprintf(['The multilinear rank should consist of %d nonnegative numbers, ' ...
                                      'separated by commas, semicolons or spaces.'], this.model.N));
                    isvalid = false; 
                else             
                    size_core = regexp(size_core, '[ ,;]+', 'split');
                    size_core(cellfun(@isempty, size_core)) = [];
                    size_core = cellfun(@str2num, size_core);
                    if numel(size_core) ~= this.model.N || any(~isnonnegint(size_core,true)) 
                        errordlg(sprintf(['The multilinear rank should consist of %d nonnegative numbers, ' ...
                                          'separated by commas, semicolons or spaces.'], this.model.N));
                        isvalid = false;
                    elseif any(size_core > this.model.size_tens)
                        errordlg('Each mode-n rank is limited by the corresponding tensor dimension.');
                        isvalid = false;
                    else 
                        [m,i] = max(size_core); 
                        if m > prod(size_core([1:i-1 i+1:end]))
                            errordlg(['The maximal mode-n rank in a given mode is (theoretically) ' ...
                                      'limited by the product of the mode-n ranks of the other ' ...
                                      'modes.']);
                            isvalid = false;
                        end
                    end 
                end 
            end 
            try 
                if isvalid 
                    this.view.setProperty('txtsizeupper', 'ForegroundColor', 'black');
                    this.view.setProperty('btncompute', 'Enable', 'on');
                    this.model.size_upper = size_core;
                else 
                    this.view.setProperty('txtsizeupper', 'ForegroundColor', 'red');
                    this.view.setProperty('btncompute', 'Enable', 'off');
                    this.model.size_upper = [];
                    this.clearGraphs();
                end 
            end 
        end
        
        function cbSizeValidation(this, h, evt)
            import tensorlab.auxiliary.*;
            
            size_core = char(h.getText());
            if isempty(size_core)
                isvalid = false;
            else 
                if regexp(size_core, '^[0-9 ,;]+$')
                    size_core = regexp(size_core, '[ ,;]+', 'split');
                    size_core(cellfun(@isempty, size_core)) = [];
                    size_core = cellfun(@str2num, size_core);
                    isvalid = this.model.isValidSizeCore(size_core);
                else 
                    isvalid = false;
                end 
            end 
            if isvalid
                this.view.setProperty(h.name, 'ForegroundColor', 'black');
            else 
                this.view.setProperty(h.name, 'ForegroundColor', 'red');
            end 
        end
        
        function cbAlgorithm(this, ~, evt, ~)
            switch lower(get(evt.NewValue, 'Tag'))
              case 'rbmlsvd',    
                if strcmpi(this.model.type, 'sparse')
                    this.model.algorithm = 'mlsvds';
                else 
                    this.model.algorithm = 'mlsvd';
                end 
              case 'rbmlsvdrsi', this.model.algorithm = 'mlsvd_rsi';
              case 'rblmlraaca', this.model.algorithm = 'lmlra_aca';
            end             
        end 

        function cbYScale(this, ~, evt, ~)
            switch lower(get(evt.NewValue, 'Tag'))
              case 'rblin', this.model.yscale = 'linear';
              case 'rblog', this.model.yscale = 'log';
            end             
        end 
        
        function cbPlotFun(this, src, evt, ~)
            switch lower(get(evt.NewValue, 'Tag'))
              case 'rbmlsv',   this.model.plotfun = 'mlsv';
              case 'rbenergy', this.model.plotfun = 'energy';
            end 
        end 

        function cbEndDrag(this)
            if this.model.isValidSizeCore(this.view.linepos)
                this.model.size_core = this.view.linepos;
            else 
                % reset plot
                errordlg('The selected multilinear rank is not valid.');
                this.model.size_core = this.model.size_core;
            end 
        end 
        
        function cbModesChanged(this, src, ~, ~)
            val = get(src, 'Value');
            % prevent recursion
            if val == round(val), this.updateGraphs();
            else set(src, 'Value', round(val)); end 
        end 
        
        function setPanelState(this, panel, state)
            import tensorlab.auxiliary.*;
            view = this.view;
            model = this.model;
            
            switch panel
              case 'tabcompute'
                algos = model.getSupportedAlgorithms();
                view.setPanelState(panel, state);
                view.setEnabled('rbmlsvdrsi', state && any(strcmpi('mlsvd_rsi', algos)));
                view.setEnabled('rbmlsvd', state && any(strcmpi('mlsvd', algos)));
                view.setEnabled('rblmlraaca', state && any(strcmpi('lmlra_aca', algos)));
              otherwise
                view.setPanelState(panel, state);
            end 
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Handle fields
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function handleModelEvents(this, src, evt)
            import tensorlab.auxiliary.*;
            view  = this.view;
            model = this.model;
            
            switch lower(src.Name)
              case 'size_upper'
                if isempty(model.size_upper) || ~model.isValidSizeCore(model.size_upper)
                    view.setProperty('btncompute', 'Enable', 'off')
                    this.setPanelState('tabchoose', false);
                    this.setPanelState('taboptions', false);
                    this.setPanelState('tabfinalize', false);
                elseif isempty(model.size_computed) || any(model.size_computed ~= model.size_upper)
                    view.setText('txtsizeupper', size2str(model.size_upper, ' ', false));
                    view.setProperty('btncompute', 'Enable', 'on')
                    this.setPanelState('tabchoose', false);
                    this.setPanelState('taboptions', false);
                    this.setPanelState('tabfinalize', false);
                end
              case 'size_core'
                if ~isempty(model.size_core) 
                    view.setText('txtsizecore', size2str(model.size_core, ' ', false));
                    this.view.setProperty('txtsizecore', 'ForegroundColor', 'black');
                    this.setPanelState('tabfinalize', model.hasValidSizeCore());
                    this.cbComputeErrors();
                    this.updateGraphs();
                end 
                if ~model.hasValidSizeCore()
                    this.setPanelState('taboptions', false);
                    this.setPanelState('tabfinalize', false);
                    this.clearGraphs();
                end
              case 'ratio'
                view.setText('txtratio', sprintf('%.2f', model.ratio));
              case 'abserr'
                view.setText('txtabserr', sprintf('%6.3e', model.abserr));
              case 'relerr'
                view.setText('txtrelerr', sprintf('%6.3e', model.relerr));
                view.setText('txtrelfit',  sprintf('%.2f %%', (1-model.relerr)*100));
              case 'algorithm'
                switch model.algorithm
                  case {'mlsvd', 'mlsvds'}, view.setProperty('rbmlsvd', 'Value', 1);
                  case 'mlsvd_rsi', view.setProperty('rbmlsvdrsi', 'Value', 1);
                  case 'lmlra_aca', view.setProperty('rblmlraaca', 'Value', 1);
                end   
                this.setPanelState('tabcompute', true);
              case 'plotfun'
                switch model.plotfun
                  case 'mlsv'
                    view.setProperty('rbmlsv', 'Value', 1);
                  case 'energy'
                    view.setProperty('rbenergy', 'Value', 1);
                end
                this.updateGraphs();
              case 'yscale'
                switch model.yscale
                  case 'linear'
                    view.setProperty('rblin', 'Value', 1);
                  case 'log'
                    view.setProperty('rblog', 'Value', 1);
                end
                this.updateGraphs();
            end 
        end 
        
    end 
    
end 