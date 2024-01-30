classdef MLSVDPresenter < handle

    properties
        model
        view 
    end 
    
    properties (SetAccess=private)
        listeners
    end 
    
    methods
        
        function this = MLSVDPresenter(model)
            this.model = model;
            this.view = tensorlab.gui.mlsvd.MLSVDView();
            this.create();
            this.model.reset();
            this.resetView();
            % Ugly hack for matlabs figure scaling behavior
            for k = 1:5
                drawnow();
                pause(0.1);
                pos = get(this.view.fh, 'Position');
                if all(round(pos(3:4) == [855 432])), break; end 
                pos(3:4) = [855 432];
                set(this.view.fh, 'Position', pos);
            end 
            movegui(this.view.fh, 'center');
            set(this.view.fh, 'Visible', 'on');
        end 
        
        function create(this, h, evt, handles)
            import tensorlab.gui.common.*;
            view = this.view;
                        
            %% Add callbacks to menus and buttons
            % Data tab
            view.setCallback('btnfromfile', @this.cbLoadFromFile);
            view.setCallback('btnfromworkspace', @this.cbLoadFromWorkspace);
            view.setCallback('btnreset', @this.cbReset);
            % Button groups 
            view.setProperty('bgalgorithm', 'SelectionChangedFcn', @this.cbChangeAlgorithm);
            view.setProperty('bgrefinement', 'SelectionChangedFcn', @this.cbChangeRefinement);
            view.setCallback('txtsizecore', @this.cbSizeCore);
            % Compute button
            view.setCallback('btncompute', @this.cbCompute);
            view.setCallback('btnhelprank', @this.cbHelpRank);
            view.setCallback('btnvisualize', @this.cbVisualize);
            view.setCallback('btnerrslice', @this.cbPlotAbsErrSlice);
            view.setCallback('btnrelerrslice', @this.cbPlotRelErrSlice);
            view.setCallback('btninspect', @this.cbInspect);
            view.setCallback('btnexport', @this.cbExport);
            % menus
            view.setCallback('mnuloadfile', @this.cbLoadFromFile);
            view.setCallback('mnuloadworkspace', @this.cbLoadFromWorkspace);
            view.setCallback('mnureset', @this.cbReset);
            view.setCallback('mnugeneratecode', @this.cbGenerateCode);
            view.setCallback('mnuexit', @this.cbExit);
            view.setCallback('mnucite', @this.cbCite);
            view.setCallback('mnuwebsite', @this.cbWebsite);
            view.setCallback('mnuabout', @this.cbAbout);
            view.setCallback('mnuuserguide', @this.cbUserguide);
            % Set data
            view.setProperty('rbmlsvd', 'UserData', 'mlsvd');
            view.setProperty('rbmlsvdrsi', 'UserData', 'mlsvd_rsi');
            view.setProperty('rblmlraaca', 'UserData', 'lmlra_aca');
            view.setProperty('rbnone', 'UserData', 'none');
            view.setProperty('rblmlraminf', 'UserData', 'minf');
            view.setProperty('rblmlranls', 'UserData', 'nls');
            % Export function
            view.setCallback('txtexportU', @this.cbChangeExport);
            view.setProperty('txtexportU', 'UserData', 'exportU');
            view.setCallback('txtexportS', @this.cbChangeExport);
            view.setProperty('txtexportS', 'UserData', 'exportS');
            view.setCallback('txtexportT', @this.cbChangeExport);
            view.setProperty('txtexportT', 'UserData', 'exportT');
            
            %% Add input validation
            try 
                jtxtsizecore = findjobj(view.fh, 'txtsizecore');
                set(jtxtsizecore, 'KeyPressedCallback', @this.cbSizeValidation);
            end 
            
            %% Remove old listeners, if any
            delete(this.listeners);
            
            %% Add listener to all relevant fields in the model
            fields = {'T', 'path', 'varname', 'N', 'type', 'size_tens', 'size_core', ...
                      'algorithm', 'refinement', 'ratio', 'abserr', 'relerr', 'maxabserr', ...
                      'maxabserrind', 'maxrelerr', 'maxrelerrind', 'exportU', 'exportS', ...
                      'exportT'}; 
            l = event.proplistener.empty();
            for f = fields
                l(end+1) = addlistener(this.model, f, 'PostSet', @this.handleModelEvents);
            end 
            l(end+1) = addlistener(this.model, 'Reset', @this.resetView);
            this.listeners = l;
        end 
        
        function varargout = output(this, h, ~, ~)
            varargout = {h};
        end 
       
        function resetView(this, ~, ~)
            this.view.setText('txtsizecore', '');
            this.view.setProperty('btnreset', 'Enable', 'off');
            this.view.setProperty('mnureset', 'Enable', 'off');
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Listeners: communication from model to view
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function handleModelEvents(this, src, evt)
            import tensorlab.auxiliary.*;
            % Check for deleted handle
            view  = this.view;
            if ~isvalid(view.fh), return; end
            model = evt.AffectedObject;

            switch lower(src.Name)
              case 'algorithm'
                view.setPanelState('tabresults', false);
                view.setPanelState('tabexport', false);
                switch model.algorithm 
                  case {'mlsvd', 'mlsvds'}, view.setProperty('rbmlsvd', 'Value', 1);
                  case 'mlsvd_rsi', view.setProperty('rbmlsvdrsi', 'Value', 1);
                  case 'lmlra_aca', view.setProperty('rblmlraaca', 'Value', 1);
                end 
              case 'refinement'
                view.setPanelState('tabresults', false);
                view.setPanelState('tabexport', false);
                switch model.refinement 
                  case 'none', view.setProperty('rbnone', 'Value', 1);
                  case 'nls', view.setProperty('rblmlranls', 'Value', 1);
                  case 'minf', view.setProperty('rblmlraminf', 'Value', 1);
                end 
              case 't'
                % check logic
                view.setPanelState('tabsettings', ~isempty(model.T));
                if ~isempty(model.T)
                    if model.hasValidSizeCore(), view.setProperty('btncompute', 'Enable', 'on');
                    else view.setProperty('btncompute', 'Enable', 'off'); end 
                end 
                view.setProperty('mnugeneratecode', 'Enable', 'off');
              case 'path'
                view.setText('txtpath', model.path);
              case 'varname'
                view.setText('txtvariable', model.varname);
                view.setProperty('btnreset', 'Enable', 'on');
                view.setProperty('mnureset', 'Enable', 'on');
              case 'type'
                view.setText('txttype', model.type);
                if strcmpi(model.type, 'sparse')
                    view.setText('rbmlsvd', 'MLSVDS');
                    view.setProperty('rbmlsvd', 'UserData', 'mlsvds');
                else 
                    view.setText('rbmlsvd', 'MLSVD');
                    view.setProperty('rbmlsvd', 'UserData', 'mlsvd');
                end                 
                if view.getProperty('rbmlsvd', 'Value') == 1
                    view.setProperty('rbmlsvd', 'Value', 1); % trigger change
                end 
              case 'n'
                view.setText('txtorder', int2str(model.N));
              case 'size_tens'
                if ~isempty(model.size_tens)
                    view.setText('txtsize', size2str(model.size_tens, 'x', false));
                else
                    view.setText('txtsize', '');
                end 
              case 'size_core'
                if numel(model.size_computed) ~= numel(model.size_core) || ...
                        any(model.size_computed  ~= model.size_core)
                    view.setPanelState('tabresults', false);
                    view.setPanelState('tabexport', false);
                end 
                if ~isempty(model.size_core)
                    view.setText('txtsizecore', size2str(model.size_core, ' ', false));
                    if model.hasValidSizeCore(), state = 'on';
                    else state = 'off'; end
                    view.setProperty('btncompute', 'Enable', state);
                    view.setProperty('mnugeneratecode','Enable', state);
                else 
                    view.setText('txtsizecore', '');
                    view.setProperty('btncompute', 'Enable', 'off');
                    view.setProperty('mnugeneratecode', 'Enable', 'off');
                end
              case 'ratio'
                view.setText('txtratio', sprintf('%.2f', model.ratio));
              case 'abserr'
                view.setText('txtabserr', sprintf('%6.3e', model.abserr));
              case 'relerr'
                view.setText('txtrelerr', sprintf('%6.3e', model.relerr));
                if ~isempty(model.relerr)
                    view.setText('txtrelfit', sprintf('%.2f %%', (1-model.relerr)*100));
                else 
                    view.setText('txtrelfit', '');
                end 
              case 'maxabserr'
                view.setText('txtmaxabserr', sprintf('%6.3e', model.maxabserr));
              case 'maxabserrind'
                if ~isempty(model.maxabserrind)
                    view.setText('txtmaxabserrind', size2str(model.maxabserrind, ' '));
                else 
                    view.setText('txtmaxabserrind', '');
                end 
              case 'maxrelerr'
                view.setText('txtmaxrelerr', sprintf('%6.3e', model.maxrelerr));
              case 'maxrelerrind'
                if ~isempty(model.maxrelerrind)
                    view.setText('txtmaxrelerrind', size2str(model.maxrelerrind, ' '));
                else 
                    view.setText('txtmaxrelerrind', '');
                end 
              case 'exportu'
                view.setText('txtexportU', model.exportU);
              case 'exports'
                view.setText('txtexportS', model.exportS);
              case 'exportt'
                view.setText('txtexportT', model.exportT);
            end 
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Callbacks: communication from GUI to model.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cbCompute(this, h, evt)
        % hObject    handle to calculate (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        % set(findall(handles.figure1,'-property','enable'),'enable','off')
            
        % Make GUI inactive during computation
            if ~this.model.hasValidSizeCore()
                errordlg('The chosen multilinear rank is not valid.');
                return;
            end 
            this.view.setPanelState('tabsettings', 'off');
            this.view.setPanelState('tabdata', 'off');
            this.view.setPanelState('tabresults', 'off');
            this.view.setPanelState('tabexport', 'off');
            this.view.setProperty('mnumain', 'enable', 'off');
            this.view.setProperty('mnuhelp', 'enable', 'off');
            set(this.view.fh, 'Pointer', 'watch');
            drawnow();
            try
                this.model.compute();
            catch e
                set(this.view.fh,'Pointer','arrow')
                drawnow();
                rethrow(e)
            end 
            % Enable all GUI elements
            set(findall(this.view.fh,'-property','enable'),'enable','on')
            set(this.view.fh,'Pointer','arrow')
            drawnow();
        end 
        
        function cbExport(this, ~, ~)
            this.model.cbExportToWorkspace();
        end 
        
        function cbSizeCore(this, h, evt)
            import tensorlab.auxiliary.*;
            
            size_core = this.view.getText('txtsizecore');
            if isempty(size_core), this.model.size_core = []; return; end
            this.view.setProperty('txtsizecore', 'ForegroundColor', 'red');
            if isempty(regexp(size_core, '^[0-9,; ]+$'))
                errordlg(sprintf(['The multilinear rank should consist of %d nonnegative numbers, ' ...
                                  'separated by commas, semicolons or spaces.'], this.model.N));
                return;
            end 
            size_core = regexp(size_core, '[ ,;]+', 'split');
            size_core(cellfun(@isempty, size_core)) = [];
            size_core = cellfun(@str2num, size_core);
            if numel(size_core) ~= this.model.N || any(~isnonnegint(size_core,true)) 
                errordlg(sprintf(['The multilinear rank should consist of %d nonnegative numbers, ' ...
                                  'separated by commas, semicolons or spaces.'], this.model.N));
                return;
            end 
            if any(size_core > this.model.size_tens)
                errordlg('Each mode-n rank is limited by the corresponding tensor dimension.');
                return;
            end 
            [m,i] = max(size_core); 
            if m > prod(size_core([1:i-1 i+1:end]))
                errordlg(['The maximal mode-n rank in a given mode is (theoretically) limited by the ' ...
                          'product of the mode-n ranks of the other modes.']);
                return;
            end
            this.view.setProperty('txtsizecore', 'ForegroundColor', 'black');
            this.model.size_core = size_core;
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
                this.view.setProperty('txtsizecore', 'ForegroundColor', 'black');
            else 
                this.view.setProperty('txtsizecore', 'ForegroundColor', 'red');
            end 
        end
        
        function cbChangeExport(this, h, evt)
            val = get(h, 'String');
            isvalid = true;
            if isempty(val)
                set(h, 'String', this.model.(get(h,'UserData')));
                errordlg('The export name cannot be empty.');
                isvalid = false;
            elseif isempty(regexp(val, '^[a-zA-Z][a-zA-Z0-9_]*$'))
                errordlg('The given name is not a valid name for a variable.');
                isvalid = false;
            end 
            if isvalid 
                switch lower(get(h, 'UserData'))
                  case 'exportU', this.model.exportU = val;
                  case 'exportS', this.model.exportS = val;
                  case 'exportT', this.model.exportT = val;
                end
            else
                % reset fields
                this.model.exportU = this.model.exportU;
                this.model.exportS = this.model.exportS;
                this.model.exportT = this.model.exportT;
            end 
        end 
        
        function cbLoadFromFile(this, ~, ~, ~)
            [file, path] = uigetfile({'*.mat';},'File Selector');
            if path ~= 0
                path = fullfile(path,file);
                name = tensorlab.gui.common.selectTensor(path, true);
                if ~isempty(name)
                    this.model.path = path;
                    this.model.varname = name;
                    this.model.loadFromFile();
                end
            end
        end 

        function cbLoadFromWorkspace(this, ~, ~, ~)
            import tensorlab.gui.common.*;
        % retrieve value from workspace
            name = selectTensor();
            if isempty(name), return; end 
            T = evalin('base', name);
            % If something actually selected, set model new tensor
            if isempty(T), return; end
            T = fmt(T);
            if isvalidtensor(T)
                if strcmp('incomplete', getstructure(T))
                    errordlg('Only full tensors are supported!')
                else
                    this.model.path    = 'Workspace';
                    this.model.T       = T;
                    this.model.varname = name;
                end
            else
                errordlg(sprintf('The selected tensor %s is invalid. See isvalidtensor(%s).', name, ...
                                 name));
            end
        end 
        
        function cbReset(this, ~, ~, ~)
            answer = questdlg('Are you sure you want to reset the GUI?','Attention','Ok', ...
                              'Cancel','Ok');
            if strcmpi(answer,'Ok')
                this.model.reset()
            end
        end
        
        function cbInspect(this, ~, ~, ~)
            import tensorlab.gui.mlsvd.*;
            model = ExplorationModel();
            model.T = this.model.T;
            model.varname = this.model.varname;
            presenter = tensorlab.gui.mlsvd.ExplorationPresenter(model, @this.cbInspectResult);
        end 
        
        function cbInspectResult(this, size_core)
            if nargin < 2, size_core = []; end
            if ~isempty(size_core), this.model.size_core = size_core; end
        end 
        
        function cbHelpRank(this, ~, ~, ~)
            web('https://www.tensorlab.net/doc/lmlra.html#choosing-the-size-of-the-core-tensor')
        end 
        
        function cbChangeAlgorithm(this, ~, evt, ~)
            this.model.algorithm = get(evt.NewValue, 'UserData');
        end 
        
        function cbChangeRefinement(this, ~, evt, ~)
            this.model.refinement = get(evt.NewValue, 'UserData');
        end 

        function cbVisualize(this, ~, ~, ~)
            this.model.plotres();
        end 
        
        function cbPlotAbsErrSlice(this, ~, ~, ~)
            this.model.plot_errslice('absolute');
        end 

        function cbPlotRelErrSlice(this, ~, ~, ~)
            this.model.plot_errslice('relative');
        end 

        %% Menu callbacks 
        function cbGenerateCode(this, ~, ~, ~)
            sug = this.model.getExportSuggestion();
            [file,path] = uiputfile('*.m','Export Code to file', [sug '.m']);
            if path ~= 0, this.model.saveCode(fullfile(path,file)); end
        end 
        
        function cbExit(this, ~, ~, ~)
            close(this.view);
        end 
        
        function cbCite(this, ~, ~, ~)
            tensorlab.gui.common.helpel(3);        
        end 
        
        function cbWebsite(this, ~, ~, ~)
            web('https://www.tensorlab.net');
        end 
        
        function cbAbout(this, ~, ~, ~)
            tensorlab.gui.common.helpel(2);
        end 
        
        function cbUserguide(this, ~, ~, ~)
            tensorlab.gui.common.helpel(1);
        end 
        
    end 
    
end

