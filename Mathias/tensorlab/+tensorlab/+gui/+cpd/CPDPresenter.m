classdef CPDPresenter < handle

    properties
        model
        view 
    end 
    
    properties (SetObservable)
        history
    end 

    
    properties (SetAccess=private)
        listeners
    end 
    
    methods
        
        function this = CPDPresenter(model)
            this.model = model;
            this.view = tensorlab.gui.cpd.CPDView();
            this.createFactors();
            this.create();
            this.model.reset();
            %this.resetView();
            drawnow();
            pause(0.1);
            pos = get(this.view.fh, 'Position');
            set(this.view.fh, 'Position', [pos(1:2) 1105 585]);
            movegui(this.view.fh, 'center');
        end 
        
        function createFactors(this)
            h = findall(this.view.fh, 'Tag', 'tabfactors');
            if isempty(this.model.T)
                delete(get(h, 'Children'))
                return;
            end 
            if ~isempty(this.model.symmetrysettings) && this.model.supportsSymmetry();
                N = numel(this.model.symmetrysettings);
                labels = this.model.symmetrysettings;
                [tmp,firstidx] = unique(labels);
                counts = histc(this.model.symmetrysettings, tmp);
                ind    = find(counts > 1);
                cidx   = zeros(1, N);
                color  = 1;
                for k = ind(:).'
                    if cidx(firstidx(k)) > 0, continue; end
                    cidx(labels == k) = color;
                    color = color + 1;
                end
            else 
                N = this.model.N;
                labels = 1:this.model.N;
                cidx = zeros(1, N);
            end
            if     N <= 5,  nbrows = 1; nbcols = N; 
            elseif N <= 8,  nbrows = 2; nbcols = 4;
            elseif N <= 10, nbrows = 2; nbcols = 5;
            elseif N <= 12, nbrows = 2; nbcols = 6;
            elseif N <= 15, nbrows = 3; nbcols = 5;
            else nbrows = ceil(N/5); nbcols = 5; end
            styles = struct;
            styles.regular = {'k', 'LineWidth', 2};
            styles.symmetric = {'LineWidth', 2};
            nonnegmodes = false(1,N);
            if this.model.supportsNonnegativity
                nonnegmodes(this.model.nonnegmodes) = true;
            end 
            delete(get(h, 'Children'));
            width = 1./(1.1*nbcols+0.1);
            height = 1./(1.1*nbrows+0.1);
            for n = 1:N
                ax = axes('parent', h); %subplot(nbrows, nbcols, n, 'parent', h);
                tmp = [0.1*width+1.1*width*mod(n-1,nbcols), 1-1.1*height*ceil(n/nbcols), ...
                       width, height];
                set(ax, 'Position', tmp);
                if cidx(n) ==  0, 
                    style = styles.regular;
                else
                    colors = get(ax, 'ColorOrder');
                    style = {styles.symmetric{:}, 'Color', colors(mod(cidx(n),size(colors,1)),:)};
                end 
                if n <= this.model.N
                    plot(ax, [0 0 1 1 0], [0 1.5 1.5 0 0], style{:});
                else 
                    plot(ax, [0 0 1 1 0], [1.2 1.5 1.5 1.2 1.2], style{:});
                end
                hold(ax, 'on');
                axis(ax, [0 1 0 1.5]);
                axis(ax, 'equal');
                ylim(ax, [0 1.5]);
                set(ax, 'Visible', 'off', 'Ytick',[], 'Xtick', []);
                txt = text(ax, 0.5, 1, sprintf('U\\{%d\\}', labels(n)));
                set(txt,'HorizontalAlignment', 'center');
                if nonnegmodes(n)
                    if n <= this.model.N
                        txt = text(ax, 0.5, 0.5, '+');
                    else 
                        txt = text(ax, 0.5, 1.1, '+');
                    end 
                    set(txt,'HorizontalAlignment', 'center');
                end 
            end 
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
            view.setProperty('bginitmethod', 'SelectionChangedFcn', @this.cbChangeInitMethod);
            view.setProperty('bginitalgorithm', 'SelectionChangedFcn', ...
                                           @this.cbChangeInitAlgorithm);
            view.setProperty('bgalgorithm', 'SelectionChangedFcn', @this.cbChangeAlgorithm);
            view.setCallback('txtrank', @this.cbRank);
            view.setCallback('txtsymmetry', @this.cbChangeSymmetry);
            view.setCallback('txtnonneg', @this.cbChangeNonneg);
            view.setCallback('chkplotcvg', @this.cbChangePlotConvergence);
            % Compute button
            view.setCallback('btninspect', @this.cbInspect);
            view.setCallback('btninitselect', @this.cbSelectInitialization);
            view.setCallback('btnhelprank', @this.cbHelpRank);
            view.setCallback('btncompute', @this.cbCompute);
            view.setCallback('btnoptions', @this.cbOptions);
            view.setCallback('btnvisualize', @this.cbVisualize);
            view.setCallback('btnerrslice', @this.cbPlotErrSlice);
            view.setCallback('btnshowcurves', @this.cbPlotConvergenceCurves);
            view.setCallback('btnshownorms', @this.cbPlotNorms);
            view.setCallback('btnshowangles', @this.cbPlotAngles);
            view.setCallback('btncomputecondition', @this.cbComputeCondition);
            view.setCallback('btnfactors', @this.cbPlotFactors);
            view.setCallback('btnexport', @this.cbExport);
            view.setCallback('btnhistrename', @this.cbHistRename);
            view.setCallback('btnhistdelete', @this.cbHistDelete);
            % history
            view.setCallback('lsthistory', @this.cbChangeHistory);
            % menus
            view.setCallback('mnuloadfile', @this.cbLoadFromFile);
            view.setCallback('mnuloadworkspace', @this.cbLoadFromWorkspace);
            view.setCallback('mnureset', @this.cbReset);
            view.setCallback('mnuloadsession', @this.cbLoadSession);
            view.setCallback('mnusavesession', @this.cbSaveSession);
            view.setCallback('mnugeneratecode', @this.cbGenerateCode);
            view.setCallback('mnuexit', @this.cbExit);
            view.setCallback('mnucite', @this.cbCite);
            view.setCallback('mnuwebsite', @this.cbWebsite);
            view.setCallback('mnuabout', @this.cbAbout);
            view.setCallback('mnuuserguide', @this.cbUserguide);
            % Set data
            view.setProperty('rbmanual', 'UserData', 'manual');
            view.setProperty('rbauto', 'UserData', 'auto');
            view.setProperty('rbgevd', 'UserData', 'cpd_gevd');
            view.setProperty('rbrandn', 'UserData', 'randn');
            view.setProperty('rbrand', 'UserData', 'rand');
            view.setProperty('rbals', 'UserData', 'cpd_als');
            view.setProperty('rbnls', 'UserData', 'cpd_nls');
            % Export function
            view.setCallback('txtexportU', @this.cbChangeExport);
            view.setProperty('txtexportU', 'UserData', 'exportU');
            
            this.view.setEnabled('btnhistrename', false);
            this.view.setEnabled('btnhistdelete', false);

            
            %% Add input validation
            try 
                % jtxtsizecore = findjobj(view.fh, 'txtsizecore');
                % set(jtxtsizecore, 'KeyPressedCallback', @this.cbSizeValidation);
            end 
            
            this.addmodellisteners();
            addlistener(this, 'history', 'PostSet', @this.handleHistoryEvents);
        end 
        
        function addmodellisteners(this)
            %% Remove old listeners, if any
            delete(this.listeners);
            
            %% Add listener to all relevant fields in the model
            fields = {'T', 'path', 'varname', 'N', 'type', 'size_tens', 'initmethod', 'initalgorithm', ...
                      'Uinitname', 'Uinitmanual', 'Uinit', 'R', 'Rmanual', 'exportU', ...
                      'relerr', 'output', 'corcondia', 'conditionnumber', 'normrank1', 'ratiosv', ...
                      'algorithm', 'algorithmoptions', 'symmetrysettings', 'nonnegmodes', ...
                      'congruence', 'name'};
            l = event.proplistener.empty();
            l(1) = addlistener(this.model, fields, 'PostSet', @this.handleModelEvents);
            l(2) = addlistener(this.model, 'Reset', @this.resetView);
            this.listeners = l;
        end 
        
        
        function varargout = output(this, h, ~, ~)
            varargout = {h};
        end 
       
        function resetView(this, ~, ~)
            this.setPanelState('tabdata', true);
            this.view.setProperty('mnureset', 'Enable', 'off');
        end 
        
        function updateresults(this)
            fields = {'Ures', 'output', 'relerr', 'corcondia', 'conditionnumber', 'ratiosv', ...
                       'normrank1', 'congruence', 'exportU'};
            for k = 1:numel(fields)
                this.model.(fields{k}) = this.model.(fields{k});
            end
        end 
        
        function updateallfields(this)
            % import tensorlab.auxiliary.*;
            % model = this.model;
            % view  = this.view;
            
            % switch model.algorithm 
            %   case 'cpd_als', view.setProperty('rbals', 'Value', 1);
            %   case 'cpd_nls', view.setProperty('rbnls', 'Value', 1);
            % end 
            % view.setText('txtpath', model.path);
            % view.setText('txttype', model.type);
            % view.setText('txtorder', int2str(model.N));
            % if ~isempty(model.size_tens)
            %     view.setText('txtsize', size2str(model.size_tens, 'x', false));
            % else
            %     view.setText('txtsize', '');
            % end 
            % tbl = struct('manual', 'rbmanual', 'auto', 'rbauto');
            % view.setProperty(tbl.(model.initmethod), 'Value', 1);
            % tbl = struct('cpd_gevd', 'rbgevd', 'randn', 'rbrandn', 'rand', 'rbrand');
            % view.setProperty(tbl.(model.initalgorithm), 'Value', 1);
            % view.setText('txtrank', num2str(model.Rmanual));
            % view.setText('txtinittensor', model.Uinitname);
            % if ~isempty(model.Uinitmanual)
            %     view.setText('txtinitmanualrank', num2str(size(model.Uinitmanual{1}, 2)));
            % end 
            % if ~isempty(model.symmetrysettings)
            %     view.setText('txtsymmetry', size2str(model.symmetrysettings, ' ', false));
            % else
            %     view.setText('txtsymmetry', '');
            % end 
            % if ~isempty(model.nonnegmodes)
            %     view.setText('txtnonneg', size2str(model.nonnegmodes, ' ', false));
            % else
            %     view.setText('txtnonneg', '');
            % end 
            % view.setProperty('chkplotcvg', 'Value', model.plotconvergence);
            % if ~isstruct(model.output)
            %     info = inf;
            % elseif isfield(model.output.Refinement, 'info')
            %     info = model.output.Refinement.info;
            % else 
            %     info = model.output.Algorithm.info;
            % end
            % switch info
            %   case 1, txt = 'Relative function value tolerance reached.';
            %   case 2, txt = 'Relative step length tolerance reached.';
            %   case 3, txt = 'Maximum number of iterations reached.';
            %   case 4, txt = 'Absolute function value tolerance reached.';
            %   case 5, txt = 'Measure tolerance reached.';
            %   otherwise, txt = '';
            % end 
            % view.setText('txtstopcriterion', txt);
            % view.setText('txtrelerr', sprintf('%6.3e', model.relerr));
            % if isempty(model.corcondia), txt = '';
            % else txt = sprintf('%.2g%%', model.corcondia); end 
            % view.setText('txtcorcondia', txt);
            % view.setText('txtcondition', sprintf('%6.g', model.conditionnumber));

            % view.setText('txtratiosv', sprintf('%6.3e', model.ratiosv));
            % if isempty(model.normrank1)
            %     txt = '';
            % else
            %     nrm = num2cell(model.normrank1);
            %     if numel(nrm) > 3
            %         txt = sprintf('%3.1e, %3.1e, ..., %3.1e', nrm{[1 2 end]});
            %     elseif numel(nrm) == 3
            %         txt = sprintf('%3.1e, %3.1e, %3.1e', nrm{:});
            %     elseif numel(nrm) == 2
            %         txt = sprintf('%3.1e, %3.1e', nrm{:});
            %     elseif numel(nrm) == 1
            %         txt = sprintf('%3.1e', nrm{:});
            %     end 
            % end 
            % view.setText('txtnorms', txt);
            % if isempty(model.congruence)
            %     txt = '';
            % else
            %     angles = acos(abs(model.congruence - eye(size(model.congruence))))/pi*180;
            %     minangle = min(angles(:));
            %     txt = sprintf('%.2f degrees (smallest)', minangle);
            % end 
            % view.setText('txtangles', txt);
            % view.setText('txtexportU', model.exportU); 
            
            fields = properties(this.model);
            for k = 1:numel(fields)
                this.model.(fields{k}) = this.model.(fields{k});
            end  
        end 
        
        function [model,index] = findmodel(this)
            for k = 1:numel(this.history)
                if isequal(this.model, this.history(k))
                    model = this.history(k);
                    index = k;
                    return;
                end 
            end 
            model = [];
            index = 0;
        end 
        
        function changemodel(this, newmodel)
            if isequal(this.model, newmodel), return; end
            delete(this.listeners);
            this.model = copy(newmodel);
            this.addmodellisteners();
            this.updateallfields();
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
                this.setPanelState('tabsettings', true);
                this.setPanelState('tabresults', false);
                this.setPanelState('tabexport', false);
                switch model.algorithm 
                  case 'cpd_als', view.setProperty('rbals', 'Value', 1);
                  case 'cpd_nls', view.setProperty('rbnls', 'Value', 1);
                end 
                this.createFactors();
                this.findandselectmodel(); 
              case 'algorithmoptions'
                this.setPanelState('tabresults', false);
                this.setPanelState('tabexport', false);
                this.findandselectmodel(); 
              case 't'
                % check logic
                this.setPanelState('tabinit', ~isempty(model.T));
                this.setPanelState('tabdata', true);
                if ~isempty(model.T)
                    this.setPanelState('tabsettings', model.hasValidInitialization());
                    view.setEnabled('btncompute', model.hasValidInitialization());
                end 
                view.setProperty('mnugeneratecode', 'Enable', 'off');
                view.setProperty('mnusavesession', 'Enable', 'off');
                % update model
                this.createFactors();
              case 'path'
                view.setText('txtpath', model.path);
              case 'varname'
                view.setText('txtvariable', model.varname);
                this.setPanelState('tabdata', true);
                view.setProperty('mnureset', 'Enable', 'on');
              case 'type'
                view.setText('txttype', model.type);
              case 'n'
                view.setText('txtorder', int2str(model.N));
              case 'size_tens'
                if ~isempty(model.size_tens)
                    view.setText('txtsize', size2str(model.size_tens, 'x', false));
                else
                    view.setText('txtsize', '');
                end 
              case 'initmethod'
                tbl = struct('manual', 'rbmanual', 'auto', 'rbauto');
                view.setProperty(tbl.(model.initmethod), 'Value', 1);
                this.setPanelState('tabinit', ~isempty(model.T));
                this.setPanelState('tabsettings', model.hasValidInitialization());
                this.findandselectmodel(); 
              case 'initalgorithm'
                tbl = struct('cpd_gevd', 'rbgevd', 'randn', 'rbrandn', 'rand', 'rbrand');
                view.setProperty(tbl.(model.initalgorithm), 'Value', 1);
                this.setPanelState('tabsettings', model.hasValidInitialization());
                if strcmpi(model.initalgorithm, 'cpd_gevd') && ~model.isValidForGEVD(model.Rmanual)
                    this.view.setProperty('txtrank', 'ForegroundColor', 'red');
                else 
                    this.view.setProperty('txtrank', 'ForegroundColor', 'black');
                end 
                % Revalidate R
                this.model.Rmanual = this.model.Rmanual;
                this.setPanelState('tabinit', ~isempty(model.T));
                this.setPanelState('tabsettings', model.hasValidInitialization());
                this.findandselectmodel();
              case 'rmanual'
                view.setText('txtrank', num2str(model.Rmanual));
                this.setPanelState('tabinit', ~isempty(model.T));
                this.setPanelState('tabsettings', model.hasValidInitialization());
                this.findandselectmodel();
              case 'uinitname'
                view.setText('txtinittensor', model.Uinitname);
              case 'uinitmanual'
                if ~isempty(model.Uinitmanual)
                    view.setText('txtinitmanualrank', num2str(size(model.Uinitmanual{1}, 2)));
                end 
                this.setPanelState('tabsettings', model.hasValidInitialization());
                this.findandselectmodel();
              case 'symmetrysettings'
                if ~isempty(model.symmetrysettings)
                    view.setText('txtsymmetry', size2str(model.symmetrysettings, ' ', false));
                else
                    view.setText('txtsymmetry', '');
                end 
                if ~model.isValidSymmetrySetting(model.symmetrysettings)
                    view.setProperty('txtsymmetry', 'ForegroundColor', 'red');
                else 
                    view.setProperty('txtsymmetry', 'ForegroundColor', 'black');
                end 
                view.setEnabled('btncompute', model.isComputable());
                this.createFactors();
                this.findandselectmodel();
              case 'nonnegmodes'
                if ~isempty(model.nonnegmodes)
                    view.setText('txtnonneg', size2str(model.nonnegmodes, ' ', false));
                else
                    view.setText('txtnonneg', '');
                end 
                if ~model.isValidNonnegativeMode(model.nonnegmodes)
                    view.setProperty('txtnonneg', 'ForegroundColor', 'red');
                else 
                    view.setProperty('txtnonneg', 'ForegroundColor', 'black');
                end 
                view.setEnabled('btncompute', model.isComputable());
                this.createFactors();
                this.findandselectmodel();
              case 'plotconvergence'
                view.setProperty('chkplotcvg', 'Value', model.plotconvergence);
              case 'output'
                if ~isstruct(model.output)
                    info = inf;
                elseif isfield(model.output.Refinement, 'info')
                    info = model.output.Refinement.info;
                else 
                    info = model.output.Algorithm.info;
                end
                switch info
                  case 1, txt = 'Relative function value tolerance reached.';
                  case 2, txt = 'Relative step length tolerance reached.';
                  case 3, txt = 'Maximum number of iterations reached.';
                  case 4, txt = 'Absolute function value tolerance reached.';
                  case 5, txt = 'Measure tolerance reached.';
                  otherwise, txt = '';
                end 
                view.setText('txtstopcriterion', txt);
                this.setPanelState('tabresults', isfinite(info));
                this.setPanelState('tabexport', isfinite(info));
                view.setEnabled('mnugeneratecode', isfinite(info));
              case 'relerr'
                view.setText('txtrelerr', sprintf('%6.3e', model.relerr));
              case 'corcondia'
                if isempty(model.corcondia), txt = '';
                else txt = sprintf('%.2g%%', model.corcondia); end 
                view.setText('txtcorcondia', txt);
              case 'conditionnumber'
                view.setText('txtcondition', sprintf('%6.g', model.conditionnumber));
              case 'ratiosv'
                view.setText('txtratiosv', sprintf('%6.3e', model.ratiosv));
              case 'normrank1'
                if isempty(model.normrank1)
                    txt = '';
                else
                    nrm = num2cell(model.normrank1);
                    if numel(nrm) > 3
                        txt = sprintf('%3.1e, %3.1e, ..., %3.1e', nrm{[1 2 end]});
                    elseif numel(nrm) == 3
                        txt = sprintf('%3.1e, %3.1e, %3.1e', nrm{:});
                    elseif numel(nrm) == 2
                        txt = sprintf('%3.1e, %3.1e', nrm{:});
                    elseif numel(nrm) == 1
                        txt = sprintf('%3.1e', nrm{:});
                    end 
                end 
                view.setText('txtnorms', txt);
              case 'congruence'
                if isempty(model.congruence)
                    txt = '';
                else
                    angles = acos(abs(model.congruence - eye(size(model.congruence))))/pi*180;
                    minangle = min(angles(:));
                    txt = sprintf('%.2f degrees (smallest)', minangle);
                end 
                view.setText('txtangles', txt);
              case 'exportu'
                view.setText('txtexportU', model.exportU); 
              case 'name'
                if ~isempty(model.name)
                    view.setProperty('tabresults', 'Title', sprintf('Results - %s', model.name));
                else 
                    view.setProperty('tabresults', 'Title',  'Results');
                end 
            end 
        end 
        
        function handleHistoryEvents(this, ~, ~)
            import tensorlab.gui.cpd.*;
            lst = findall(this.view.fh, 'Tag', 'lsthistory');
            if isempty(this.history)
                set(lst, 'String', {});
                this.view.setEnabled('btnhistrename', false);
                this.view.setEnabled('btnhistdelete', false);
                this.view.setEnabled('mnusavesession', false);
            else 
                names = {this.history.name};
                [names, ind] = naturalsort(names);
                set(lst, 'String', names, 'UserData', ind);
                [~,index] = this.findmodel();
                if index > 0, 
                    set(lst, 'Value', find(ind==index,1)); 
                else 
                    set(lst, 'Max', 2, 'Value', []);
                end 
                this.view.setEnabled('btnhistrename', index > 0);
                this.view.setEnabled('btnhistdelete', index > 0);
                this.view.setEnabled('mnusavesession', true);
            end 
        end 
        
        function findandselectmodel(this)
            [hist,ind] = this.findmodel();
            lst = findall(this.view.fh, 'Tag', 'lsthistory');
            if ~isempty(hist)
                this.changemodel(hist);
                set(lst, 'Value', ind, 'Max', 1);
            elseif isempty(hist) && ~isempty(get(lst, 'Value'))
                this.model.resetresults();
                set(lst, 'Max', 2, 'Value', []);
            end 
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% GUI Logic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function setPanelState(this, panel, state)
            import tensorlab.auxiliary.*;
            view = this.view;
            model = this.model;
            
            switch panel
              case 'tabdata'
                view.setPanelState('tabdata', state);
                view.setEnabled('btnfromfile', state && isempty(this.model.T));
                view.setEnabled('btnfromworkspace', state && isempty(this.model.T));
                view.setEnabled('btnreset', state && ~isempty(this.model.T));
                view.setEnabled('btninspect', state && ~isempty(this.model.T));
              case 'tabinit'
                view.setPanelState('tabinit', state);
                view.setPanelState('tabinitmanual', state && strcmpi(model.initmethod, 'manual'));
                view.setPanelState('tabinitauto',   state && strcmpi(model.initmethod, 'auto'));
                view.setEnabled('rbgevd', state && strcmpi(model.initmethod, 'auto') && ...
                                          model.isValidForGEVD());
              case 'tabsettings'
                view.setPanelState('tabsettings', state);
                view.setEnabled('lblsymmetry', state && model.supportsSymmetry());
                view.setEnabled('txtsymmetry', state && model.supportsSymmetry());
                % if model.supportsSymmetry() && ~isempty(model.nonnegmodes)
                %     view.setText('txtsymmetry', size2str(model.symmetrysettings, ' ', false));
                % else 
                %     view.setText('txtsymmetry', '')
                % end 
                view.setEnabled('lblnonneg', state && model.supportsNonnegativity());
                view.setEnabled('txtnonneg', state && model.supportsNonnegativity());
                % if model.supportsNonnegativity() && ~isempty(model.nonnegmodes)
                %     view.setText('txtnonneg', size2str(model.nonnegmodes, ' ', false));
                % else 
                %     view.setText('txtnonneg', '')
                % end 
                view.setEnabled('btncompute', state && model.hasValidInitialization() && ...
                                              model.hasValidSettings());         
              case 'tabresults'
                view.setPanelState('tabresults', state);
                tmp = state && ~isempty(model.output) && isempty(model.conditionnumber);
                view.setEnabled('btncomputecondition', tmp);
                if ~model.isConditionAutoComputed()
                    view.setProperty('btncomputecondition', 'Visible', 'on');
                end
              otherwise
                this.view.setPanelState(panel, state);
            end 
        end 
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Callbacks: communication from GUI to model.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cbCompute(this, h, evt)
            if ~this.model.hasValidInitialization()
                errordlg('The chosen initialization is not valid.');
                return;
            end 
            this.setPanelState('tabinit', false);
            this.setPanelState('tabsettings', false);
            this.setPanelState('tabdata', false);
            this.setPanelState('tabresults', false);
            this.setPanelState('tabexport', false);
            this.view.setProperty('mnumain', 'enable', 'off');
            this.view.setProperty('mnuhelp', 'enable', 'off');
            set(this.view.fh, 'Pointer', 'watch');
            drawnow();
            pause(0.1);
            try
                this.model.compute();
            catch e
                set(this.view.fh,'Pointer','arrow')
                drawnow();
            end 
            % Add to history
            this.model.name = this.model.getSessionSuggestion();
            if isempty(this.history)
                this.history = tensorlab.gui.cpd.CPDModel.empty();
            end 
            this.history(end+1) = this.model.copy(); 
            
            % Enable all GUI elements
            this.setPanelState('tabinit', true);
            this.setPanelState('tabsettings', true);
            this.setPanelState('tabdata', true);
            this.setPanelState('tabresults', true);
            this.setPanelState('tabexport', true);
            this.setPanelState('tabhistory', true);
            this.view.setEnabled('mnumain', true);
            this.view.setEnabled('mnuhelp', true);
            this.model.initmethod = this.model.initmethod;
            this.model.initmethod = this.model.initmethod;
            set(this.view.fh,'Pointer','arrow')
            drawnow();
        end 
        
        function cbComputeCondition(this, h, evt)
            if isempty(this.model.Ures)
                errordlg('No decomposition is computed.');
                return;
            end 
            set(this.view.fh, 'Pointer', 'watch');
            drawnow();
            pause(0.1);
            this.model.computeCondition();
            set(this.view.fh,'Pointer','arrow')
            drawnow();                
        end
        
        function cbExport(this, ~, ~)
            this.model.cbExportToWorkspace();
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
        
        function cbChangeInitMethod(this, src, evt, ~)
            this.model.initmethod = get(evt.NewValue, 'UserData');
        end 
        
        function cbSelectInitialization(this, ~, ~, ~) 
            import tensorlab.gui.common.*;
            import tensorlab.auxiliary.*;
            
            % retrieve value from workspace
            options.SelectionFcn = @(T) strcmpi(getstructure(T),'cpd') && isvalidtensor(T) && ...
                tensorlab.auxiliary.sizeequals(getsize(T),this.model.size_tens);
            options.Name = 'initialization';
            options.RankProperty = true;
            
            name = selectTensor([],false,options);
            if isempty(name), return; end 
            U = evalin('base', name);
            % If something actually selected, set initialization
            if isempty(U), return; end
            if strcmpi(getstructure(U), 'cpd') && isvalidtensor(U)
                if sizeequals(this.model.size_tens, getsize(U))
                    this.model.Uinitname   = name;
                    this.model.Uinitmanual = U;
                else
                    errordlg(sprintf(['The dimensions of the tensor (%s) and the selected ' ...
                                      'initialization (%s) do not match.'], ...
                                     size2str(this.model.size_tens), size2str(getsize(U))));
                end
            else
                errordlg(sprintf('The selected variable %s is not a valid CPD.', name));
            end
        end 
        
        function cbRank(this, src, ~, ~)
            import tensorlab.auxiliary.*;

            R = this.view.getText('txtrank');
            isvalid = true;
            if ~isempty(regexp(R, '^\s*$')), return; end
            if isempty(regexp(R, '^[0-9]+$'))
                isvalid = false;
                errordlg('The chosen rank R should be a nonnegative number.');
            else
                R = str2num(R);
                if ~isnonnegint(R, true)
                    isvalid = false;
                    errordlg('The chosen rank R should be a nonnegative number.');
                else
                    this.model.Rmanual = R;
                    if strcmpi(this.model.initalgorithm, 'cpd_gevd') && ...
                        ~this.model.isValidForGEVD(R)
                        isvalid = false;
                        errordlg(['The chosen rank R is not compatible with the GEVD initialization. ' ...
                                  '(At least two tensor dimensions should be equal to or larger than ' ...
                                  'R.)']);
                    end 
                end 
            end 
            if isvalid
                this.model.Rmanual = R;
                this.view.setProperty('txtrank', 'ForegroundColor', 'black');
            else 
                this.view.setProperty('txtrank', 'ForegroundColor', 'red');
            end 
        end
        
        function cbChangeInitAlgorithm(this, src, evt, ~)
            this.model.initalgorithm = get(evt.NewValue, 'UserData');
        end 
            
        function cbReset(this, ~, ~, ~)
            answer = questdlg(['Are you sure you want to reset the GUI? This erases all session ' ...
                               'data and history.'], 'Attention','Ok', 'Cancel','Ok');
            if strcmpi(answer,'Ok')
                this.model.reset();
                this.history = [];
            end
        end
        
        function cbInspect(this, ~, ~, ~)
            import tensorlab.gui.mlsvd.*;
            model = ExplorationModel();
            model.T = this.model.T;
            model.varname = this.model.varname;
            presenter = tensorlab.gui.mlsvd.ExplorationPresenter(model);
            presenter.cbResult = @this.cbInspectResult;
        end 
        
        function cbChangePlotConvergence(this, src, ~, ~)
            this.model.plotconvergence = get(src, 'Value') == 1;
        end 
        
        function cbPlotConvergenceCurves(this, ~, ~, ~)
            figure('name', sprintf('Convergence curves - %s', this.model.name));
            plot_convergence(this.model.output, this.model.algorithmoptions)
        end 
        
        function cbChangeSymmetry(this, src, ~, ~)
            setting = str2num(get(src, 'String'));
            [isvalid, msg] = this.model.isValidSymmetrySetting(setting);
            if ~isvalid, errordlg(msg); end 
            this.model.symmetrysettings = setting;
        end
        
        function cbChangeNonneg(this, src, ~, ~)
            nonnegmodes = str2num(get(src, 'String'));
            [isvalid, msg] = this.model.isValidNonnegativeMode(nonnegmodes);
            if ~isvalid, errordlg(msg); end 
            this.model.nonnegmodes = nonnegmodes;
        end 
        
        function cbInspectResult(this, size_core)
            if nargin < 2, size_core = []; end
            if ~isempty(size_core), this.model.size_core = size_core; end
        end 
        
        function cbOptions(this, ~, ~, ~)
            import tensorlab.gui.common.*;
            switch this.model.algorithm
              case 'cpd_nls', name = 'NLS';
              case 'cpd_als', name = 'ALS';
            end 
            
            values = this.model.algorithmoptions.(this.model.algorithm);
            defaults = this.model.getDefaultAlgorithmOptions(this.model.algorithm);
            
            options = struct;
            options(1).name       = 'tolfun';
            options(1).string     = 'TolFun';
            options(1).tooltip    = 'Relative function tolerance';
            options(1).value      = values.TolFun;
            options(1).default    = defaults.TolFun;
            options(1).validation = @(v) v >= 0;
            options(1).error      = 'Function tolerance should be a nonnegative number.';

            options(2).name       = 'tolx';
            options(2).string     = 'TolX';
            options(2).tooltip    = 'Relative step tolerance';
            options(2).value      = values.TolX;
            options(2).default    = defaults.TolX;
            options(2).validation = @(v) v >= 0;
            options(2).error      = 'Step tolerance should be a nonnegative number.';
            
            options(3).name       = 'maxiter';
            options(3).string     = 'Maximum iterations';
            options(3).tooltip    = 'Maximum number of iterations';
            options(3).value      = values.MaxIter;
            options(3).default    = defaults.MaxIter;
            options(3).validation = @(v) tensorlab.auxiliary.isnonnegint(v);
            options(3).error      = 'The maximum number of iterations should be a nonnegative integer.';

            values = inputoptionsdlg(name, options);
            if ~isempty(values)
                this.model.algorithmoptions.(this.model.algorithm).TolFun  = values{1};
                this.model.algorithmoptions.(this.model.algorithm).TolX    = values{2};
                this.model.algorithmoptions.(this.model.algorithm).MaxIter = values{3};
            end 
        end 
        
        function cbHelpRank(this, ~, ~, ~)
            web('https://www.tensorlab.net/doc/cpd.html#choosing-the-number-of-rank-one-terms');
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
        
        function cbPlotErrSlice(this, ~, ~, ~)
            this.model.ploterrslice('absolute');
        end 
        
        function cbPlotFactors(this, ~, ~, ~)
            plot_rank1terms(this.model.Ures, 'ModeSettings', this.model.modesettings, 'TermSettings', ...
                            this.model.termsettings, 'ExportName', this.model.getExportSuggestion, ...
                            'Title', this.model.name);
        end 

        function cbPlotNorms(this, ~, ~, ~)
            figure('name', sprintf('Norms of rank-1 terms - %s', this.model.name));
            stem(this.model.normrank1, 'filled');
            xlabel('r');
            ylabel('Norm rank-1 term');
        end 
        
        function cbPlotAngles(this, ~, ~, ~)
            figure('name', sprintf('Angles between rank-1 terms - %s', this.model.name));
            plot_angles(this.model.Ures);
        end 
        
        function cbChangeHistory(this, src, ~, ~)
            if isempty(this.history)
                this.view.setEnabled('btnhistrename', false);
                this.view.setEnabled('btnhistdelete', false);
                return
            end
            idx = get(src, 'Value');
            if numel(idx) > 1 
                idx = idx(1);
                set(src, 'Value', idx);
            end 
            set(src, 'Max', 1);
            if idx < 1 || isempty(idx)
                this.view.setEnabled('btnhistrename', false);
                this.view.setEnabled('btnhistdelete', false);
                return
            end 
            data = get(src, 'Userdata');
            val  = data(idx);
            if ~isempty(this.model) && isequal(this.model, this.history(val))
                return; 
            end 
            this.changemodel(this.history(val));
            this.view.setEnabled('btnhistrename', true);
            this.view.setEnabled('btnhistdelete', true);
        end 
        
        function cbHistRename(this, src, ~, ~)
            import tensorlab.gui.common.*;
            lst = findall(this.view.fh, 'Tag', 'lsthistory');
            if isempty(get(lst, 'Value')), return; end
            newname = inputrenamedlg('session', this.model.name, 'Rename session');
            if ~isempty(newname)
                % Update current model and the history model
                this.model.name = newname;
                hist = this.findmodel();
                if ~isempty(this), 
                    hist.name = newname;
                    this.handleHistoryEvents();
                end
            end 
        end 
        
        function cbHistDelete(this, src, ~, ~)
            import tensorlab.gui.common.*;
            lst = findall(this.view.fh, 'Tag', 'lsthistory');
            if isempty(get(lst, 'Value')), return; end
            [~,ind] = this.findmodel();
            this.history(ind) = [];
            this.model.resetresults();
        end 
        
        %% Menu callbacks 
        function cbGenerateCode(this, ~, ~, ~)
            sug = this.model.getExportSuggestion();
            [file,path] = uiputfile('*.m','Export Code to file', [sug '.m']);
            if path ~= 0, this.model.saveCode(fullfile(path,file)); end
        end 
        
        function cbSaveSession(this, ~, ~, ~)
            sug = this.model.getExportSuggestion();
            [file, path] = uiputfile('*.mat','Export session', [sug '.mat']);
            if path ~= 0
                history = this.history;
                savedata = false;

                if strcmpi(history(1).path, 'Workspace'), default = 'Yes';
                else default = 'No'; end
                    
                answer = questdlg('Do you want to export the data as well?', 'Export session', ...
                                  'Yes', 'No', default);
                if strcmpi(answer, 'yes')
                    for k = 1:numel(history)
                        history(k).path = fullfile(path,file);
                    end 
                    savedata = true;
                end                 
                save(fullfile(path,file), 'history');
                if savedata
                    s = struct(history(1).varname, history(1).T);
                    save(fullfile(path,file), '-append', '-struct', 's');
                end 
            end
        end

        function cbLoadSession(this, ~, ~, ~)
            [file, path] = uigetfile('*.mat', 'Load session');
            if path ~= 0
                try 
                    tmp = load(fullfile(path,file), 'history');
                catch 
                    errordlg('The selected file does not contain a valid CPD GUI session.');
                    return;
                end 
                orig = copy(tmp.history(1));
                tmp.history(1)
                if strcmpi(tmp.history(1).path, 'Workspace')
                    tmp.history(1).loadFromWorkspace();
                else 
                    tmp.history(1).loadFromFile();
                end
                if isempty(tmp.history(1).T)
                    return;
                end 
                if ~isequal(tmp.history(1), orig)
                    errordlg(['The loaded tensor data does not correspond to the data in the ' ...
                              'session.']);
                    return;
                end 
                this.history = tmp.history;
                for k = 2:numel(this.history)
                    this.history(k).T = this.history(1).T;
                end 
                this.changemodel(this.history(1));
                this.updateallfields();
            end
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

