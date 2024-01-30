classdef CPDModel < handle & matlab.mixin.Copyable

    properties (SetObservable)
        %% Data properties
        T
        path
        varname
        N
        type
        size_tens
        isreal
        
        %% Initialization settings
        initmethod
        initalgorithm
        Uinitname
        Uinitmanual
        Uinit
        R
        Rmanual
        
        %% Computation settings
        algorithm
        algorithmoptions
        symmetrysettings
        nonnegmodes
        plotconvergence
        
        %% Computation results
        Ures
        output
        corcondia
        conditionnumber
        ratiosv
        abserr
        normrank1
        congruence
        relerr

        %% Export names
        exportU
        
        %% Model name
        name
        
        %% Plot options
        modesettings
        termsettings
    end 
    
    events
        Reset
    end 
    
    methods
        
        function this = CPDModel()
            this.reset();
        end
        
        function sug = getExportSuggestion(this)
            if strcmpi(this.initmethod, 'manual')
                init = 'man';
            else 
                switch this.initalgorithm
                  case 'cpd_gevd', init = 'gevd';
                  otherwise, init = this.initalgorithm;
                end 
            end 
            switch this.algorithm
              case 'cpd_nls', algo = 'nls';
              case 'cpd_als', algo = 'als';
            end 
            sug = sprintf('cpd_%s_R%d_%s_%s', this.varname, this.R, init, algo);
        end 

        function sug = getSessionSuggestion(this)
            if strcmpi(this.initmethod, 'manual')
                init = 'man';
            else 
                switch this.initalgorithm
                  case 'cpd_gevd', init = 'gevd';
                  otherwise, init = this.initalgorithm;
                end 
            end 
            switch this.algorithm
              case 'cpd_nls', algo = 'nls';
              case 'cpd_als', algo = 'als';
            end 
            sug = sprintf('%s_R%d_%s_%s', this.varname, this.R, init, algo);
        end
   
        function set.T(this,T)
            this.T = T;
            this.N = getorder(T);
            this.size_tens = getsize(T); 
            this.type = getstructure(T);
            this.isreal = isempty(T) || isreal(ful(T,1));
            
            % Push tensor to workspace
            if ~isempty(T) && ~strcmpi('Workspace',this.path)
                assignin('base','T', this.T)
            end
        end

        function reset(this)
            this.path             = '';
            this.varname          = '';
            this.T                = [];
            this.N                = [];
            this.type             = '';
            this.size_tens        = [];
            this.algorithm        = 'cpd_nls';
            this.algorithmoptions = this.getDefaultAlgorithmOptions();
            this.symmetrysettings = [];
            this.nonnegmodes      = [];
            this.plotconvergence  = false;
            this.initmethod       = 'auto';
            this.initalgorithm    = 'cpd_gevd';
            this.Uinitname        = [];
            this.Uinitmanual      = [];
            this.R                = [];
            this.Rmanual          = [];
            this.modesettings     = [];
            this.termsettings     = [];
            this.resetresults();
            
            notify(this, 'Reset');
        end 
        
        function resetresults(this)
            this.Uinit            = {};
            this.Ures             = {};
            this.output           = [];
            this.relerr           = [];
            this.corcondia        = [];
            this.conditionnumber  = [];
            this.ratiosv          = [];
            this.normrank1        = [];
            this.congruence       = [];
            this.exportU          = 'Ures';            
        end 
        
        function loadFromFile(this)
            try 
                vars = load(this.path);
                try 
                    this.T = vars.(this.varname);
                catch 
                    errordlg(sprintf('The variable %s is not found in %s.', this.varname, this.path));
                    this.reset();
                end 
            catch e
                errordlg(sprintf('The file %s is not found.', this.path));
                this.reset();
            end 
        end
        
        function loadFromWorkspace(this)
            try
                this.T = evalin('base', ['eval(''' this.varname ''')']);
            catch
                errordlg(sprintf('The tensor %s is not found.', this.varname));
                this.reset();
            end
        end       
                
        function compute(this)
            options = struct;
            switch lower(this.initmethod)
              case 'manual'
                Uinit = this.Uinitmanual;
                this.R = size(Uinit{1},2);
              case 'auto'
                Uinit = this.Rmanual;
                this.R = this.Rmanual;
                switch lower(this.initalgorithm)
                  case 'cpd_gevd'
                    options.Initialization = @cpd_gevd;
                  case 'randn'
                    options.Initialization = @cpd_rnd;
                    options.InitializationOptions = struct('Real', @randn);
                  case 'rand'
                    options.Initialization = @cpd_rnd;
                    options.InitializationOptions = struct('Real', @rand);
                end 
            end 
            
            switch lower(this.algorithm)
              case 'cpd_nls'
                options.Algorithm = @cpd_nls;
              case 'cpd_als'
                options.Algorithm = @cpd_als;
            end 
            options.AlgorithmOptions = this.algorithmoptions.(this.algorithm);
            if this.plotconvergence
                figure;
                options.AlgorithmOptions.ShowCurves = this.plotconvergence;
            end 
            
            if ~isempty(this.nonnegmodes) && this.supportsNonnegativity()
                options.NonnegativeModes = this.nonnegmodes;
            end 
            
            if ~isempty(this.symmetrysettings) && this.supportsSymmetry()
                options.VariablesInMode = this.symmetrysettings;
            end 
            
            options.Display = 10;
            
            [Ures,out] = cpd(this.T, Uinit, options);
            
            % Sort terms
            nrm = cellfun(@(u) sqrt(sum(abs(u).^2,1)), Ures, 'UniformOutput', false);
            nrm = prod(vertcat(nrm{:}),1);
            [this.normrank1, ind] = sort(nrm, 2, 'descend');
            Ures = cellfun(@(u) u(:,ind), Ures, 'UniformOutput', false);
            
            % Compute angles
            nrm  = cellfun(@(u) sqrt(sum(abs(u).^2,1)), Ures, 'UniformOutput', false);
            Utmp = cellfun(@(u,n) bsxfun(@rdivide, u, n), Ures, nrm, 'UniformOutput', false);
            W = cellfun(@(u) real(u'*u), Utmp, 'UniformOutput', false);
            this.congruence = prod(cat(3,W{:}),3);

            this.Ures = Ures;
            this.output = out;
            
            this.corcondia = corcondia(this.T, Ures);
            if isfield(out.Refinement, 'fval')
                this.relerr = sqrt(2*out.Refinement.fval(end))/frob(this.T);
            else 
                this.relerr = sqrt(2*out.Algorithm.fval(end))/frob(this.T);
            end 
            UHU = cellfun(@(u) u'*u, Ures, 'UniformOutput', false);
            [~,D] = eig(prod(cat(3,UHU{:}),3));
            D = sqrt(abs(diag(D)));
            this.ratiosv = max(D)/min(D);
            
            sz = cellfun('size', Ures, 1);
            if this.isConditionAutoComputed()
                this.conditionnumber = cpd_cond(Ures);
            else 
                this.conditionnumber = [];
            end 
        end
        
        function res = isConditionAutoComputed(this)
            res = (isempty(this.size_tens) || isempty(this.R));
            res = res || sum(min(this.size_tens,this.R)*this.R) < 3000;
        end 
        
        function computeCondition(this)
            this.conditionnumber = cpd_cond(this.Ures);
        end 
        
        function sub = ind2sub(this, ind)
            sub = cell(1, this.N);
            [sub{:}] = ind2sub(this.size_tens, ind);
            sub = horzcat(sub{:});
        end 
        
        function cbExportToWorkspace(this)
            assignin('base', this.exportU,  this.Ures);
            fprintf('Results saved as %s\r\n', this.exportU);
        end
        
        function plotres(this)
            if ~isempty(this.Ures)
                figure('name', sprintf('Low multilinear rank model versus original tensor - %s', this.name));
                visualize(this.Ures, 'original', this.T) 
            end
        end
        
        function ploterrslice(this, errtype)
        % errtype can be 'absolute' or 'relative'
            if nargin < 2, errtype = 'relative'; end
            figure('name', sprintf('Errors per slice - %s', this.name));
            plot_errslice(this.T, this.Ures, errtype);
        end
        
        function isvalid = isValidForGEVD(this, R)
            if nargin < 2, R = this.Rmanual; end
            if isempty(R) || isempty(this.N), isvalid = true; return; end
            size_tens = this.size_tens;
            isvalid = this.N > 2 && sum(size_tens>=R) >= 2 && sum(size_tens>=2) >= 3;
        end 
        
        function isvalid = hasValidInitialization(this)
            import tensorlab.auxiliary.*;
            if strcmpi(this.initmethod, 'manual')
                U = this.Uinitmanual;
                isvalid = ~isempty(U); 
                isvalid = isvalid && strcmpi(getstructure(U), 'cpd');
                isvalid = isvalid && isvalidtensor(U);
                isvalid = isvalid && sizeequals(this.size_tens, getsize(U));
            else 
                R = this.Rmanual;
                isvalid = ~isempty(R) && isnonnegint(R, true);
                isvalid = isvalid && (~strcmpi(this.initalgorithm, 'cpd_gevd') || this.isValidForGEVD(R));
            end 
        end 
        
        function isvalid = hasValidSettings(this)
            isvalid = true;
        end
        
        function [isvalid,msg] = isValidSymmetrySetting(this, setting)
            import tensorlab.auxiliary.*;
            msg = '';
            isvalid =  true;
            if isempty(this.T) || isempty(setting), return; end
            isvalid = false;
            if numel(setting) < this.N, 
                msg = sprintf(['Either %d (the order) or %d modes should be selected. In the case ' ...
                               'mode %d is selected, an additional mode with dimension 1 is added.'], ...
                               this.N, this.N+1, this.N+1);
                return; 
            end
            if ~all(isnonnegint(setting,true)), 
                msg = sprintf('The selected modes should be nonnegative integers.');
                return; 
            end
            if strcmpi('initmethod', 'manual')
                if max(setting) > this.N + 1, 
                    msg = sprintf(['The selected modes should be nonnegative integers ' ...
                                   'between 1 and %d.'], numel(this.Uinitmanual) + 1);
                    return;
                end
                sz = [cellfun('size', this.Uinitmanual, 1) 1];
                if ~sizeequals(sz(setting), this.size_tens)
                    msg = sprintf(['The size of the selected symmetric CPD (%s) does not match the ' ...
                                   'size of the data (%s)'], size2str(sz(setting), 'x', false), ...
                                  size2str(this.size_tens, 'x', false));
                    return; 
                end
            else 
                if max(setting) > this.N+1, 
                    msg = sprintf(['The selected modes should be nonnegative integers ' ...
                                   'between 1 and %d.'], this.N + 1);
                    return;
                end
                sz = [this.size_tens 1];
                if ~sizeequals(sz(setting), this.size_tens),
                    msg = sprintf(['The size of the selected symmetric CPD (%s) does not match the ' ...
                                   'size of the data (%s)'], size2str(sz(setting), 'x', false), ...
                                  size2str(this.size_tens, 'x', false))
                    return;
                end
            end 
            if ~isempty(this.nonnegmodes)
                tmp = zeros(1, numel(setting));
                tmp(this.nonnegmodes) = 1;
                counts = accumarray(setting(:), tmp(:), [max(setting) 1]).';
                if any(counts > 0 & counts < histc(setting,1:max(setting)))
                    ind = find(counts > 0 & counts < histc(setting,1:max(setting)), 1);
                    modes = find(setting == ind);
                    nnmodes = modes(find(tmp(modes)));
                    modestr = regexprep(size2str(modes, ', ', false), '(.*), ', '$1 and ');
                    if numel(nnmodes) == 1
                        msg = sprintf(['Modes %s are identical through symmetry, but only mode %d ' ...
                                       'is nonnegative.'], modestr, nnmodes(1));
                    else 
                        nnmodestr = regexprep(size2str(nnmodes, ', ', false), '(.*), ', '$1 and ');
                        msg = sprintf(['The modes %s are identical through symmetry, but only ' ...
                                       'modes %s are nonnegative.'], modestr, nnmodestr);
                    end 
                    return; 
                end
            end 
            if strcmpi(this.algorithm, 'cpd_als')
                msg = 'ALS currently does not support symmetry.';
                return;
            end 
            isvalid = true;
        end 
        
        function [isvalid, msg] = isValidNonnegativeMode(this, nonnegmodes)
            import tensorlab.auxiliary.*;
            msg = '';
            if isempty(this.T) || isempty(nonnegmodes), isvalid = true; return; end
            isvalid = false;
            if numel(nonnegmodes) > this.N || max(nonnegmodes) > this.N
                msg = sprintf('Only modes between 1 and %d (the order) can be selected', ...
                              this.N);
                return;
            end
            if ~all(isnonnegint(nonnegmodes, true))
                msg = sprintf(['The selected nonnegative modes should be nonnegative integers ' ...
                               'between 1 and %d.'], this.N);
                return;
            end
            if any(diff(sort(nonnegmodes))==0)
                msg = sprintf('Each mode can only appear once.');
                return;
            end
            settings = this.symmetrysettings;
            if ~isempty(settings)
                tmp = zeros(1, numel(settings));
                tmp(nonnegmodes) = 1;
                counts = accumarray(settings(:), tmp(:), [max(settings) 1]).';
                if any(counts > 0 & counts < histc(settings,1:max(settings)))
                    ind = find(counts > 0 & counts < histc(settings,1:max(settings)), 1);
                    modes = find(settings == ind)
                    nnmodes = modes(find(tmp(modes)))
                    modestr = regexprep(size2str(modes, ', ', false), '(.*), ', '$1 and ');
                    if numel(nnmodes) == 1
                        msg = sprintf(['Modes %s are identical through symmetry, but only mode %d ' ...
                                       'is nonnegative.'], modestr, nnmodes(1));
                    else 
                        nnmodestr = regexprep(size2str(nnmodes, ', ', false), '(.*), ', '$1 and ');
                        msg = sprintf(['The modes %s are identical through symmetry, but only ' ...
                                       'modes %s are nonnegative.'], modestr, nnmodestr);
                    end 
                    return; 
                end
            end 
            isvalid = true;
        end 
        
        function isvalid = isComputable(this)
            isvalid = this.hasValidInitialization() && ...
                      this.isValidSymmetrySetting(this.symmetrysettings) && ...
                      this.isValidNonnegativeMode(this.nonnegmodes);
        end 
        
        function isvalid = supportsSymmetry(this)
            isvalid = strcmpi(this.algorithm, 'cpd_nls');
        end 

        function isvalid = supportsNonnegativity(this)
            isvalid = this.isreal && strcmpi(this.algorithm, 'cpd_nls');
        end 

        function saveCode(this, resultfile)
            import tensorlab.auxiliary.*;

            templatefile = fullfile(which('tensorlab.gui.cpd.cpd_template'));
            
            conditions = struct;
            conditions.fromfile    = ~strcmpi('Workspace',this.path);
            conditions.showcurves  = this.plotconvergence;
            conditions.usenonneg   = ~isempty(this.nonnegmodes) && this.supportsNonnegativity();
            conditions.usesymmetry = ~isempty(this.symmetrysettings) && this.supportsSymmetry();
            conditions.autoinit    = strcmpi('auto', this.initmethod);
            conditions.randinit    = strcmpi('auto', this.initmethod) && ...
                ~strcmpi('cpd_gevd', this.initalgorithm);
            
            replacements           = struct;
            replacements.timestamp = datestr(datetime);
            replacements.tensor    = this.varname;
            replacements.algorithm = ['@', this.algorithm];
            replacements.TolFun    = num2str(this.algorithmoptions.(this.algorithm).TolFun);
            replacements.TolX      = num2str(this.algorithmoptions.(this.algorithm).TolX);
            replacements.MaxIter   = num2str(this.algorithmoptions.(this.algorithm).MaxIter);
            replacements.file      = strrep(strrep(this.path, '\', '\\'), '\\', '\\\\');
            replacements.var       = this.varname;
            replacements.U         = this.exportU;
            replacements.dist      = ['@', this.initalgorithm];
            replacements.modelname = this.name;
            replacements.nonnegmodes = size2str(this.nonnegmodes, ' ');
            replacements.varinmode = size2str(this.symmetrysettings, ' ');
            
            if strcmpi(this.initmethod, 'auto') 
                replacements.Uinit = 'R';
                replacements.R     = num2str(this.R);
            else 
                replacements.Uinit = this.Uinitname; 
                replacements.R     = sprintf('size(%s{1},2)', this.Uinitname);
            end 
            if strcmpi(this.initalgorithm, '@cpd_gevd') 
                replacements.initalgorithm = '@cpd_gevd'; 
            else 
                replacements.initalgorithm = '@cpd_rnd'; 
            end 

            success = filltemplate(templatefile, resultfile, replacements, conditions)
            
            if success, edit(resultfile); end
        end
        
        function opt = getDefaultAlgorithmOptions(this, algo)
            opt = struct;
            opt.cpd_nls = struct('TolFun', 1e-12, 'TolX', 1e-8, 'MaxIter', 200, 'CGMaxIter', 15);
            opt.cpd_als = struct('TolFun', 1e-12, 'TolX', 1e-8, 'MaxIter', 200, 'CGMaxIter', 15);
            if nargin >= 2
                if isfield(opt, algo)
                    opt = opt.(algo);
                else 
                    error('Unknown algorithm ''%s''.', algo);
                end 
            end 
        end

        function res = isequal(this, that)
            fields = {'path', 'varname', 'type', 'N', 'size_tens', 'initmethod', 'algorithm', ...
                      'algorithmoptions'};
            res = false;
            if ~isvalid(that), return; end
            for k = 1:numel(fields)
                f1 = this.(fields{k});
                f2 = that.(fields{k});
                if ~strcmpi(class(f1), class(f2)), return; end
                if ~isequal(f1, f2), return; end
            end
            if this.supportsSymmetry() && ~isequal(this.symmetrysettings, that.symmetrysettings)
                return; 
            end 
            if this.supportsNonnegativity() && ~isequal(this.nonnegmodes, that.nonnegmodes)
                return; 
            end 
            if strcmpi(this.initmethod, 'manual')
                if ~isequal(this.Uinitname, that.Uinitname), return; end
                if ~isequal(this.R, that.R), return; end
            else                     
                if ~isequal(this.Rmanual, that.Rmanual), return; end
                if ~isequal(this.initalgorithm, that.initalgorithm), return; end
            end 
            res = true;
        end 
        
        function s = saveobj(this)
            fields = properties(this);
            fields(strcmpi(fields, 'T')) = [];
            s = struct();
            for k = 1:numel(fields)
                s.(fields{k}) = this.(fields{k});
            end 
        end 
        
    end
    
    methods (Static)
        
        function this = loadobj(s)
            this = tensorlab.gui.cpd.CPDModel();
            fields = fieldnames(s);
            for k = 1:numel(fields)
                if ~isprop(this, fields{k}), continue; end
                this.(fields{k}) = s.(fields{k});
            end 
        end 
        
    end 
    
end
