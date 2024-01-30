classdef MLSVDModel < handle

    properties (SetObservable)
        %% Data properties
        T
        path
        varname
        N
        type
        size_tens
        
        %% Computation settings
        size_core
        size_computed
        algorithm
        refinement
        
        %% Computation results
        U
        S
        sv
        ratio
        abserr
        relerr
        maxabserr
        maxabserrind
        maxrelerr
        maxrelerrind
        
        %% Export names
        exportU
        exportS
        exportT
    end 
    
    events
        Reset
    end 
    
    methods
        
        function this = MLSVDModel()
            this.reset();
        end

        function sug = getExportSuggestion(this)
            str = tensorlab.auxiliary.size2str(this.size_core, 'x', false);
            sug = sprintf('%s_%s_%s_mlr%s', this.varname, this.algorithm, this.refinement, str);
        end
   
        function set.T(this,T)
            this.T = T;
            this.N = getorder(T);
            this.size_tens = getsize(T); 
            this.type = getstructure(T);
            
            % Push tensor to workspace
            if ~isempty(T) && ~strcmpi('Workspace',this.path)
                assignin('base','T', this.T)
            end
        end

        function reset(this)
            this.path         = '';
            this.varname      = '';
            this.T            = [];
            this.N            = [];
            this.type         = '';
            this.size_tens    = [];
            this.size_core    = [];
            this.algorithm    = 'mlsvd_rsi';
            this.refinement   = 'none';
            this.U            = {};
            this.S            = [];
            this.ratio        = [];
            this.abserr       = [];
            this.relerr       = [];
            this.maxabserr    = [];
            this.maxabserrind = [];
            this.maxrelerr    = [];
            this.maxrelerrind = [];
            this.exportU      = 'Uc';
            this.exportS      = 'Sc';
            this.exportT      = 'Tc';
            
            notify(this, 'Reset');
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
                throw(e)
                errordlg(sprintf('The file %s is not found.', this.path));
                this.reset();
            end 
        end
        
        function loadFromWorkspace(this)
            try
                this.T = evalin('base',['eval(''' this.variable_name ''')']);
            catch
                errordlg(sprintf('The tensor %s is not found.', this.varname));
                this.reset();
            end
        end       
        
        function res = hasValidSizeCore(this)
            res = this.isValidSizeCore(this.size_core);
        end

        function res = isValidSizeCore(this, size_core)
            res = false;
            if isempty(size_core), return; end
            if numel(size_core) ~= this.N, return; end
            if any(~tensorlab.auxiliary.isnonnegint(size_core,true)), return; end 
            if any(size_core > this.size_tens), return; end
            [m,i] = max(size_core); 
            if m > prod(size_core([1:i-1 i+1:end])), return; end
            res = true;
        end 
        
        function compute(this)
            % Algorithm step 
            switch lower(this.algorithm)
              case 'mlsvd',     fun = @mlsvd; 
              case 'mlsvds',    fun = @mlsvds; 
              case 'lmlra_aca', fun = @lmlra_aca; 
              case 'mlsvd_rsi', fun = @mlsvd_rsi; 
            end 
            fprintf('Performing %s on %s... ', this.algorithm, this.varname);
            try 
                if strcmpi(this.algorithm, 'lmlra_aca')
                    [this.U, this.S, this.sv] = fun(this.T, this.size_core, 'FillCore', true);
                    this.sv = this.sv.sv;
                else 
                    [this.U, this.S, this.sv] = fun(this.T, this.size_core);
                end 
                fprintf('finished\r\n');
            catch e
                fprintf('error\r\n');
                rethrow(e)
            end 

            % Refinement step 
            switch lower(this.refinement)
              case 'nls',  fun = @lmlra_nls; 
              case 'minf', fun = @lmlra_minf; 
              otherwise,   fun = [];
            end 
            if ~isempty(fun) 
                try 
                    fprintf('Performing lmlra_%s on %s... \r\n', this.refinement, this.varname);
                    [this.U, this.S] = fun(this.T, this.U, this.S, 'Display', 10);
                    fprintf('Computation finished\r\n');
                catch e
                    fprintf('error\r\n');
                    rethrow(e)
                end                 
            end  
            
            this.size_computed = this.size_core;
            
            % Global statistics
            residual = lmlrares(this.T, this.U, this.S); 
            this.abserr = frob(residual);
            this.relerr = this.abserr/frob(this.T);
            this.ratio  = prod(this.size_tens) / (sum(cellfun(@numel, this.U)) + numel(this.S));
            
            % Maximal error 
            [this.maxabserr, ind] = max(abs(residual(:)));
            this.maxabserrind = this.ind2sub(ind);
            residual = abs(residual(:));
            if strcmpi(this.type, 'sparse')
                tmp = residual(this.T.ind);
                residual = inf(this.T.size);
                residual(this.T.ind) = tmp./this.T.val;
            else 
                residual = residual(:)./this.T(:);
            end 
            [this.maxrelerr, ind] = max(abs(residual));
            this.maxrelerrind = this.ind2sub(ind(1));
        end
        
        function sub = ind2sub(this, ind)
            sub = cell(1, this.N);
            [sub{:}] = ind2sub(this.size_tens, ind);
            sub = horzcat(sub{:});
        end 
        
        function cbExportToWorkspace(this)
            assignin('base', this.exportU,  this.U);
            assignin('base', this.exportS,  this.S);
            assignin('base', this.exportT, lmlragen(this.U, this.S));
            fprintf('Results saved as %s, %s and %s\r\n', ...
                    this.exportU, this.exportS, this.exportT);
        end
        
        function plotres(this)
            if ~isempty(this.U)
                f = figure('Name','Original data and low multilinear rank model');
                % set(f,'WindowStyle','modal');
                visualize({this.U, this.S},'original',this.T) % visualize difference
            end
        end
        
        function plot_errslice(this, errtype)
        % errtype can be 'absolute' or 'relative'
            figure();
            plot_errslice(this.T, this.U, this.S, errtype);
        end
        
        function saveCode(this, resultfile)
            import tensorlab.auxiliary.filltemplate;
            
            templatefile = fullfile(which('tensorlab.gui.mlsvd.mlsvd_template'));
            
            conditions = struct;
            conditions.refinement = ~strcmpi('none',this.refinement);
            conditions.fromfile   = ~strcmpi('Workspace',this.path);
                        
            replacements = struct;
            replacements.timestamp  = datestr(datetime);
            replacements.tensor     = this.varname;
            replacements.sizecore   = mat2str(this.size_core);
            replacements.algorithm  = this.algorithm;
            replacements.refinement = ['lmlra_', this.refinement];
            replacements.file       = strrep(strrep(this.path, '\', '\\'), '\\', '\\\\');
            replacements.var        = this.varname;
            replacements.U          = this.exportU;
            replacements.S          = this.exportS;
            replacements.Tc         = this.exportT;            

            success = filltemplate(templatefile, resultfile, replacements, conditions);
            
            if success, edit(resultfile); end
            
        end
    end
end
