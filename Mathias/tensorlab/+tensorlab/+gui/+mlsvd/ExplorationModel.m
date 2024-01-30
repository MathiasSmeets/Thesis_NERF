classdef ExplorationModel < handle

    properties (SetObservable)
        %% Data properties
        T
        varname
        N
        type
        size_tens
        expanded
        
        %% Computation settings
        size_upper
        size_core
        size_computed
        algorithm
        
        %% Computation results
        U
        S
        sv
        ratio
        abserr
        relerr
        
        %% Plot options
        yscale
        plotfun
    end 
    
    methods
        
        function this = ExplorationModel(T, name)
            this.algorithm = 'mlsvd_rsi';
            this.expanded  = false;
            if nargin > 0
                if nargin < 2, name = 'unknown'; end 
                this.T       = T;
                this.varname = name;
            end 
        end
        
        function set.T(this,T)
            this.T         = T;
            this.N         = getorder(T);
            this.size_tens = getsize(T); 
            this.type      = getstructure(T);
            if ~any(strcmpi(this.type, {'full', 'sparse'}))
                if prod(this.size_tens) <= 1e6
                    this.T = ful(T);
                    this.expanded = true;
                else 
                    this.algorithm = 'lmlra_aca';
                end 
            end 
        end
        
        function algorithms = getSupportedAlgorithms(this)
            if ~any(strcmpi(this.type, {'full', 'sparse'})) && ~this.expanded
                algorithms = {'lmlra_aca'};
            else 
                algorithms = {'mlsvd_rsi', 'lmlra_aca', 'mlsvd'};
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
        
    end
end