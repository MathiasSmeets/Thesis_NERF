function algorithm = validate_algorithm(type, algorithm, default)
%VALIDATE_ALGORITHM Checks if the algorithm matches the optimization type.
%   ALGORITHM = VALIDATE_ALGORITHM(TYPE) checks if the given type is known by
%   Tensorlab and selects the appropriate algorithm. If a valid type is given ,
%   ALGORITHM contains a suitable optimization routine. 
%
%   ALGORITHM = VALIDATE_ALGORITHM(TYPE, ALGORITHM) also checks if the given
%   algorithm is valid for the given optimization type. If not, an error is
%   thrown. 
%   
%   ALGORITHM = VALIDATE_ALGORITHM(TYPE, ALGORITHM, DEFAULT) allows the
%   default algorithms for NLS and MINF to be set. DEFAULT should be a struct
%   containing two fields 'nls' and 'minf', each containing a suitable
%   function handle. If one of the fields is missing, the Tensorlab defaults
%   are chosen, which are @nls_gndl for NLS and @minf_lbfgsdl for MINF.
%   Alternatively, DEFAULT can be a function handle.
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/04/04   NV      Initial version
    
    if nargin == 1, algorithm = ''; end    

    % Known algorithms
    known = struct;
    known.nls  = {@nls_gndl,@nls_gncgs,@nls_lm,@nlsb_gndl};
    known.minf = {@minf_lbfgsdl,@minf_lbfgs,@minf_ncg,@minf_sr1cgs};
    known.nls  = cellfun(@func2str, known.nls, 'UniformOutput', false);
    known.minf = cellfun(@func2str, known.minf, 'UniformOutput', false);
    
    % Process type 
    if ~ischar(type)
        error('validate_algorithm:unknownOptimizationType', ...
              ['The optimization type should be a string (''nls'' or ' ...
               '''minf'').']);
    end
    type = lower(type);
    if ~any(strcmpi(type, {'nls', 'minf'}))
        error('validate_algorithm:unknownOptimizationType', ...
              ['The optimization type %s is unknown. Only ''nls'' and ''minf'' ' ...
               'are supported now'], type);
    end
    
    if strcmpi(type, 'nls'), othertype = 'minf'; 
    else othertype = 'nls'; end     
    
    % Process default argument 
    if nargin <= 2, default = struct; end
    if isa(default, 'function_handle')
        f = default;
        default = struct;
        default.(type) = f;
    elseif ~isstruct(default)
        error('validate_algorithm:default', ['Default should be a struct or a ' ...
        'function handle']); 
    end
    if ~isfield(default, 'nls'),  default.nls  = @nls_gndl; end
    if ~isfield(default, 'minf'), default.minf = @minf_lbfgsdl; end
    if ~all(cellfun(@(s) any(strcmp(s, {'nls','minf'})), ...
                    fieldnames(default)))
        error('validate_algorithm:default', ['Default can only contain a ''nls'' ' ...
                            'and/or a ''minf'' field (both lowercase).']);
    end

    % Process algorithm
    if ~isa(algorithm, 'function_handle')
        algorithm = default.(type);
    end 
    if ~any(strcmpi(getfunctionname(algorithm), known.(type)))
        if any(strcmpi(getfunctionname(algorithm), known.(othertype)))
            error('validate_algorithm:incompatibleAlgorithm', ...
                  ['The %s method is incompatible with the %s optimization ' ...
                   'type.'], getfunctionname(algorithm), type);
        else 
            warning('validate_algorithm:unknownAlgorithm', ...
                    'The %s method is not known by Tensorlab.', ...
                    getfunctionname(algorithm));
        end
    end
end

function name = getfunctionname(algorithm)
% Convert an (anonymous) function handle to a string
    name = func2str(algorithm);
    name = regexp(name, '(?:@\([^)]*\))?([^()]+)(?:\([^)]*\))?', 'tokens');
    name = name{1}{1};
end
