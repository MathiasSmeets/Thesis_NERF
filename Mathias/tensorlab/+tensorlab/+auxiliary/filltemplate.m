function success = filltemplate(templatefile, outfile, replacements, conditions)
%FILLTEMPLATE Create file using template. 
%   FILLTEMPLATE(TEMPLATEFILE, OUTFILE, REPLACEMENTS, CONDITIONS) creates a new file OUTFILE by
%   filling in the template defined in TEMPLATEFILE by replacing the placeholders defined in the
%   struct REPLACEMENTS, i.e., each occcurrence of ${name} in the template is replaced by
%   REPLACEMENTS.('name'). Conditional output is supported by the following syntax in template: 
%
%       $IF{condition} true code $ELSE false code $END
%   
%       $IF{condition}
%       true code  (multiline possible)
%       $ELSE 
%       false code (multiline possible)
%       $END   
%
%  The conditions are fields in CONDITIONS that evaluate to true or false, e.g.,
%  CONDITIONS.('name').  

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/08/06   NV      Initial version

    try 
        tid = fopen(templatefile, 'r');
        rid = fopen(outfile, 'w+');

        prevempty = false;
        while 1
            tline = fgetl(tid);
            % End of file
            if ~ischar(tline), break, end
            % Ignore additional empty line for conditional statements
            if prevempty && isempty(tline)
                prevempty = false;
                continue;
            end 
            prevempty = false;
            % Ignore template comment
            if regexp(tline, '^\s*%#'), continue; end
            % Escape some characters
            tline = regexprep(tline, '\\', '\\\\');
            % Check for inline conditionals
            pat = '\$IF{(?<cond>[a-zA-Z]+)}\s+(?<val>.+)\s+\$END';
            conds = regexp(tline, pat, 'names');
            for k = 1:numel(conds)
                rep = regexp(conds(k).val, '\s+\$ELSE\s+', 'split');
                if numel(rep) == 1, rep = {rep, ''}; end
                if isfield(conditions, conds(k).cond) && conditions.(conds(k).cond)
                    rep = rep{1};
                else 
                    rep = rep{2};
                end
                rep   = regexprep(rep, '\$', '\\\$');
                tline = regexprep(tline, pat, rep, 'once');
            end
            % Check for multiline conditionals
            matches = regexp(tline, '^\s*\$IF{(?<cond>[a-zA-Z]+)}\s*$', 'names');
            if ~isempty(matches)
                append = isfield(conditions, matches.cond) && conditions.(matches.cond);
                tline = '';
                while 1
                    tmp = fgetl(tid);
                    if ~ischar(tmp), break; end
                    tmp = regexprep(tmp, '\\', '\\\\');
                    if regexp(tmp, '^\s*\$ELSE\s*$'), append = ~append; continue, end
                    if regexp(tmp, '^\s*\$END\s*$'), break; end
                    if append && ~isempty(tline), tline = [tline '\r\n']; end
                    if append, tline = [tline tmp]; end
                end 
                if isempty(tline), prevempty = true; continue; end
            end 
            % Replace patterns
            matches = regexp(tline, '\${([0-9a-zA-Z-]+)}', 'match');
            matches = cellfun(@(s) s(3:end-1), matches, 'UniformOutput', false);
            for k = 1:numel(matches)
                if isfield(replacements, matches{k})
                    rep = replacements.(matches{k});
                else 
                    rep = ['UNKNOWN FIELD(' matches{k} ')'];
                end 
                tline = regexprep(tline, '\${[0-9a-zA-Z-]+}', rep, 'once');
            end 
            tline = regexprep(tline, '\$', '\\\$');
            tline = regexprep(tline, '%', '%%');
            fprintf(rid, [tline '\r\n']);
        end
        success = true;
    catch 
        success = false;
    end 
    
    fclose(rid);
    fclose(tid);
end 
