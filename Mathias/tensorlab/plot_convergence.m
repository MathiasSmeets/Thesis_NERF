function plot_convergence(output, options)
%PLOT_CONVERGENCE short discription
%   long description
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/23   NV      Initial version
    

    defaults = getdefaults();
    if nargin < 2, options = struct; end

    % Create figure with top radio button pane and axes
    fh = gcf;
    set(fh, 'Visible', 'off');
    styles = tensorlab.gui.common.getDefaultStyles();

    ptop    = uipanel(fh, styles.place);
    pbottom = uipanel(fh, styles.place);
        
    if ~isfield(output, 'fval')
        ind    = structfun(@(f) isfield(f, 'fval'), output);
        names  = fieldnames(output);
        output = rmfield(output, names(~ind));
        names  = names(ind);
        N      = numel(fields(output));
        % manual label positioning, yay matlab.
        uicontrol(ptop, styles.label, 'String', 'Step ', 'position', [10 0 50 18]);
        % Create radio buttons
        bg = uibuttongroup(ptop, styles.buttongroup, 'Position', [55 0 50 30], 'SelectionChangedFcn', ...
                           @draw);
        rb = cell(1, numel(output));
        for k = 1:N
            rb{k} = uicontrol(bg,'Style', 'radiobutton', 'String', names{k}, 'UserData', names{k}, ...
                              'Units', 'normalized', 'Position', [(k-1)/N 0 1/N 0.8]);
        end 
        set(rb{1}, 'Value', 1);
    end 
    
    axes('parent', pbottom);
    
    movegui(fh, 'center');
    set(fh, 'Visible', 'on');
    tensorlab.gui.common.setminsize(fh);
    set(fh, 'SizeChangedFcn', @cbResize);
    cbResize();
    draw();
    
    function draw(src,evt)
        userdata = [];
        if nargin >= 2, userdata = get(evt.NewValue, 'UserData'); end 
        if isempty(userdata) 
            if isfield(output, 'fval')
                out = output;
            else 
                names = fieldnames(output);
                out = output.(names{1});
            end 
        else 
            out = output.(userdata);
        end 
        algorithmname = strsplit(out.Name, '_');
        algorithmname = algorithmname{end};

        fields = {'fval', 'relfval', 'relstep', 'measure'};
        props  = {'TolAbs', 'TolFun', 'TolX', 'TolMeasure'};
        
        cla; 
        lgd  = {}; % legend entries 
        lgdh = []; % legend handles
        cidx = 1;  % color index
        if out.iterations <= 25, style = '.-'; 
        else style = '-'; end
        
        % Plot curves and tolerances
        for k = 1:numel(fields)
            if ~isfield(out, fields{k}) || isempty(out.(fields{k})), continue; end
            v = out.(fields{k});
            x = out.iterations - numel(v)+1:out.iterations;
            x(isnan(v)) = [];
            v(isnan(v)) = [];
            % plot and legend
            lgdh(end+1) = semilogy(x.', v.', style);
            lgd{end+1} = fields{k};
            % tolerances
            if cidx == 1; hold on; end
            if isstruct(options) && isfield(options, props{k})
                tol = options.(props{k});
            elseif isfield(defaults.(algorithmname), props{k})
                tol = defaults.(algorithmname).(props{k});
            else 
                tol =  [];
            end 
            if ~isempty(tol) && isfinite(tol)
                set(gca,'ColorOrderIndex',cidx);
                semilogy([0 out.iterations],[tol tol],'--');
            end 
            cidx = cidx + 1;
        end 
        switch out.info
          case 1, txt = 'Relative function value tolerance reached.';
          case 2, txt = 'Relative step length tolerance reached.';
          case 3, txt = 'Maximum number of iterations reached.';
          case 4, txt = 'Absolute function value tolerance reached.';
          case 5, txt = 'Measure tolerance reached.';
        end 
        title(txt);
        legend(lgdh, lgd{:}, 'location', 'SouthWest');
        xlabel('Iteration k');
        drawnow;
    end 
    
    function cbResize(h,evt)
        if isfield(output, 'fval'), height = 0; 
        else, height = 30; end 
        pos = get(fh, 'position');
        set(ptop,     'position', [0, pos(4)-height, pos(3), height]);
        set(pbottom,  'position', [20, 0, pos(3)-20, pos(4) - height]);
        if exist('bg', 'var')
            set(bg,   'Position', [100, 0, pos(3)-50, height]);
        end 
    end

    function defaults = getdefaults()
        defaults.als.TolFun   = 1e-8;
        defaults.als.TolX     = 1e-6;
        defaults.als.TolAbs   = -inf;
        defaults.als.MaxIter  = 500;
        defaults.minf.TolFun  = 1e-12;
        defaults.minf.TolX    = 1e-8;
        defaults.minf.TolAbs  = -inf;
        defaults.minf.MaxIter = 500;
        defaults.nls.TolFun   = 1e-12;
        defaults.nls.TolX     = 1e-8;
        defaults.nls.TolAbs   = -inf;
        defaults.nls.MaxIter  = 200;
    end 
    
    
end
