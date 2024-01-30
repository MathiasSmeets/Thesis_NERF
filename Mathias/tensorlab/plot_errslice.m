function output = plot_errslice(T, U, varargin)
%PLOT_ERRSLICE Plot the error in each slice.
%   TODO
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/06   NV      Initial version
    
    import tensorlab.auxiliary.*;
    
    if numel(varargin) >= 1 && isnumeric(varargin{1}) && iscellmat(U)
        U = {U,varargin{1}};
        varargin = varargin(2:end);
    end 
    abserr = false;
    if numel(varargin) >= 1 && ischar(varargin{1})
        abserr = strncmpi(varargin{1}, 'absolute', max(3, min(numel(varargin{1}), 8)));
        if abserr || strncmpi(varargin{1}, 'relative', max(3, min(numel(varargin{1}), 8)));
            varargin = varargin(2:end);
        end 
    end 
    
    p = inputParser;
    p.addOptional('NumberOfAxes', 3);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    T = ful(T);
    residual = abs(T-ful(U)).^2;
    T = abs(T).^2;
    sz = [1 getsize(T) 1];
    N  = getorder(T);
    nrm_slice = cell(1, N);
    err_slice = cell(1, N);
    for n = 1:N
        tmp = reshape(T, [prod(sz(1:n)) sz(n+1) prod(sz(n+2:end))]);
        nrm_slice{n} = sqrt(sum(sum(tmp,1),3)).';
        tmp = reshape(residual, [prod(sz(1:n)) sz(n+1) prod(sz(n+2:end))]);
        err_slice{n} = sqrt(sum(sum(tmp,1),3)).';
    end 
    
    useSlider = N > p.Results.NumberOfAxes;
    
    % Create figure with top radio button pane and axes
    f = gcf;
    set(f, 'Visible', 'off', ...
           'Units', 'pixels', ...
           'Position', [0 0 500 500]);
    ptop = uipanel('Units','pixels','Position',[0 0 10 10], 'Tag','Head', ...
                   'BorderType','none', 'Parent', f);
    pbottom = uipanel('Units','pixels','Position',[0 0 10 10], 'Tag','Head', ...
                      'BorderType','none', 'Parent', f);
    paxes = uipanel('Units','pixels','Position',[0 0 10 10], 'Tag','Head', ...
                    'BorderType','none', 'Parent', pbottom);
    
    % Create radio buttons
    bg = uibuttongroup(ptop, 'Visible','on',...
                       'Position',[0 0 1 1],...
                       'SelectionChangedFcn',@draw, 'BorderType', 'none');
    % manual label positioning, yay matlab.
    uicontrol(ptop,'style','text','String','Error: ', 'position', [10 -7 60 30]);
    
    r1 = uicontrol(bg,'Style',...
                   'radiobutton',...
                   'String','Absolute',...
                   'Tag', 'abs', ...
                   'Position', [70 0 100 30],...
                   'HandleVisibility','on');
    
    r2 = uicontrol(bg,'Style','radiobutton',...
                   'String','Relative',...
                   'Tag', 'rel', ...
                   'Position',[170 0 100 30],...
                   'HandleVisibility','on');
    
    % Create axes
    for n = 1:min(N,p.Results.NumberOfAxes)
        ax(n) = subplot(min(N,p.Results.NumberOfAxes), 1, n, 'parent', paxes);
    end

    if useSlider
        styles = tensorlab.gui.common.getDefaultStyles();
        pslider = uipanel('Units','pixels','Position',[0 0 10 10], 'Tag','Head', ...
                          'BorderType','none', 'Parent', pbottom);
        sldaxes = uicontrol(pslider, styles.slider);
        set(sldaxes, 'Min', 0, 'Max', N-p.Results.NumberOfAxes, 'Value', 0, ...
                     'SliderStep', [1 1]./(N-p.Results.NumberOfAxes));
        set(sldaxes, 'Value', N-p.Results.NumberOfAxes);
        set(sldaxes, 'CallBack', @cbModesChanged);
    end 

    
    movegui(f, 'center');
    set(f, 'Visible', 'on');
    tensorlab.gui.common.setminsize(f);
    set(f, 'ResizeFcn', @cbResize);
    cbResize();
    
    if abserr, set(r1, 'Value', 1); 
    else set(r2, 'Value', 1); end
    draw();
    
    function draw(~,evt)
        if nargin > 0
            abserr = strcmpi(get(evt.NewValue, 'Tag'), 'abs');
        end 
        if useSlider, offset = N-p.Results.NumberOfAxes-get(sldaxes, 'Value');
        else offset = 0; end
        for n = offset+(1:min(N, p.Results.NumberOfAxes));
            h = ax(n-offset);
            err = err_slice{n};
            if ~abserr, err = err./nrm_slice{n}; end
            plot(h, 1:numel(err), err, '.');
            xlim(h, [1 numel(err)]);
            xlabel(h, sprintf('mode-%d slice number', n));
            if abserr, ylabel(h, 'Abs. error');
            else, ylabel(h, 'Rel. error'); end 
        end 
    end 
    
    function cbResize(h,evt)
        height = 30;
        pos = get(f, 'position');
        set(ptop,    'position', [0, pos(4)-height, pos(3), height]);
        set(pbottom, 'position', [20, 0, pos(3)-20, pos(4) - height]);
        pos = get(pbottom, 'Position');
        if useSlider
            width  = 20;
            set(paxes, 'Position', [0,0,pos(3)-width,pos(4)]);
            set(pslider, 'Position', [pos(3)-width,0,width,pos(4)]);
            set(sldaxes, 'Position', [0,10,width,pos(4)-10]);
        else 
            set(paxes, 'Position', [0,0,pos(3),pos(4)]);
        end 
    end
    
    function cbModesChanged(src, ~, ~)
        val = get(src, 'Value');
        % prevent recursion
        if val == round(val), draw();
        else set(src, 'Value', round(val)); end 
    end 

end

