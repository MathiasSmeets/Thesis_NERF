function [plot_options, f] = plot_rank1terms(U, varargin)
% plot_rank1_term short discription
%   long description
%

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Rob Zink            (Rob.Zink@esat.kuleuven.be)
%
% Version History:
% - 2018/04/26   RZ      Second version
% - 2018/04/25   NV      Initial version

    
    p = inputParser();
    p.addOptional('Title', '');
    p.addOptional('ModeSettings', []);
    p.addOptional('TermSettings', []);
    p.addOptional('ExportName', 'CPD');
    p.addOptional('MaxModes', 4);
    p.addOptional('MaxTerms', 5);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    
    N = numel(U);
    R = size(U{1},2);
    
    modesettings = p.Results.ModeSettings;
    if isempty(modesettings)
        names = arrayfun(@(n) sprintf('Mode %d', n), 1:N, 'UniformOutput', false);
        modesettings = struct('name', names, 'type', 'line');
    end 
    termsettings = p.Results.TermSettings;
    if isempty(termsettings)
        names = arrayfun(@(n) sprintf('Term %d', n), 1:R, 'UniformOutput', false);
        termsettings = struct('name', names, 'button', 0, 'panel', 0);
    end 
    
    nbcol = min(R, options.MaxTerms);
    nbrow = min(N, options.MaxModes);
    selectedterms = zeros(1,R);

    f = figure('Units', 'pixels', 'Position', [0, 0, 700, 500]);
    if isempty(options.Title), set(f, 'Name', 'CPD factors');
    else, set(f, 'Name', ['CPD factors - ' options.Title]); end 

    movegui(f, 'center');
    set(f, 'SizeChangedFcn', @cbResize);

    styles = tensorlab.gui.common.getDefaultStyles();
    
    pt      = uipanel(f, styles.panel);
    pb      = uipanel(f, styles.place);
    plegend = uipanel(f, styles.place);
    
    if R > options.MaxTerms
        psldR = uipanel(f, styles.place);
        sldR = uicontrol(psldR, styles.slider, 'min', 0, 'max', R-nbcol, 'Callback', @update, ...
                         'SliderStep', [1 1]/(R-nbcol));
    end
    if N > options.MaxModes
        psldN = uipanel(f, styles.place);
        sldN = uicontrol(psldN, styles.slider, 'min', 0, 'max', N-nbrow, 'Callback', @update, ...
                         'SliderStep', [1 1]/(N-nbrow), 'Value', (N-options.MaxModes)); 
    end 
    
    btnsettings = uicontrol(pt, styles.button, 'String', 'Plot settings',...
                            'Position', [3, 3, 100, 25], 'Callback', @cbPlotOptions);
    btnsavefigure = uicontrol(pt, styles.button, 'String', 'Save figure',...
                              'Position', [108, 3, 100, 25], ...
                              'Callback', @export);
    chkcollapsed = uicontrol(pt, styles.checkbox, 'String', 'Combine terms', 'Value', 0, ...
                             'position', [213 5 130, 20], 'Callback', @cbCollapse);
    chkyfixed = uicontrol(pt, styles.checkbox, 'String', 'Fix y axes', ...
                          'Position', [348, 5, 100, 20], 'Value', 0, ...
                          'Callback', @cbYFixed, 'Visible', 'On');
    txtselectedterms = uicontrol(pt, styles.edit, 'Callback', @cbSelectTerms, ...
                                 'Visible', 'off', 'Position', [468 5 100 20]);

    ax = [];
    createaxes(nbrow,nbcol);
    cbResize();
    update();
    
    function createaxes(nbrow,nbcol)
        delete(get(pb, 'Children'));
        ax = [];
        for n = 1:nbrow
            for r = 1:nbcol
                ax(n,r) = subplot(nbrow, nbcol, (n-1)*nbcol + r, 'parent', pb);
            end
        end
        if nbcol == 1
            for n = 1:nbrow
                pos = get(ax(n,1), 'Position');
                pos(3) = 0.95-pos(1);
                set(ax(n,1), 'Position', pos);
            end 
        end 
    end

    function cbResize(~,~)
        if R > options.MaxTerms && ~get(chkcollapsed, 'Value'), sldheight = 25;
        else, sldheight = 0; end 
        if N > options.MaxModes, sldwidth = 25;
        else, sldwidth = 0; end
        if get(chkcollapsed, 'Value'), lgdwidth = 150;
        else lgdwidth = 0; end 
        
        height = 30;
        pos = get(f, 'position');
        set(pt, 'position', [3, pos(4)-height, pos(3)-4, height]);
        set(pb, 'position', [0, sldheight, max(pos(3)-lgdwidth-sldwidth,1), max(pos(4)-height- ...
                                                          sldheight,1)]);
        tmp(1) = pos(3)-lgdwidth-sldwidth;
        tmp(2) = pos(4) - height - sldheight - 40 -  R*20 ;
        tmp(3) = max(lgdwidth,1);
        tmp(4) = R*20 + 10;
        set(plegend, 'Position', tmp);
        if R > options.MaxTerms 
            set(psldR, 'Position', [0 0 pos(3) sldheight]);
            set(sldR,  'Position', [5 5 pos(3)-5-sldwidth, max(sldheight-5,1)]);
        end 
        if N > options.MaxModes 
            set(psldN, 'Position', [pos(3)-sldwidth, sldheight, sldwidth, pos(4)-height-sldheight]);
            set(sldN,  'Position', [5, 5, max(sldwidth-5,1), pos(4)-height-sldheight-10]);
        end 
        set(txtselectedterms, 'Position', [468, 5, max(pos(3)-468-15,1), 20]);
    end

    function cbSelectTerms(src,evt)
        str = get(src,'String');
        if ~all(ismember(str, '0123456789+-.,: '))
            set(src, 'ForegroundColor', 'red');
            return;
        end 
        set(src, 'ForegroundColor', 'black');
        idx = unique(eval(sprintf('[%s]', str)));
        selectterms = false(1,R);
        selectterms(idx) = true;
    end

    function update(src,~)
        if nargin >= 1 && ~isempty(src)
            val = get(src, 'Value');
            if round(val) ~= val, set(src, 'Value', round(val)); end
        end 
        if R <= options.MaxTerms, termoffset = 0;
        else, termoffset = round(get(sldR, 'Value')); end
        
        if N > options.MaxModes
            modeoffset = max(N - options.MaxModes - round(get(sldN, 'Value')), 0);
        else 
            modeoffset = max(N - options.MaxModes, 0);
        end 
        
        if get(chkcollapsed, 'Value')
            nbrow = min(options.MaxModes, N);
            if nbcol ~= 1, createaxes(nbrow, 1); end 
            nbcol = 1;
        else  
            nbrow = min(options.MaxModes, N);
            tmp   = min(options.MaxTerms, R);
            if nbcol ~= tmp, createaxes(nbrow, tmp); end 
            nbcol = tmp;
        end            
        yfixed = get(chkyfixed, 'Value');
       
        for r = 1:nbcol
            if get(chkcollapsed, 'Value'), 
                idx = find(selectedterms);
            else 
                idx = min(r+termoffset,R);
            end 
            
            for n = 1:nbrow
                data = U{modeoffset+n}(:,idx);
                axes(ax(n,r));
                switch lower(modesettings(modeoffset+n).type)
                  case 'bar'
                    bar(data);
                    xlim([0.5 size(data,1)+0.5]);
                  case 'line'
                    plot(data);
                    xlim([1 size(data,1)]);
                  case 'dot' 
                    plot(data, 'o');
                    xlim([1 size(data,1)]);
                end
                % Fix colors
                if get(chkcollapsed, 'Value')
                    colors = get(gca, 'ColorOrder');
                    lines  = get(gca, 'Children');
                    lines  = lines(end:-1:1);
                    for k = 1:numel(idx)
                        c = colors(mod(selectedterms(idx(k))-1,size(colors,1))+1,:);
                        switch lower(modesettings(modeoffset+n).type)
                          case 'bar'
                            set(lines(k), 'FaceColor', c);
                          otherwise 
                            set(lines(k), 'Color', c);
                        end 
                    end 
                end 
                if yfixed
                    ylim([min(U{modeoffset+n}(:)), max(U{modeoffset+n}(:))]);
                end 
                if r == 1, 
                    ylabel(modesettings(modeoffset+n).name, 'FontWeight', 'bold'); 
                end 
                if n == nbrow
                    xlabel('Index');
                end 
                if n == 1 && ~get(chkcollapsed, 'Value')
                    title(termsettings(termoffset+r).name);
                end 
            end
        end
    end
    
    function drawlegend()
        import tensorlab.gui.common.*;
        if ~isempty(get(plegend, 'Children')), return; end
        set(plegend,'BorderType','line','HighlightColor',[0 0 0],'BackgroundColor','white');
        for r = 1:R
            p = uipanel(plegend, styles.place, 'Position', [0 (R-r)*20+4 145 25], 'BackgroundColor', 'white');
            h = axes(p, 'Visible', 'off', 'Units', 'pixels', 'Position', [5 7 25 10], ...
                     'Ytick', [], 'Xtick', []);
            hold(h, 'on');
            box(h, 'on');
            if selectedterms(r) > 0
                set(h, 'ColorOrderIndex', selectedterms(r));
                plot(h, [0, 1], [0.5 0.5], 'LineWidth', 2);
            else 
                plot(h, [0, 1], [0.5 0.5], 'LineWidth', 2, 'Color', 'white');
            end 
            set(h, 'Visible', 'off')
            h = uicontrol(p, styles.checkbox, 'String', termsettings(r).name, ...
                          'Position', [38 3 100 22], 'BackgroundColor', 'white', ...
                          'Callback', {@selectterm, r}, 'Tag', 'chkterm');
            termsettings(r).button = h;
            termsettings(r).panel  = p;
        end 
    end 

    function cbCollapse(h,evt)
        if get(chkcollapsed, 'Value'), 
            %set(txtselectedterms,'Visible','On');
            drawlegend();
        else 
            %set(txtselectedterms,'Visible','Off');
        end
        cbResize();
        update();
    end

    function export(src, evt)
        zipped_files = {};
        
        [extensions, zipped] = tensorlab.gui.common.inputfigexportdlg();
        
        if ~isempty(extensions)
            if ~iscell(extensions), extensions = {extensions}; end;
            if zipped, outfmt = '*.zip';
            else outfmt = '*.*'; end
            
            [name, path] = uiputfile(outfmt, 'Please enter the name for your figures', ...
                                     options.ExportName);
            if name ~= 0
                pos = get(f, 'Position');
                ptmp = uipanel(f, styles.place, 'BackgroundColor', 'white', ...
                               'Position', [0 0 pos(3:4)]);
                uistack(plegend, 'top');
                uistack(pb, 'top');
                if get(chkcollapsed, 'Value')
                    lgdwidth = 150;
                    tmp(1) = pos(3) - lgdwidth;
                    tmp(2) = pos(4) - 40 - R*20 ;
                    tmp(3) = lgdwidth;
                    tmp(4) = R*20 + 10;
                    set(plegend, 'Position', tmp);
                else
                    lgdwidth = 0;
                end 
                set(pb, 'Position', [0 0 pos(3)-lgdwidth, pos(4)]);

                bgcolor = get(pb, 'BackgroundColor');
                set(pb, 'BackgroundColor', 'white');
                set(plegend, 'BackgroundColor', 'white');
                drawnow();
                
                try 
                    [~,name,~]  = fileparts(name);
                    filename    = [path, name];                            
                    zipfilename = [filename, '.zip'];
                    
                    for k = 1:length(extensions)
                        if strcmp(extensions{k},'pdf')
                            set(f, 'PaperOrientation', 'landscape');
                            print(f, filename, '-dpdf', '-bestfit');
                        else
                            saveas(f, filename, extensions{k});
                        end
                    end
                    
                    if zipped
                        files = cellfun(@(ext) sprintf('%s.%s',filename,ext), extensions, ...
                                        'UniformOutput', false);
                        zip(zipfilename, files)
                        for k = 1:length(files), 
                            delete(files{k});
                        end
                    end
                end 
                
                delete(ptmp);                
                set(pb, 'BackgroundColor', bgcolor);
                set(plegend, 'BackgroundColor', bgcolor);
                cbResize();
            end
        end
    end

    function cbPlotOptions(src,evt)
        [ms,ts] = tensorlab.gui.common.plotsettingsdlg(modesettings, termsettings);
        if ~isempty(ms)
            modesettings = ms;
            termsettings = ts;
            
            if get(chkcollapsed, 'Value')
                for k = 1:numel(termsettings)
                    h = findall(termsettings(k).panel, 'Tag', 'chkterm');
                    set(h, 'String', termsettings(k).name);
                end 
            end 
            
        end 
        
        update();
    end

    function cbYFixed(src, evt)
        update();
    end

    function selectterm(src,evt,r)
        if get(termsettings(r).button, 'Value') == 1
            ind = histc(selectedterms(selectedterms > 0), 1:size(get(gca, 'ColorOrder'), 1));
            [~, selectedterms(r)] = min(ind);
        else 
            selectedterms(r) = 0;
        end
        colors = get(gca, 'ColorOrder');
        for k = 1:R
            h = findall(termsettings(k).panel, 'Type', 'Line');
            if selectedterms(k) > 0
                set(h, 'Color', colors(mod(selectedterms(k)-1,size(colors,1))+1,:));
            else 
                set(h, 'Color', 'white');
            end 
        end 
        update();
    end

end