function [f, psdmean] = lfp_comparison_recordings(varargin)

def_chs = [];
def_filt = [];
def_fs = 30000;
def_labs = {};
def_plt = 1;

%inputs control
parser = inputParser;
addRequired(parser, 'data', @iscell)
addParameter(parser, 'select_channels',def_chs);
addParameter(parser, 'filter_signal', def_filt);
addParameter(parser, 'fs', def_fs);
addParameter(parser, 'label_names', def_labs);
addParameter(parser, 'plot_results', def_plt);

parse(parser,varargin{:})
data = parser.Results.data;
sel_chs = parser.Results.select_channels;
filt_sign = parser.Results.filter_signal;
label_names = parser.Results.label_names;
fs = parser.Results.fs;
plt = parser.Results.plot_results;

if isempty(sel_chs)
    sel_chs = 1:size(data{1},1);
end


%iterate on different data periods
sizemean = [];
sizemean_down = [];
for i = 1:length(data)
    %select channels
    templfp = data{i}(sel_chs,:);
    
    %filter %data
    if ~isempty(filt_sign)
        templfp = data_filter(templfp, 'bandpass', filt_sign, 3, 'fs', fs);
    end
%     templfp = resample(templfp,1,20,2);
    [fall{i},P] = psd(templfp, fs, 1:size(templfp,1), 'welch', 0);
    psdmean_all{i} = P;
    sizemean = [sizemean size(psdmean_all{i},2)];
    
%     templfp_down = downsample(templfp,300);
%     [fall_down{i},P_down] = psd(templfp_down, fs/300, 1:size(templfp_down,1), 'welch', 0);
%     psdmean_all_down{i} = P_down;
%     sizemean_down = [sizemean_down size(psdmean_all_down{i},2)];
end

[maxL, indmax] = max(sizemean);
% [maxL_down, indmax_down] = max(sizemean_down);
for i = 1:length(psdmean_all)
    if sizemean(i)<maxL
        psdmean(i,:) = interp1(fall{i}, psdmean_all{i}, fall{indmax});
    else
        psdmean(i,:) = psdmean_all{i};
    end
%     if sizemean_down(i)<maxL_down
%         psdmean_down(i,:) = interp1(fall_down{i}, psdmean_all_down{i}, fall_down{indmax_down});
%     else
%         psdmean_down(i,:) = psdmean_all_down{i};
%     end
end
f = fall{indmax};
% f_down = fall_down{indmax_down};

%<<<<<<< Updated upstream
% for i = 1:size(psdmean,1)
%     f_down = linspace(f(),f(end),round(length(f)/30));
%     psdmean_down(i,:) = interp1(f, psdmean(i,:), f_down);
% end
%=======
% colors = colormap(jet(size(psdmean,1)));
%colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
%>>>>>>> Stashed changes

if plt
    %colors = colormap(parula(size(psdmean,1)));
    colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
    figure, hold on
    if isempty(label_names)
        flag = 1;
    else 
        flag = 0;
    end
    for i = 1:size(psdmean,1)
        %normalize for representation
        psdmean_norm = psdmean(i,:)./max(psdmean(i,:));
        f_down = linspace(f(1),f(end),round(length(f)/100));
    %     psdmean_norm_down = interp1(f, psdmean_norm, f_down, 'nearest');
        flim = find(f<=110);
        psdmean_norm_smooth = smooth(f(1:flim(end)),psdmean_norm(1:flim(end)),0.01,'rloess');
        plot(f,10*log10(psdmean_norm), 'color', piuchiaropiuscuro(colors(i,:),.5));
        plot(f(1:flim(end)),10*log10(psdmean_norm_smooth),'color',colors(i,:),'LineWidth',4);
        if flag
            label_names{i} = sprintf('Condition %d',i);
        end
    end
    xlabel('f [Hz]')
    ylabel('PSD [db/Hz]')
    xlim([0 100])
    legend(label_names)
end

