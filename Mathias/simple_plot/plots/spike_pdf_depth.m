function [pdfs, cdfs, ampBins, depthBins, ypl] = spike_pdf_depth(varargin)

%plot spike PDF and CDF by depth
%
% required input:
%   - Data = can be a char or a cell of chars
% optional inputs:
%   - plot_results = boolean (default 1). If 1 results are plotted
%   - amp_binning = 1 value (in uV units); defines the binning step for the spike amplitude
%                   (default 10)
%   - depth_binning = 1 value (in um units); defines the binning step for the recording
%                   depth (default 10)
%   - depth_lims = array of 2 values [min max] in un units); defines the depth limits to be shown
%                   (default [0 1200])
%   - condition_names = cellarray; defines the title of the different
%                       plots. If empty, plots will benumbered sequentially
%                       ('condition n.1', ... , 'condition n. n')
%
%   N.B. to work it needs in the path the spike-master code from cortex-lab: https://github.com/cortex-lab/spikes
%
% Mattia D'Andola, Feb 2021

def_plt = 1;
def_ampB = 10;
def_depB = 10;
% def_deplim = [0 1200];
def_condnames = {};
def_Doffset = 140;
def_folder = {};

% Added by Simon 23/08
def_recording_definition = 'default';

checkdata = @(x) iscell(x) || ischar(x) || isstruct(x);
%inputs control
parser = inputParser;
addRequired(parser, 'data',checkdata)
addParameter(parser, 'plot_results', def_plt);
addParameter(parser, 'amp_binning', def_ampB);
addParameter(parser, 'depth_binning', def_depB);
% addParameter(parser, 'depth_lims', def_deplim);
addParameter(parser, 'condition_names', def_condnames);
addParameter(parser, 'depth_offset', def_Doffset);
addParameter(parser, 'folder', def_folder);

% Added by Simon 23/08
addParameter(parser, 'recording', def_recording_definition);

parse(parser,varargin{:})
data = parser.Results.data;
plt = parser.Results.plot_results;
ampB = parser.Results.amp_binning;
depB = parser.Results.depth_binning;
% deplims = parser.Results.depth_lims;
condnames = parser.Results.condition_names;
Doffset = parser.Results.depth_offset;
folder = parser.Results.folder;

% Added by Simon 23/08
recording_definition = parser.Results.recording;

%check inputs
if ischar(data)
    data = {data};
elseif isstruct(data) && isempty(folder)
    warning('If data input is a spike struc, "folder" optional input is needed')
    return;
end

if ~isempty(folder) && ischar(folder)
    folder = {folder};
end

%iterate on data
for i = 1:length(data)
    if isstruct(data)
        spk = data(i);
        load(fullfile(folder{i},'chanMap.mat'), 'ycoords');
    else
        spk = loadKSdir(data{i});
        load(fullfile(data{i},'chanMap.mat'), 'ycoords');
    end
    
    if  contains(recording_definition,'early_phase') ...
            | contains(recording_definition,'mid_phase') ...
            | contains(recording_definition,'late_phase')
        
        %It has to be SORT_HORRIDGE
        load(fullfile(data{i},'events.mat'),'events');
        evnames = {events.type};
        evtid = find(contains(evnames,'horridge'));
        tdist = [];
        for j = 2:length(events(evtid).onsets)
            tdist(j) = events(evtid).onsets(j) - events(evtid).onsets(j-1);
        end
        dstim = smoothdata(tdist);
        
        a = 0;
        for k=1:length(dstim)
            if dstim(k)<1
                idx(k)=1;
            else
                idx(k)=2;
                a = 1;
            end
        end
        
        if a == 0
            idmax = length(dstim);
        else
            idmax=find(idx==2,1,'first');
        end
        
        tbin = floor([events(evtid).onsets(1) events(evtid).onsets(idmax); events(evtid).onsets(idmax) events(evtid).offsets(end); events(evtid).offsets(end) (events(evtid).onsets(1)+600)])
        
        if tbin(3,2)<tbin(3,1)         
        tbin(3,2)=tbin(2,2);
        tbin(3,1)=tbin(3,2)-30;
        end
    
    end
    
    
    
    deplims = [0 ycoords(Doffset(i))];
    %create depth binning vector
    depthBins{i} = deplims(1):depB:deplims(2);
    
    idmax=[];
    idmin=[];
    
    
    %Find Cluster name that are SU
    spk.cids(find(spk.cgs==2));
    [tf,idx] = ismember(spk.spikeTemplates,spk.cids(find(spk.cgs==2)));
    spk.st=spk.st(tf==true);
    
    if contains(recording_definition,'before')
        %Include long recording
        %Only consider below 900s
        idx=spk.st>900;
        idmax=find(idx==1,1,'first');
        inter=1:idmax;
        idmin=0;
        
    elseif contains(recording_definition,'horridge')
        %Only for long recording
        %Only consider above 900 and below 1500
        idx_min=spk.st>900;
        idmin=find(idx_min==1,1,'first');
        idx_max=spk.st>1500;
        idmax=find(idx_max==1,1,'first');
        inter=idmin:idmax;
        
    elseif contains(recording_definition,'early_phase')
        %Only for long recording
        %Only consider above 900 and below 1500
        
        idx_min=spk.st>tbin(1,1);
        idmin=find(idx_min==1,1,'first');
        idx_max=spk.st>tbin(1,2);
        idmax=find(idx_max==1,1,'first');
        inter=idmin:idmax;
        
    elseif contains(recording_definition,'mid_phase')
        %Only for long recording
        %Only consider above 900 and below 1500
        
        idx_min=spk.st>tbin(2,1);
        idmin=find(idx_min==1,1,'first');
        idx_max=spk.st>tbin(2,2);
        idmax=find(idx_max==1,1,'first');
        inter=idmin:idmax;
        
        
    elseif contains(recording_definition,'late_phase')
        %Only for long recording
        %Only consider above 900 and below 1500
        
        idx_min=spk.st>tbin(3,1);
        idmin=find(idx_min==1,1,'first');
        idx_max=spk.st>tbin(3,2);
        idmax=find(idx_max==1,1,'first');
        inter=idmin:idmax;
        
    elseif contains(recording_definition,'after')
        %Only for long recording
        %Only consider above 1500
        idx_min=spk.st>1500;
        idmin=find(idx_min==1,1,'first');
        inter=idmin:length(idx_min);
        idmax=0;
        
    elseif contains(recording_definition,'default')
        %For short recordibg
        %Consider all the recording
        idx_min=spk.st>0;
        inter=1:length(idx_min);
        idmax=0;
        
    end
    
    if isempty(idmax) | isempty(idmin)
        idx_min=spk.st>0;
        inter=1:length(idx_min);
    end
    
    
    
    %spk.spikeTemplates(tf==true);
    %spk.tempScalingAmps(tf==true);
    
    [spikeAmps, spikeDepths] = templatePositionsAmplitudes(spk.temps, ...
        spk.winv, spk.ycoords, spk.spikeTemplates(tf==true), spk.tempScalingAmps(tf==true));
    
    %     [spikeAmps2, spikeDepths2] = templatePositionsAmplitudes(spk.temps, ...
    %         spk.winv, spk.ycoords, spk.spikeTemplates, spk.tempScalingAmps);
    
    
    
    %recordingDur = spk.st(end);
    recordingDur = spk.st(inter(end))-spk.st(inter(1));
    
    ampBins{i} = 0:ampB:min(max(spikeAmps),800);
    [pdfs{i}, cdfs{i}] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins{i}, depthBins{i}, recordingDur);
    if plt
        [h1, h2, h3, h4] = plot_spk_cdf_pdf(pdfs{i}, cdfs{i}, ampBins{i}, depthBins{i});
        g = get(h1, 'ytick');
        set(h1, 'yticklabel', abs(g-ycoords(Doffset(i))))
        g = get(h2, 'ytick');
        set(h2, 'yticklabel', abs(g-ycoords(Doffset(i))));
        g = get(h3, 'ytick');
        set(h3, 'yticklabel', abs(g-ycoords(Doffset(i))));
        hs = [h1 h2 h3 h4];
        if ~isempty(condnames)
            suptitle(condnames{i})
        else
            suptitle(sprintf('Condition n.%d', i))
        end
    end
    ypl = smoothdata(sum(cdfs{i},2),'gaussian',10);
end
end

function [h1 h2 h3 h4] = plot_spk_cdf_pdf(pdfs, cdfs, ampBins, depthBins)
f = figure('Position', [200 100 660 890]);
depthX = depthBins(1:end-1)+mean(diff(depthBins))/2;
ampX = ampBins(1:end-1)+mean(diff(ampBins))/2;

h1 = subplot(2,3,1);
h1.Position(2) = 0.35;
h1.Position(4) = 0.6;

daspect([1 .7 1])
imagesc(ampX, depthX, pdfs)
xlabel('spike amplitude (µV)');
ylabel('depth (µm)');
title('pdf');
set(gca, 'YDir', 'normal');
%makepretty

h2 = subplot(2,3,2);
h2.Position(2) = 0.35;
h2.Position(4) = 0.6;
daspect([1 .7 1])
imagesc(ampX, depthX, cdfs)
xlabel('spike amplitude (µV)');
ylabel('depth (µm)');
title('inverse cdf');
set(gca, 'YDir', 'normal');
colorbar
colormap(colormap_greyZero_blackred)
caxis([0 20]);
%makepretty
ch = get(f, 'Children');
chTypes = get(ch, 'Type');
cbar = ch(strcmp(chTypes, 'colorbar'));
cbar.Label.String = 'firing rate (sp/s)';

h3 = subplot(2,3,3);
h3.Position(4) = h2.Position(4);
h3.Position(2) = h2.Position(2);
ypl = smoothdata(sum(cdfs,2),'gaussian',10);
plot(ypl, depthX, '-k', 'linewidth',1)
set(gca, 'YDir', 'normal')
ylim([0 depthX(end)])

pos = [h2.Position(1), 0.15, h2.Position(3), 0.15];
h4 = subplot('Position',pos);
xpl = smoothdata(sum(cdfs),'gaussian',10);
plot(ampX, xpl, '-k', 'linewidth', 1)

end