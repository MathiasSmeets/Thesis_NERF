function [clusters, cludepth, spk_binning, spiketimes] = neurons_spike_map(varargin)

% Plot spike distribution along time for each neuron, given its depth along
% the probe. Calculate also percentage of spikes during a specific range of
% time in which events occur (for example, during horridge learning)
%
%   Inputs:
%
%   
% Mattia D'Andola, May 2021
def_clusters = [];
def_include_mua = 0;
% def_maxdepth = 1200;
def_eventrange = 0;
def_evtname = '';
def_plt = 1;
def_binning = [];
def_depth_offset = 0;

%inputs control
parser = inputParser;
addRequired(parser, 'sorted_path', @ischar);
addParameter(parser, 'clusters', def_clusters);
addParameter(parser, 'include_mua', def_include_mua);
% addParameter(parser, 'max_depth', def_maxdepth);
addParameter(parser, 'add_events_range', def_eventrange);
addParameter(parser, 'events_name', def_evtname);
addParameter(parser, 'plot_results', def_plt);
addParameter(parser, 'time_binning', def_binning);
addParameter(parser, 'depth_offset', def_depth_offset); %last channel inside

parse(parser,varargin{:})
spk_path = parser.Results.sorted_path;
clusters = parser.Results.clusters;
% maxdepth = parser.Results.max_depth;
add_events = parser.Results.add_events_range;
evt_name = parser.Results.events_name;
include_mua = parser.Results.include_mua;
plt = parser.Results.plot_results;
tbins = parser.Results.time_binning;
Doffset = parser.Results.depth_offset;

%load spikes
spk = loadKSdir(spk_path);

if add_events == 1 && isempty(evt_name)
    warning('you need to specify the ''events_name'' option');
    return;
end

%find only SU or SU+MUA
if isempty(clusters)
    if include_mua
        clusters = spk.cids;
    else
        clusters = spk.cids(spk.cgs==2);
    end
end

%get cluster probe depth
depths = get_cluster_depth(spk_path, clusters);

for i = 1:length(clusters)
    spiketimes{i} = spk.st(spk.clu==clusters(i));
    cludepth(i) = depths.depth(depths.cluster==clusters(i));
end
%plot spike map
if plt
    figure;
    load(fullfile(spk_path,'chanMap.mat'), 'ycoords');
    colors = jet(length(clusters));
    axes('NextPlot','add')
    for i = 1:length(clusters)
        plot(spiketimes{i}', ones(1,length(spiketimes{i}))*cludepth(i), '.', 'color', colors(i,:), 'DisplayName',['Neu ', num2str(clusters(i))]);
    end
    ylim([0 ycoords(Doffset)]);
    g = get(gca, 'ytick');
    set(gca,'color', [.75 .75 .75], 'yticklabel', abs(g-ycoords(Doffset)))
    xlabel('Time [s]')
    ylabel('Depth [\mu m]')
    legend('show')
    
    if add_events
        load(fullfile(spk_path,'events.mat'),'events');
        evnames = {events.type};
        evtid = find(contains(evnames,evt_name));
        %plot events onset and offset
        plot(ones(1,2)*events(evtid).onsets(1),[0 1200],'--b', 'DisplayName',[evnames{evtid}, ' onset'])
        plot(ones(1,2)*events(evtid).offsets(end),[0 1200],'--k', 'DisplayName',[evnames{evtid}, ' offset'])
    end
end

%compute response ratio during different time bins
if ~isempty(tbins)
    if size(tbins,2) ~= 2
        disp('the option ''time_binning'' must be a matrix Nx2, where N is the number of time bins, and the columns are respectively onset and offset of time bins')
        return
    end
    spk_binning = zeros(length(clusters), size(tbins,1));
    for i = 1:length(clusters)
        for k = 1:size(tbins,1)
            spk_binning(i,k) = numel(find(spiketimes{i}>=tbins(k,1) & spiketimes{i}<tbins(k,2)))/length(spiketimes{i});
        end
    end
else
    spk_binning = [];
end
