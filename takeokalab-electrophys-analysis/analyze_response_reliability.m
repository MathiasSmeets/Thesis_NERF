function t = analyze_response_reliability(varargin)

% Analyze response to stimuli in terms of reliability: reliable and similar (latency, FR) response to
% all trials/repetitions.
%
%INPUTS REQUIRED:
%   1. path to sorting (kilosort) folder (char)
%   2. path to event file (char)
%   3. name of event type (as reported in event(n).type) (char)
%OPTIONAL INPUTS:
%   - 'include_MUA': check all SUs and MUs: default 0, i.e. checks only SU
%           (int boolean)
%   - 'window': time window to consider around the stim, in s: default
%           [-0.5 0.1] (1x2 vector of float/doubles)
%   - 'clusters': to be speficied if one wants to check only a subset of
%           clusters: default empty, i.e. checks all of them (1xN vector of int) 
%
% Mattia D'Andola, June 2021

def_muainc = 0;
def_win = [-.05 .1];
def_clusters = [];

%inputs control
parser = inputParser;
addRequired(parser, 'sorted_path', @ischar);
addRequired(parser, 'event_file', @ischar);
addRequired(parser, 'event_type', @ischar);
addParameter(parser, 'include_MUA', def_muainc);
addParameter(parser, 'window', def_win);
addParameter(parser, 'clusters', def_clusters);

parse(parser,varargin{:})
spk_path = parser.Results.sorted_path;
evt_path = parser.Results.event_file;
evt_type = parser.Results.event_type;
muainc = parser.Results.include_MUA;
window = parser.Results.window;
clusters = parser.Results.clusters;

%upload spike file
spk = loadKSdir(spk_path);

if isempty(clusters)
    if muainc
        clusters = spk.cids; %select all groups (SU + MU)
    else
        clusters = spk.cids(spk.cgs == 2); %select only SU
    end
end

%get the events onsets
events = load(evt_path);
events = events.events;
idev = find(strcmp({events.type}, evt_type));
if isempty(idev)
    warning(sprintf('no event type named %s',evt_type));
    return
end
evt_onsets = events(idev).onsets;

for i = 1:length(clusters)
    % select the spikes belonging to the unit i
    st = spk.st(spk.clu==clusters(i));
    
    %compute psth (get only binned array output and bins - bin size 1 ms)
    [~, bins, ~, ~, ~, ba] = psthAndBA(st, evt_onsets, window, 0.001);
    %smooth binned array with gaussian filter (window 5ms) and construct
    %psth
    smoothw = 5;
    gw = gausswin(round(smoothw*6),3);
    smWin = gw./sum(gw);
    baSm{i} = conv2(smWin, 1, ba', 'same')' ./ 0.001;
    psth(i,:) = mean(baSm{i});
    
    %normalize on baseline
    psth(i,:) = psth(i,:) / mean(psth(i,bins<-0.002),2);
    
    %find latency (ignore everything before 2ms past the stim to avoid artifacts)
    idx = find(bins>0.002);
    for j = 1:size(ba,1) %first on all trials
        ll = find(ba(j,idx)==1);
        if ~isempty(ll)
            latt(j) = bins(idx(ll(1)));
        else
            latt(j) = NaN;
        end
    end
    lati{i} = latt;    
    latmean = find(psth(i,bins>0.002) > (max(psth(i,bins>0.002)) * 0.2)); %then on the averaged psth
    if ~isempty(latmean)
        latency(i) = bins(idx(latmean(1)));
    else
        latency(i) = NaN;
    end
    
    %find response FR / baseline FR (ignore everything 2ms before and past the stim to avoid artifacts)
    for j = 1:size(baSm{i},1) %first on all trials
        ff(j) = max(baSm{i}(j, bins>0.002)) / mean(baSm{i}(j, bins<-0.002));
    end
    rsFRi{i} = ff;
    rsFR(i) = max(psth(i, bins>0.002)) / mean(psth(i, bins<-0.002)); %then on the averaged psth
    
    %check reliability
    noresp = numel(find(~isnan(lati{i})));
    reliability(i) = noresp/size(baSm{i},1);
end

%construct output table
t = table;
t.clusters = clusters';
for i = 1:length(clusters)
    t.all_responses{i} = baSm{i};
    t.psth(i,:) = psth(i,:);
end
for i = 1:length(clusters)
    t.latencies{i} = lati{i};
end
t.average_latency = latency';
for i = 1:length(clusters)
    t.FRs{i} = rsFRi{i};
end
t.average_FR = rsFR';
t.reliability = reliability';



