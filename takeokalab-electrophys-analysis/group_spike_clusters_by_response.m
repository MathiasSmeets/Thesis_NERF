function [clusters, psth, bins, latency, offi, rsFR] = group_spike_clusters_by_response(varargin)

% Group slike clusters based on their psth 
%
% INPUTS:
%   required:
%       1. path to sorting foldesr
%       2. path to event file
%       3. event_type = type of events, as reported in event.type
%   options:
%       - 'recording_system': 'SpikeGLX' (default), 'Plexon', 'OpenEphys'
%       - 'nch': number of channels (default 385 - Neuropixels 3A)
%       - 'include_MUA': boolean, include/exclude MUA clusters (default 0)
%       - 'window': window around the stim (default [-.1 .2]);
%
% Mattia D'Andola
% Apr 2021


def_nch = 385;
def_muainc = 0;
def_win = [-.05 .1];
def_clusters = [];

%inputs control
parser = inputParser;
addRequired(parser, 'sorted_path', @ischar);
addRequired(parser, 'event_file', @ischar);
addRequired(parser, 'event_type', @ischar);
addParameter(parser, 'nch', def_nch);
addParameter(parser, 'include_MUA', def_muainc);
addParameter(parser, 'window', def_win);
addParameter(parser, 'clusters', def_clusters);

parse(parser,varargin{:})
spk_path = parser.Results.sorted_path;
evt_path = parser.Results.event_file;
evt_type = parser.Results.event_type;
nch = parser.Results.nch;
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
    baSm = conv2(smWin, 1, ba', 'same')' ./ 0.001;
    psth(i,:) = mean(baSm);
    
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
%     lati{i} = latt; 
    latmean = find(psth(i,bins>0.002) > (max(psth(i,bins>0.002)) * 0.2)); %then on the averaged psth
    if ~isempty(latmean)
        latency(i) = bins(idx(latmean(1)));
    else
        latency(i) = NaN;
    end
    
    %get response offset
    idxoff = find(diff(latt)~=1);
    if isempty(idxoff)
        idxoff(1) = length(latit);
    end
    try
        offi(i) = bins(idx(latt(idxoff(1))));
    catch
        offi(i) = NaN;
    end
    %find response FR / baseline FR (ignore everything 2ms before and past the stim to avoid artifacts)
    rsFR(i) = max(psth(i, bins>0.002)) / mean(psth(i, bins<-0.002));
    
%     %construct a binarized response (0 -> baseline, [+1, -1] -> [>, <] baseline+/-std 
%     binpos = find(psth(i, bins>0.002) > (1+std(psth(i,bins<-0.002))));
%     binneg = find(psth(i, bins>0.002) < (1-std(psth(i,bins<-0.002))));
%     binaryresp(i,:) = zeros(1, length(psth));
%     binaryresp(i,idx(binpos)) = ones(1,length(binpos));
%     binaryresp(i,idx(binneg)) = ones(1,length(binneg))*-1;
%     %use the bin response to classify the cluster as inhibitory/excitatory
%     %response to the stimulus: define which is the first response lasting
%     %for at least 20ms (decrease in firing VS increase in firing)
%     idxdisc = find(diff(binaryresp(i,:)~=0)); %find discontinuity points
%     nbins = [];
%     if length(idxdisc)==1
%         nbins(1) = length(binaryresp)-idxdisc(1);
%     else
%         for k = 2:length(idxdisc)
%             nbins(k-1) = idxdisc(k)-idxdisc(k-1);
%         end
%     end
%     idxmin = find(nbins>=20);
%     excludezeros = binaryresp(i,idxdisc(idxmin(:))+1)~=0; 
%     if ~isempty(excludezeros)
%         signflag(i) =  binaryresp(i,idxdisc(idxmin(find(excludezeros,1,'first')))+1);
%     else
%         signflag(i) = 0;
%     end
end

%check cross-correlation of all psth couples
%combs = nchoosek(1:length(clusters),2);
%C = corrcoef(psth');
    
    