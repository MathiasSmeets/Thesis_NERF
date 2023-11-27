function [data, events] = preprocess_plexon_recording(varargin)

% Upload plexon files and save infos and data separately
%
% N.B The function needs the readPLXFile MEX file, please compily before
% usage by running the script 'build_readPLXFileC';
%
% INPUT NOTES: 
%   - filepath = complethe path to the .plx file. If you recorded into
%                 .pl2 format, you can use the 'PlX Utilities' app to convert it into .plx
%   - event_channels: this must be a matrix [n x 4], with n = number of
%   different type of events. In Plexon recordings, every event is
%   represented in 4 different channels for: start (absolute time), onsets,  stop
%   (absolute time), offsets. Usually those 
%
% Mattia D'Andola, 2019
def_continuous = 1; 
def_spikes = 0;
def_waves = 0;
def_evt = [];
def_range = [];
def_ach = [];
def_export = 1;
def_evtch = [17 18 19 20];
def_evtnames = {};

%inputs control
parser = inputParser;
addRequired(parser, 'filepath', @ischar)
addParameter(parser, 'continuous', def_continuous);
addParameter(parser, 'spikes', def_spikes);
addParameter(parser, 'waves', def_waves);
addParameter(parser, 'events', def_evt);
addParameter(parser, 'range', def_range);
addParameter(parser, 'active_channels', def_ach);
addParameter(parser, 'export_binary', def_export);
addParameter(parser, 'events_channels', def_evtch);
addParameter(parser, 'events_names', def_evtnames);

parse(parser,varargin{:})
filepath = parser.Results.filepath;
cont_flag = parser.Results.continuous;
spk_flag = parser.Results.spikes;
wave_flag = parser.Results.waves;
evt_flag = parser.Results.events;
trange = parser.Results.range;
ach = parser.Results.active_channels;
export = parser.Results.export_binary;
evtch = parser.Results.events_channels;
evtnames = parser.Results.events_names;

%use MEX functions to upload the data info
if cont_flag 
    cont_in = 'continuous';
else
    cont_in = 'nocontinuous';
end
if spk_flag
    spk_in = 'spikes';
else
    spk_in = 'nospikes';
end
if wave_flag
    wv_in = 'waves';
else
    wv_in = 'nowaves';
end
if ~isempty(evt_flag)
    evt_in = 'events';
else
    evt_in = 'noevents';
end

if isempty(trange)
    plx = readPLXFileC(filepath,cont_in,spk_in,wv_in,evt_in);
else
    plx = readPLXFileC(filepath,cont_in,spk_in,wv_in,evt_in,'range',trange);
end

data.source = filepath;
data.Fs = plx.ADFrequency;
data.timerange = trange;
if cont_flag
    for k = 1:length(plx.ContinuousChannels)
        names{k} = plx.ContinuousChannels(k).Name;
    end
    for i = 1:length(ach)
        iach = find(contains(names,['WB',num2str(ach(i),'%02.f')]));
        mdata(i,:) = plx.ContinuousChannels(iach).Values;
        mnames{i} = names(iach);
    end
    data.continuous.values = mdata;
    if export
       datanew = reshape(mdata,[1, size(mdata,1)*size(mdata,2)]);
        disp('Writing new data file for spike sorting...')
        fid = fopen('recording.dat', 'w');
        fwrite(fid, datanew, 'int16');
        fclose(fid); 
    end
    data.continuous.names = mnames;
    ts = 1/data.Fs;
    data.continuous.time = 1*ts:ts:size(data.continuous.values,2)*ts;
end
if spk_flag
    data.spikes(:) = plx.SpikeChannels;
end
if ~isempty(evt_flag)
    for i = 1:size(evtch,1)
        events(i).type = evtnames{i};
        events(i).start = double(plx.EventChannels(evtch(i,1)).Timestamps)/data.Fs;
        events(i).stop = double(plx.EventChannels(evtch(i,3)).Timestamps)/data.Fs;
        events(i).onsets = double(plx.EventChannels(evtch(i,2)).Timestamps)/data.Fs;
        events(i).offsets = double(plx.EventChannels(evtch(i,4)).Timestamps)/data.Fs;
    end
end

    