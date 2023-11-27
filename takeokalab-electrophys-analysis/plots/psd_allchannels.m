function [allPowerEst, F, allPowerVar] = psd_allchannels(varargin)

def_plt = 1;
def_probe = 'Neuropixels3A';
def_interval = [0 100]; 
def_tr = [1 100];
def_chsel = [];

parser = inputParser;
addRequired(parser, 'data', @ischar)
addParameter(parser, 'plot_results', def_plt);
addParameter(parser, 'probe', def_probe);
addParameter(parser, 'interval', def_interval);
addParameter(parser, 'time_range', def_tr); %in seconds
addParameter(parser, 'chsel', def_chsel); %in seconds

parse(parser,varargin{:})
data = parser.Results.data;
plt = parser.Results.plot_results;
probe = parser.Results.probe;
interval = parser.Results.interval;
tr = parser.Results.time_range;
chsel = parser.Results.chsel;

if interval(2)<=150
    lfpD = dir(fullfile(data, '*.lf.bin')); % LFP file from spikeGLX specifically
    fileflag = 1;
else
    lfpD = dir(fullfile(data, '*.ap.bin')); % LFP file from spikeGLX specifically
    fileflag = 2;
end
lfpFilename = fullfile(data, lfpD(1).name);

switch probe
    case 'Neuropixels3A'
        if fileflag == 1
            lfpFs = 2500;  % neuropixels phase3a
        elseif fileflag == 2
            lfpFs = 30000;
        end
        nChansInFile = 385; % neuropixels phase3a, from spikeGLX
end

[allPowerEst, F, allPowerVar] = ...
    power_spectrum_density(lfpFilename, 'Fs', lfpFs, 'nch', nChansInFile, 'time_range', tr, 'chsel', chsel);

% chanMap = readNPY(fullfile(data, 'channel_map.npy'));
% nC = length(chanMap);
% 
% allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq
% 
% % plot LFP power
% marginalChans = [10:50:nC];
% freqBands = {[0.5 4], [4 8], [8 13], [13 30], [30 60], [60 100]};

% plotLFPpower(F, allPowerEst, interval, marginalChans, freqBands);
