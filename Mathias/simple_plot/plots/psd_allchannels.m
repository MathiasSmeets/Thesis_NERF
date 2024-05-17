function [allPowerEst, F, allPowerVar] = psd_allchannels(varargin)
%Initialisers for the parser
def_plt = 1;
def_probe = 'Neuropixels3A';
def_lfp = true;
def_tr = [1 100];
def_chsel = [];

lfpFs = 2500;  % neuropixels phase3a sampling rate for lfp
nChansInFile = 385; % neuropixels phase3a, from spikeGLX
%parse the input
parser = inputParser;
addRequired(parser, 'data', @ischar)
addParameter(parser, 'plot_results', def_plt);
addParameter(parser, 'probe', def_probe);
addParameter(parser, 'lfp', def_lfp);
addParameter(parser, 'time_range', def_tr); %in seconds
addParameter(parser, 'chsel', def_chsel); %in seconds

parse(parser,varargin{:})
data = parser.Results.data;
plt = parser.Results.plot_results;
probe = parser.Results.probe;
lfp = parser.Results.lfp;
tr = parser.Results.time_range;
chsel = parser.Results.chsel;
%Load data
lfpD = dir(fullfile(data, '*.lf.bin')); % LFP file from spikeGLX specifically, get data information
lfpFilename = fullfile(data, lfpD(1).name);%extract data location
%Estimate the PSD
[allPowerEst, F, allPowerVar] = ...
    power_spectrum_density(lfpFilename, 'Fs', lfpFs, 'nch', nChansInFile, 'time_range', tr, 'chsel', chsel);


