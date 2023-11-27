function analyze_spiking_phase(varargin)

% analyze spike timing phase in relation to LFP oscillations at resonance
% frequencies
%
% INPUTS
%   required:
%       1. Path to folder containing raw/preprocessed data
%       2. Path to sorting folder
%   options:
%       'recording_system': 'SpikeGLX' (default), 'Plexon', 'OpenEphys'
%       'nch': number of channels (default 385 - Neuropixels 3A)
%       'clusters': select spiking clusters (default all of them)
% Mattia D'Andola
% April 2021

def_rawpath = ''; 
def_spkpath = '';
def_recsys = 'SpikeGLX';
def_nch = 385;

%inputs control
parser = inputParser;
addRequired(parser, 'raw_path', @ischar) 
addRequired(parser, 'sorted_path', @ischar);
addParameter(parser, 'recording_system', def_recsys);
addParameter(parser, 'nch', def_nch);

parse(parser,varargin{:})
rawpath = parser.Results.raw_path;
spk_path = parser.Results.sorted_path;
recsys = parser.Results.recsys;
nch = parser.Results.nch;

%get memory pointer to raw data file
switch recsys
    case 'SpikeGLX'
        rawfile = dir(fullfile(rawpath,'*.lf.bin'));
    case {'Plexon','OpenEphys'}
        rawfile = dir(fullfile(rawpath, '*.dat'));
end
data = memmapfile(rawfile.name, 'Format', {'int16' [nch rawfile.bytes/2/nch] 'mapped'});

%upload spike file
spk = loadKSdir(pwd);
