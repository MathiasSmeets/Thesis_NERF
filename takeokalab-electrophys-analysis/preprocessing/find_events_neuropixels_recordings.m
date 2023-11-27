function events = find_events_neuropixels_recordings(varargin)

def_pinbits = 1;
def_names = {'opto'};
def_fs = 30000;
def_filtlow = 0;
def_tc = 10;

%inputs control
parser = inputParser;
addRequired(parser, 'filepath', @ischar)
addParameter(parser, 'pin_bits',def_pinbits);
addParameter(parser, 'Fs', def_fs);
addParameter(parser, 'pin_names', def_names);
addParameter(parser, 'filter_lowpass', def_filtlow);
addParameter(parser, 'Tc', def_tc);

parse(parser,varargin{:})
filepath = parser.Results.filepath;
names = parser.Results.pin_names;
pinbits = parser.Results.pin_bits;
Fs = parser.Results.Fs;
filtlow = parser.Results.filter_lowpass;
Tc = parser.Results.Tc;


d = dir(filepath);
nsamps = d.bytes/2/385;
mm = memmapfile(filepath, 'Format', {'int16', [385 nsamps], 'x'});
stim = mm.Data.x(385,:);
for i = 1:length(stim)
    stimbins(:,i) = bitget(stim(i),pinbits,'int16')';
end

for i = 1:size(stimbins,1)
    if filtlow
        stimbins(:,i) = lowpass(stimbins, Fc, Fs);
    end
    events(i).type = names{i};
    %find transitions
    trpoints = diff(stimbins(i,:));
    [~, onsets] = find(trpoints == 1);
    [~, offsets] = find(trpoints == -1);

    % transform to seconds
    onsets = (onsets/Fs)';
    offsets = (offsets/Fs)';
    events(i).onsets = onsets;
    events(i).offsets = offsets;   
end

   