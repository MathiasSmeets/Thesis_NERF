function plot_raw_example(varargin)

def_trange = [10 20];
def_nch = 385;
def_pltype = [1 2];

parser = inputParser;
addRequired(parser, 'raw_path', @ischar);
addRequired(parser, 'channel', @isnumeric);
addParameter(parser, 'trange', def_trange);
addParameter(parser, 'nch', def_nch);
addParameter(parser, 'plot_type', def_pltype);

parse(parser,varargin{:})
raw_path = parser.Results.raw_path;
ch = parser.Results.channel;
trange = parser.Results.trange;
nch = parser.Results.nch;
pltype = parser.Results.plot_type;

figure;

%low freq 
Fs = 2500;
d = dir(fullfile(raw_path,'*.lf.bin*')); 
nSamps = d.bytes/2/nch;
mmf = memmapfile(fullfile(d.folder,d.name), 'Format', {'int16', [nch nSamps], 'x'});
data = double(mmf.Data.x(ch, Fs*trange(1):Fs*trange(2)));
if pltype == [1 2]
    subplot(2,1,1)
    plot(trange(1):1/Fs:trange(2), data)
    xlabel('Time [s]')
    ylabel('LFP (Fs 2500 Hz)')
elseif pltype == 1
    plot(trange(1):1/Fs:trange(2), data)
    xlabel('Time [s]')
    ylabel('LFP (Fs 2500 Hz)')
end
%high freq
Fs = 30000;
d = dir(fullfile(raw_path,'*.ap.bin*')); 
nSamps = d.bytes/2/nch;
mmf = memmapfile(fullfile(d.folder,d.name), 'Format', {'int16', [nch nSamps], 'x'});
data = double(mmf.Data.x(ch, Fs*trange(1):Fs*trange(2)));
if pltype == [1 2]
    subplot(2,1,2)
    plot(trange(1):1/Fs:trange(2), data)
    xlabel('Time [s]')
    ylabel('LFP (Fs 30 kHz, filtered highpass 300 Hz)')
elseif pltype == 2
    plot(trange(1):1/Fs:trange(2), data)
    xlabel('Time [s]')
    ylabel('LFP (Fs 30 kHz, filtered highpass 300 Hz)')
end




