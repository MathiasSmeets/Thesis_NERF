function data = data_filter(varargin)

% Filters data using parameters reported in the structure params, to be
% created as follow:
%   - type of filter: params.type = 'notch'; ('bandpass', 'lowpass','highpass')
%   - cutting frequency: params.fc = [50]; (for notch, lowpass or highpass), or .fc = [f1 f2]
%                        for bandpass
%   - filter width: params.wd = 35 (for notch, measure in db)
%   - sampling frequency: params.fs = 30000;
%   - channels to be filtered: params.chs = [i, j, ...]; if empty, all
%   channels are filtered
% Mattia D'Andola, May 2020



%defineinputs
def_fs = 30000;
def_chs = [];
parser = inputParser;
addRequired(parser, 'data', @ismatrix)
addRequired(parser, 'type', @ischar)
addRequired(parser, 'fc', @isnumeric)
addRequired(parser, 'wd', @isnumeric)
addParameter(parser, 'fs',def_fs);
addParameter(parser, 'chs', def_chs);
parse(parser,varargin{:})

data = parser.Results.data;
type = parser.Results.type;
fc = parser.Results.fc;
wd = parser.Results.wd;
chs = parser.Results.chs;
fs = parser.Results.fs;

if isempty(chs)
    chs = 1:size(data,1);
end


%create filter function
switch type
    case 'notch'
        w0 = fc/(fs/2);
        bw = w0/wd(1);
        [b,a] = iirnotch(w0,bw,-wd(2));
        for i = chs
            disp(sprintf('Filtering ch %d',i))
            dtf = data(i,:);
            data(i,:) = filter(b,a,dtf);
        end
        
    case 'harmonics'
        fcuts = [];
        mags = 1;
        devs = 0.05;
        for i = 1:wd(2)
            fcuts = [fcuts (fc*i)-wd(1)*2 (fc*i)-wd(1) (fc*i)+wd(1) (fc*i)+wd(1)*2]; %Frequencies
%             if i == 1
%                 mags =  [mags 0];  % Main Stopband
%                 devs = [devs 0.05];
%             else
                mags = [mags 0 1]; % following passband + harmonic stopband
                devs = [devs 0.01 0.05];
%             end
        end
%         fcuts = [fcuts fs/2-1];
        [n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fs);          % Kaiser Window FIR Specification
        n = n + rem(n,2);
        hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');           % Filter Realisation
        figure
        freqz(hh,1,fs,fs)
        set(subplot(2,1,1), 'XLim',[0 100]);                        % Zoom Frequency Axis
        set(subplot(2,1,2), 'XLim',[0 100]);                        % Zoom Frequency Axis
        title('filter representation')
        
        for i = chs
            disp(sprintf('Filtering ch %d',i))
            dtf = data(i,:);
            data(i,:) = filter(hh,1,dtf);
        end
    case 'bandpass'
        [z,p,k] = butter(wd,fc/(fs/2),'bandpass');
        [sos,g]=zp2sos(z,p,k);
        for i = chs
            disp(sprintf('Filtering ch %d',i))
            dtf = data(i,:);
            data(i,:) = filtfilt(sos,g,dtf);
        end
end




