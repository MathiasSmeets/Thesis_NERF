function [allPowerEst, F, allPowerVar] = power_spectrum_density(varargin)

    def_nch = 385;
    def_fs = 1000;
    def_tr = [];
    def_clipdur = 10;
    def_freqBand = {[0 100]};
    def_chsel = [];
    
    %inputs control
    parser = inputParser;
    addRequired(parser, 'data_path', @ischar);
    addParameter(parser, 'nch', def_nch);
    addParameter(parser, 'Fs', def_fs);
    addParameter(parser, 'time_range', def_tr); %in seconds
    addParameter(parser, 'clips_duration', def_clipdur); %in seconds
    addParameter(parser, 'frequency_band', def_freqBand);
    addParameter(parser, 'chsel', def_chsel);
    parse(parser,varargin{:})
    data_path = parser.Results.data_path;
    nch = parser.Results.nch;
    Fs = parser.Results.Fs;
    tr = parser.Results.time_range;
    clipDur = parser.Results.clips_duration;
    freqBand = parser.Results.frequency_band;
    chsel = parser.Results.chsel;
        
    %get pointer to data in memory
    d = dir(data_path); 
    nSamps = d.bytes/2/nch;
    mmf = memmapfile(data_path, 'Format', {'int16', [nch nSamps], 'x'});

    % split time range into 10s long clips to reduce memory demand
    tClips = tr(1):10:tr(2);
    nClips = length(tClips)-1;
    sampStarts = tClips*Fs;%round(linspace(Fs, nSamps, nClips+1)); 
    nClipSamps = round(Fs*clipDur);
    if isempty(chsel)
        chsel = 1:nch;
    end
    
    for n = 1:nClips
        fprintf(1, 'clip %d of %d\n', n, nClips); 
        % pull out the data
        data = double(mmf.Data.x(chsel, (1:nClipSamps)+sampStarts(n)));
    
        % mean (or median?) subtract 
    %     data = bsxfun(@minus, data, median(thisDat));
        data = bsxfun(@minus, data, mean(data,2));
        
        %compute PSD with Welch method (all frequencies up to Nyquist freq)
        [F, Pxx] = mpsd(data, Fs, 1:size(data,1), 'welch', 0, Fs);
        F = F';
        
        if n==1
            allPowerEst = zeros(nClips, size(Pxx,1), size(Pxx,2));
        end
        allPowerEst(n,:,:) = Pxx;     
    end
    
    allPowerVar = squeeze(var(allPowerEst,1));
    allPowerEst = squeeze(mean(allPowerEst, 1));
end


% function [Pxx, F] = myTimePowerSpectrumMat(x, Fs)
%     L = size(x,1);
%     NFFT = 2^nextpow2(L);
%     [Pxx,F] = pwelch(x,[],[],NFFT,Fs);
% end