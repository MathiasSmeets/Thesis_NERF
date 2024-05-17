%Same as power_spectrum_density but now you filter your results to not have
%any DC component (used for the localPSD and gifPSD)
function [allPowerEst, F, allPowerVar] = power_spectrum_density_filter(varargin)

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

%% SmoothData
%
%  data_to_smooth=smoothdata(mmf.Data.x,'SmoothingFactor',.1);
%     
%     [a b]=findpeaks(-data_to_smooth,'MinPeakProminence',20);
% 
%     for i=1:numel(a)
%         
%         data_smooth(1,b(i)-20:b(i)+20)= data_smooth(1,b(i)-20);
%         
%     end
%         
%%   
   
    sampStarts = tr(1)*Fs;%round(linspace(Fs, nSamps, nClips+1)); 
    nClipSamps = Fs*clipDur; %round(Fs*clipDur);
    if isempty(chsel)
        chsel = 1:nch;
    end
    
    % pull out the data
    data = double(mmf.Data.x(chsel, (1:nClipSamps)+sampStarts));
        
    % mean (or median?) subtract 
%     data = bsxfun(@minus, data, median(thisDat));
    data = bsxfun(@minus, data, mean(data,2));
    
    
    
%%     
% data_to_smooth=double(data(1,:));
% data_smooth=data;
% [a b]=findpeaks(-data_to_smooth,'MinPeakProminence',20);
% 
% for i=2:numel(a)
%     for j=1:chsel(end)
%         
%     data_smooth(j,b(i)-20:b(i)+20)= data_smooth(j,b(i)-20);
%     
%     end
% end
    data_smooth=data;
    
%%
    %Fs = 1/(100/size(data,2))
    %t = 0:1/Fs:100;
    %data = sin(t*10*2*pi);
    %compute PSD with Welch method (all frequencies up to Nyquist freq)
    butter_bp_filter = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',30,'SampleRate',Fs,'DesignMethod','butter');
        
    filtered_data = (filtfilt(butter_bp_filter,double(data')))';%filter
    



    [F, Pxx] = mpsd(filtered_data, Fs, 1:size(filtered_data,1), 'welch', 0, Fs);
    F = F';
    
    allPowerEst = zeros(1, size(Pxx,1), size(Pxx,2));%if bins are separated, make for loop to loop through different bins and calculate multiple PSDs (for each of the bins), then initialise the output here if n == 1 and USE N IN THE FOLLOWING LINE
    
    allPowerEst(1,:,:) = Pxx;     
    
    
    allPowerVar = squeeze(var(allPowerEst,1));
    allPowerEst = squeeze(mean(allPowerEst, 1));
end


% function [Pxx, F] = myTimePowerSpectrumMat(x, Fs)
%     L = size(x,1);
%     NFFT = 2^nextpow2(L);
%     [Pxx,F] = pwelch(x,[],[],NFFT,Fs);
% end