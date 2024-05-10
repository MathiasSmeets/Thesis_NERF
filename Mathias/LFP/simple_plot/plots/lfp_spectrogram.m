function [spec,f,t] = lfp_spectrogram(data,time,chs,fs,timerange, SmoothingWindow, SubsamplingRate, nfft, FreqRange, fg)

%corticonic_plot_lfp_spectrogram(anal,opt) plot a spectrogram of all or
%selected lfp objects contained in the structure anal. Select specific
%trials or channels using the opt.isect and opt.chsel options.
%Use the opt input also to specify all the parameters for the spectrogram
%(smoothing window, subsampling rate, etc...)
%
%HARDLY RECOMENDED TO USE THIS FUNCTION FOR PATTERN LONG LESS 200s!! (THE
%FUNCTION IS REALLY MEMORY DEMANDING)
%
%Mattia D'Andola, Jun 2013
%(based on function plotLFPspectrogram.m of Lorena Perez Mendez)

if isempty(fs)
    fs = 30000;
end
if isempty(time)
    time = 1:1/fs:size(data,2)/fs;
end
if nargin<6 || isempty(SmoothingWindow)
    SmoothingWindow = 1/300;
end
if nargin<7 || isempty(SubsamplingRate)
    SubsamplingRate = 1000;
end
if nargin<8 || isempty(nfft)
    nfft = 2500;
end
if nargin<9 || isempty(FreqRange)
    FreqRange = [20 6000];
end
if nargin<10
    fg = [];
end



% if nargin<2
%     opt = options();
% elseif ~isa(opt,'options')
%     opt = options(opt);
% end
% defopt.subsample = 1;
% 
% defopt.method = 'chronox'; %'specgram'
% defopt.MuaFreqBand = [200 1500];
% defopt.SmoothingWindow = []; % s...
% defopt.SubsamplingRate = 1000;   % Hz...
% defopt.WelchWinLength = .5;
% defopt.timerange = [];
% defopt.filter = 0;
% defopt.HPFCutOffFreq = 0.3; % Hz
% 
% defopt.FreqRange = [15 90];     % Hz...
% defopt.nfft = 2500;
% 
% defopt.plot = 1;
% defopt.isect = [];
% defopt.chsel = [];
% 
% opt = setdef(opt,defopt);
% 
% if isempty(opt.SmoothingWindow)
%     opt.SmoothingWindow = 1/opt.MuaFreqBand(1);
% end
% if isempty(opt.isect)
%     opt.isect = corticonic_info_find(anal.info,'selected',1);
% end
% if isempty(opt.chsel)
%     opt.chsel = 1:nch(anal.lfp);
% end

lfpdatatot = [];
timetot = [];
for ch = chs
%         disp(sprintf('analyzing channel %d, trial %d',ch,i))
%         lfpsel = select(anal.lfp(i),ch);
    range = timerange*fs;
    if ~isempty(timerange)
        datafil = data(ch,range(1):range(2));
    end
    dt = 1/fs;
%     if opt.subsample
            % Subsample LFP...
%         lfpsel = moving_average(lfpsel,round(opt.SmoothingWindow/dt),round(1/opt.SubsamplingRate/dt));
        datafil = moving_average(datafil, time, round(SmoothingWindow/dt), round(1/SubsamplingRate/dt),1, fs);
%         datafil = datafil(ch,:);
%     end
%         % High-pass filtering of raw signal (LFP), if required...
%         dt = tsamp(lfpsel);
%         if opt.filter
%             lfpolddata = getdata(lfpsel);
%             lfpsel = moving_average(lfpsel,round(1/opt.HPFCutOffFreq/dt),round(1/opt.HPFCutOffFreq/dt/10));
%             lfpdata = getdata(lfpsel);
%             lfpdatatot = [lfpdatatot (lfpolddata.value - interp1(lfpdata.time,lfpdata.value,lfpolddata.time,'cubic'))];
%             dt = tsamp(lfpsel);
%         else
%             lfpdata = getdata(lfpsel);
%             lfpdatatot = [lfpdatatot lfpdata.value];
%         end     
    
    % compute spectrogram and plot
%     switch method
%         case 'specgram'
            window = hanning(nfft);
            [B,f,t] = spectrogram(datafil,window,fix(length(window)-1),nfft,1/dt);
            spec = 10*log10(abs(B));
%              figure;
%                 set(gcf,'Units','normalized','position',[0.1   0.3   0.8   0.5]) % [left, bottom, width, height]
%                 contourf(t+lfpdata.time(1),f,spec);
%                 shading flat
             h = fspecial('gaussian',5,3);
             specG = imfilter(spec,h);
             ind = find(f>=FreqRange(1) & f<=FreqRange(2));
%              if ~exist(fg)
%                  fg = figure;
%              end
             pcolor(t+timerange(1),f(ind),specG(ind,:));
             shading interp
                
                cb = colorbar;
                set(gca,'YLim',FreqRange);
                xlabel('Time (s)')
                ylabel('\omega/2\pi (Hz)')
                cb.Label.Interpreter = 'latex';
                cb.Label.String = 'Power $$ [10 \cdot log(mV^2/Hz)] $$';
                
%         case 'chronox'
%              disp('for chronox is necessary a subsampling. did you do it?')
%              opt.nfft = round(opt.WelchWinLength * (1/dt));
%              movingwin =  [opt.WelchWinLength 0.1];
%              tapers = [opt.WelchWinLength*10 5];
%              err=[2 0.05];
%              params =  struct('tapers', tapers, 'Fs', opt.SubsamplingRate,'err',err);  
%              [S,t,f,Serr] = mtspecgramc(lfpdatatot',movingwin,params);
%              if opt.plot
%                 figure
%                 plot_matrix(S,t+lfpdata.time(1),f,'l') ;
%                 set(gcf,'Units','normalized','position',[0.1   0.3   0.8   0.5])
%              end
%             spec = 10*log10(abs(S));
%     end
end

