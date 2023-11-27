function [C] = channel_crosscorr(varargin)

%check cross-correlation between channels in LFP recordings.
% Inputs:
%   - data (required) = lfp matrix (chs x time)
%   - select_channels (optional) = specify to use only a subset of channels;
%       default all channels are selected
%   - filter_signal (optional) = specify a frequency band to band-pass filter
%       the signal; default [0 100] Hz
%   - select_interval (optional) = use only a part of the recording; default the entire recording is used 
%   - fs (optional) = sampling frequency; default 30000
%
%   Note: filtering the signal for the entire recording for all channels is
%   time consuming. Considering prefiltering and save the filtered signal,
%   or use only smaller portions of the recording. To avoid filtering
%   within this function, specify an empty vector [] for filter option
%
%   Mattia D'Andola, Aug 2020

def_chs = [];
def_filt = [1 100];
def_fs = 30000;
def_time = [];

%inputs control
parser = inputParser;
addRequired(parser, 'data', @ismatrix)
addParameter(parser, 'select_channels',def_chs);
addParameter(parser, 'filter_signal', def_filt);
addParameter(parser, 'select_interval', def_time);
addParameter(parser, 'fs', def_fs);

parse(parser,varargin{:})
Data = parser.Results.data;
sel_chs = parser.Results.select_channels;
filt_sign = parser.Results.filter_signal;
sel_time = parser.Results.select_interval;
fs = parser.Results.fs;

if isempty(sel_chs)
    sel_chs = 1:size(Data,1);
end

%cut data interval
if ~isempty(sel_time)
    Data = Data(:, sel_time(1)*fs:sel_time(2)*fs);
end

%filter signal
if ~isempty(filt_sign)
    Data = data_filter(Data, 'bandpass', filt_sign, 3, 'fs', 30000, 'chs', sel_chs);
end

%select channels
% anmtx = double(Data(sel_chs,:));
anmtx = Data(sel_chs,:);

% [C, Clags] = xcorr(anmtx','coeff');

C = corrcoef(double(anmtx'));
figure
imagesc(C)
colormap(jet)
colorbar

if isnan(anmtx)
    [~, ancol] = find(isnan(anmtx)==1);
    anmtx = anmtx(:,1:ancol(1)-1);
    [~, aicol] = find(isinf(anmtx)==1);
    anmtx = anmtx(:,1:aicol(1)-3);
end
combs = nchoosek(1:size(anmtx,1),2);
for i = 1:size(combs,1)
%     [Cxy(i,:),F] = mscohere(anmtx(combs(i,1),:),anmtx(combs(i,2),:),hamming(100*100),80*100,100*100,fs);
    % [Cxy,F] = mscohere(anmtx(5,:)',anmtx(6,:)');
    [Cxy(i,:),F] = mscohere(anmtx(combs(i,1),:),anmtx(combs(i,2),:));
    
    [Pxy(i,:),F] = cpsd(double(anmtx(combs(i,1),:)),double(anmtx(combs(i,2),:)));
%     [Pxy(i,:),F] = cpsd(double(anmtx(combs(i,1),:)),double(anmtx(combs(i,2),:)),hamming(100*100),80*100,100*100,fs);
    % [Pxy,F] = cpsd(double(anmtx(5,:)),double(anmtx(6,:)));
    Pxy(Cxy < 0.2) = 0;
end

plot(F,angle(Pxy(120,:)/pi))
title('Cross Spectrum Phase')
xlabel('Frequency (Hz)')
ylabel('Lag (\times\pi rad)')
grid

