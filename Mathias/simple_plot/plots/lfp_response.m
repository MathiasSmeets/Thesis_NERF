function lfp_response(data, evs, interval, chs, fs)

% Plots lfp response to the stimulation events contained in the structure
% evs (which is the output of open-ephys recording preprocessing - see
% function preprocess_openephys_recording.m). 
%   Inputs:
%   - data: data matrix [nchs x nsamples]
%   - evs = events structure
%   - interval = seconds before and after stim
%   - chs specifies channels to be plotted
% plotted
%
% Mattia D'Andola,
% March 2020

if nargin < 4 || isempty(chs)
    chs = 1:size(data,1);
end

if nargin < 5
    fs = 30000;
end

[rc, n] = numSubplots(length(chs));
figure
lfptemp = [];
nt = interval*fs;
for i = 1:length(chs)
    for j = 1:size(evs.onsets,1)
%         disp(sprintf('j=%d, index = [%d %d]',j,  evs.onsets(j)*fs - nt(1),evs.onsets(j)*fs + nt(2)))

        if (evs.onsets(j)*fs - nt(1))>1 && (evs.onsets(j)*fs + nt(2))<=size(data,2)
            lfptemp(j,:) = data(chs(i), evs.onsets(j)*fs - nt(1):evs.onsets(j)*fs + nt(2));
        end     
    end
    lfpmean(i,:) = mean(lfptemp);
%     lfptemp = [];
% end

% x = linspace(nt(1),nt(2),size(lfpmean,2));
% y = bsxfun(@plus, lfpmean, (0:size(lfpmean,1)-1)'*2);
% figure
% plot(x, y);
% grid on
% colormap copper
    
    subplot(rc(1), rc(2), i)
    hold on
    for j = 1:size(lfptemp,1)
        h = plot(lfptemp(j,:),'-b', 'LineWidth', 0.5); 
        h.Color(4) = 0.2;
    end
    plot(lfpmean(i,:), '-k', 'LineWidth', 2)
    title(sprintf('Channels %d',i))
    lfptemp = [];
    
end
