function [xscale, yscale] = waveform_examples_and_correlograms(dataPath, clusters, chsel)

% dataPath: path to the spike sorted data folder
% clusters: cluster id of the clusters to plot
% chsel = channels to plot (must make sense given the selected clusters)


Fs = 30000;
nWf = 500;
dataType = 'int16';
nch = 384;
wfWin = [-50 51];
xscale = wfWin/Fs;
bin_interval = 0.0005; %bin size for correlograms (in seconds)
max_lag = 0.05; %max lag for correlograms, in seconds

load(fullfile(dataPath,'events.mat'))


%create data pointer to read the .dat file in memory
apD = dir(fullfile(dataPath, '*.dat')); 
nSamp = apD(1).bytes/2/nch;
mmf = memmapfile(fullfile(apD(1).folder,apD(1).name), 'Format', {dataType, [nch nSamp], 'x'});

cmap = parula;
colors = cmap(1:floor(size(cmap,1)/length(clusters)):end,:);

% load data and identify principal channel of each cluster
sp = loadKSdir(dataPath);
[data, header, raw] = tsvread(fullfile(dataPath,'cluster_info.tsv'));

% for i = 1:length(clusters)
%     idx = find(data(:,1)==clusters(i));
%     chsel(i) = data(idx,6);
%     chsel = [chsel_scat-floor(n_channels/2):chsel_scat+ceil(n_channels/2)];
%     idxrange = chsel >=0 & chsel<=nch;
%     chsel = chsel(idxrange);


for i = 1:length(clusters)
%     disp(sprintf('cluster %d',i))
    idx = find(data(:,1)==clusters(i));
    chsel_scat = data(idx,6);
    % get spiketimes
    spikeTimes = ceil(sp.st(sp.clu==clusters(i))*Fs);
    if ~isempty(chsel_scat) && chsel_scat~=0 && numel(spikeTimes)>=nWf
        spikeTimesRP = spikeTimes(randperm(size(spikeTimes,1)));
        spikeTimeKeeps = sort(spikeTimesRP(1:nWf));
        %upload wavefors at spike times
        for s = 1:nWf
            tmpWf = mmf.Data.x(chsel,spikeTimeKeeps(s)+wfWin(1):spikeTimeKeeps(s)+wfWin(end));
            waveForm(1:length(chsel),:,s) = tmpWf;
        end
        %average
        waveFormsMean(1:length(chsel),:,i) = nanmean(waveForm,3);
        disp(['Completed unit ' int2str(i) ' of ' int2str(length(clusters)) '.']);
    else
        disp(['Cluster ', int2str(i) ' has less than the required ' int2str(nWf) ' waforms'])
        waveFormsMean(i,length(chsel),:,i) = NaN(1,numel(wfWin(1):wfWin(2)));
    end
%     [~,idm] = min(waveFormsMean(i,:));
%     [~,idM] = max(waveFormsMean(i,:));
%     
%     L(i) = abs(idM - idm)/Fs; %wrong, use the calculation in do_group_neural_response
%     L(i) = abs(idM - idm)/Fs; %wrong, use the calculation in do_group_neural_response
end

yscale = [min(min(min(waveFormsMean))), max(max(max(waveFormsMean)))];
% figure 1: example waveforms
figure
%cont = 0;
tt = tiledlayout(length(chsel),length(clusters));
for i = 1:length(chsel)
    for j = 1:length(clusters)
        %cont = cont+1;
        %subplot(length(chsel),length(clusters), cont)
        nexttile
        plot(waveFormsMean(i,:,j), 'LineWidth', 2, 'color', colors(j,:))
        ylim(yscale)
        set(gca, 'visible', 'off')
    end
end
tt.TileSpacing = 'none';
tt.Padding = 'tight';

%compute  and plot correlograms
figure
cont = 0;
M = zeros(length(clusters),length(clusters));
for i = 1:length(clusters)
    for j = 1:length(clusters)
        cont = cont+1;
        M(i,j) = cont;
    end
end
for i = 1:length(clusters)
    for j = i:length(clusters)
        subplot(length(clusters),length(clusters),M(i,j))

        st = sp.st((sp.clu==clusters(i) | sp.clu==clusters(j)) & sp.st<events.onsets(1));
        if i == j
            sv = ones(1, length(st));
            [ccg, bins] = correlogram(st, sv, bin_interval, max_lag);
            col = colors(i,:);
            bar(bins,ccg(:,1,1),'EdgeColor', col, 'FaceColor', col)
        else
            sv = sp.clu(sp.clu==clusters(i) | sp.clu==clusters(j));
            sv = sv(1:length(st));
            sv(sv==clusters(i)) = 1;
            sv(sv==clusters(j)) = 2;
            [ccg, bins] = correlogram(st, sv, bin_interval, max_lag);
            col = 'k';
            bar(bins,ccg(:,2,1),'EdgeColor', col, 'FaceColor', col)
        end
    end
end
end


