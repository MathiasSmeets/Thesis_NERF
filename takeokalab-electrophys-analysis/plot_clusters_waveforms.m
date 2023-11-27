function plot_clusters_waveforms(varargin)

% def_chsel = [];
def_nch = 384;
def_Fs = 30000;

%inputs control
parser = inputParser;
addRequired(parser, 'sorted_path', @ischar);
addRequired(parser, 'cluster', @isnumeric);
% addParameter(parser, 'chsel', def_chsel);
% addParameter(parser, 'principal_channel', def_chsel_scat);
addParameter(parser, 'nch', def_nch);
addParameter(parser, 'Fs', def_Fs);

parse(parser,varargin{:})
dataPath = parser.Results.sorted_path;
clusel = parser.Results.cluster;
% chsel = parser.Results.chsel;
% chsel_scat = parser.Results.principal_channel;
nch = parser.Results.nch;
Fs = parser.Results.Fs;

sp = loadKSdir(dataPath);
[data, header, raw] = tsvread(fullfile(dataPath,'cluster_info.tsv'));
idx = find(data(:,1)==clusel);
chsel_scat = data(idx,6);
chsel = [chsel_scat-5:chsel_scat+5];
idxrange = chsel >=0 & chsel<=nch;
chsel = chsel(idxrange);

gwfparams.dataDir = dataPath;
apD = dir(fullfile(dataPath, '*.dat')); 
gwfparams.fileName = apD(1).name;     
gwfparams.dataType = 'int16';
gwfparams.nCh = nch;
gwfparams.wfWin = [-40 41];
gwfparams.nWf = 2000;
gwfparams.spikeTimes = ceil(sp.st(sp.clu==clusel)*Fs);
gwfparams.spikeClusters = sp.clu(sp.clu==clusel);

wf = getWaveForms(gwfparams);

wfc = squeeze(wf.waveFormsMean);
if isempty(chsel)
    chsel = 1:size(wfc,1);
end

figure
imagesc(wfc(chsel,:));
set(gca, 'YDir', 'normal','Ytick',5:5:chsel(end), 'Yticklabel',chsel(5:5:end)); 
xt = get(gca,'xtick');
xt = xt/30;
set(gca, 'Xticklabels',num2str(xt(:),'%.2f'))  
xlabel('time (ms)'); 
ylabel('channel number');
colormap(colormap_BlueWhiteRed);
box off
title(sprintf('Unit n. %d', clusel))
ci = colorbar;
ylabel(ci, 'Amplitude (\muV)')

if isempty(chsel_scat)
    [~,idm] = max(max(pr, [], 2));
else
    idm = chsel_scat;
end
figure, plot(wfc(idm,:))
ylabel('Amplitude (\muV)')
xlabel('time (samples)')
xt = get(gca,'xtick');
xt = xt/30;
set(gca, 'Xticklabels',num2str(xt(:),'%.2f'))  
xlabel('time (ms)'); title(sprintf('Unit n. %d, channel %d',clusel,idm))
