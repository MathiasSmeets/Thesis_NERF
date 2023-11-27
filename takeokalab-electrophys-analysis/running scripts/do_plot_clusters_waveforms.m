clusel = 1;
chsel = [];
chsel_scat = 2;
myKsDir = pwd;

sp = loadKSdir(myKsDir);
gwfparams.dataDir = myKsDir;
apD = dir(fullfile(myKsDir, '*.dat')); 
gwfparams.fileName = apD(1).name;     
gwfparams.dataType = 'int16';
gwfparams.nCh = 384;
gwfparams.wfWin = [-40 41];
gwfparams.nWf = 2000;
gwfparams.spikeTimes = ceil(sp.st(sp.clu==clusel)*30000);
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
    [m,idm] = max(max(pr, [], 2));
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
