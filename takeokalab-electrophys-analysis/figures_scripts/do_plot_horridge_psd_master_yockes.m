clear all

recordings_database;

ycoords = repelem(20:20:3840,2); %Coordinate in um along the neuropix probe
maxdepth = 1200; %Default, not the right depth
depth_binning = 50;
evt_name = 'horridge';

% Default depth settings
yplot = ycoords(ycoords<=maxdepth); %Default, can be conditionned to determined depth
maxch = find(ycoords<=maxdepth,1,'last');
ybin = round(0:depth_binning:maxdepth);


%%
%find mean PSD along depth during learning and after learning
for i = 1:length(raw_masters)
    disp(sprintf('master %d',i))
    load(fullfile(raw_masters{i},'events.mat'));
    evnames = {events.type};
    evtid = find(contains(evnames,evt_name));
%     trange1 = [events(evtid).onsets(1) events(evtid).offsets(end)];
    tbin = [0 60; 60 events(evtid).offsets(end); events(evtid).offsets(end) 600];
    
    [psdM1(i,:,:), FM1(i,:), ~] = psd_allchannels(raw_masters{i},'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(1,:));
    
    [psdM2(i,:,:), FM2(i,:), ~] = psd_allchannels(raw_masters{i},'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(2,:));
    
    [psdM3(i,:,:), FM3(i,:), ~] = psd_allchannels(raw_masters{i},'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(3,:));
end

%find mean PSD along depth during learning and after learning
for i = 1:length(raw_yokes)
    disp(sprintf('yokes %d',i))
    load(fullfile(raw_yokes{i},'events.mat'));
    evnames = {events.type};
    evtid = find(contains(evnames,evt_name));
%     trange1 = [events(evtid).onsets(1) events(evtid).offsets(end)];
    tbin = [0 60; 60 events(evtid).offsets(end); events(evtid).offsets(end) 600];

    [psdY1(i,:,:), FY1(i,:), ~] = psd_allchannels(raw_yokes{i},'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(1,:));
    
%     trange2 = [events(evtid).offsets(end) 600];
    [psdY2(i,:,:), FY2(i,:), ~] = psd_allchannels(raw_yokes{i},'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(2,:));
    
    [psdY3(i,:,:), FY3(i,:), ~] = psd_allchannels(raw_yokes{i},'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(3,:));
end

%average data across recordings
psdmeanM1 = squeeze(mean(psdM1));
psdmeanM2 = squeeze(mean(psdM2));
psdmeanM3 = squeeze(mean(psdM3));
stdM1 = squeeze(std(psdM1));
stdM2 = squeeze(std(psdM2));
stdM3 = squeeze(std(psdM3));
psdmeanY1 = squeeze(mean(psdY1));
psdmeanY2 = squeeze(mean(psdY2));
psdmeanY3 = squeeze(mean(psdY3));
stdY1 = squeeze(std(psdY1));
stdY2 = squeeze(std(psdY2));
stdY3 = squeeze(std(psdY3));

%smooth psd profiles with gaussian-weighted moving average + convert to db
psdplotM1 = smoothdata(10*log10(psdmeanM1),2,'gaussian');
psdplotM2 = smoothdata(10*log10(psdmeanM2),2,'gaussian');
psdplotM3 = smoothdata(10*log10(psdmeanM3),2,'gaussian');
psdplotY1 = smoothdata(10*log10(psdmeanY1),2,'gaussian','SmoothingFactor',.15); %different level of noise
psdplotY2 = smoothdata(10*log10(psdmeanY2),2,'gaussian','SmoothingFactor',.15); %different level of noise
psdplotY3 = smoothdata(10*log10(psdmeanY3),2,'gaussian','SmoothingFactor',.15); %different level of noise

%plot results for all channels
figure
F = FM1(1,:);
subplot(2,3,1)
imagesc(F,yplot,psdplotM1)
caxis([-5 20])
colormap('jet')
xlabel('F [Hz]')
ylabel('depth on probe (µm)')
set(gca, 'ydir', 'normal', 'xlim', [0 100])
title('Masters during adaptation')

subplot(2,3,2)
imagesc(F,yplot,psdplotM2)
caxis([-5 20])
colormap('jet')
xlabel('F [Hz]')
ylabel('depth on probe (µm)')
set(gca, 'ydir', 'normal', 'xlim', [0 100])
title('Masters during consolidation')

subplot(2,3,3)
imagesc(F,yplot,psdplotM3)
caxis([-5 20])
colormap('jet')
xlabel('F [Hz]')
ylabel('depth on probe (µm)')
set(gca, 'ydir', 'normal', 'xlim', [0 100])
title('Masters during retaining')

subplot(2,3,4)
imagesc(F,yplot,psdplotY1)
caxis([-5 20])
colormap('jet')
xlabel('F [Hz]')
ylabel('depth on probe (µm)')
set(gca, 'ydir', 'normal', 'xlim', [0 100])
title('Yokes during adaptation')

subplot(2,3,5)
imagesc(F,yplot,psdplotY2)
caxis([-5 20])
colormap('jet')
xlabel('F [Hz]')
ylabel('depth on probe (µm)')
set(gca, 'ydir', 'normal', 'xlim', [0 100])
title('Yokes during consolidation')

subplot(2,3,6)
imagesc(F,yplot,psdplotY3)
caxis([-5 20])
colormap('jet')
xlabel('F [Hz]')
ylabel('depth on probe (µm)')
set(gca, 'ydir', 'normal', 'xlim', [0 100])
title('Yokes during retaining')

hp4 = get(subplot(2,3,6), 'Position');
cbar = colorbar('Position',  [hp4(1)+hp4(4)+0.01  hp4(3)  0.03  hp4(3)+hp4(4)*2.1]);
cbar.Label.String = 'PSD [db/Hz]';
