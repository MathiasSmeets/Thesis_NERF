clear;clc;close all;
addpath(genpath("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Thesis_NERF\Mathias"))
myKsDir = 'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220725\20220725_long_long_recording_g0';

ycoords = repelem(20:20:3840,2);
maxdepth = 1500;
yplot = ycoords(ycoords<=maxdepth);
maxch = find(ycoords<=maxdepth,1,'last');
freqBands = {[0.5 4], [4 8], [8 12], [12 35], [35 200]};
binSize = 10;%how many seconds you want each bin to be
endrest = 7443;
starttime = 4438;
endtime = 7443;
counter = 1;
for i=floor(starttime/binSize):floor(endtime/binSize) % to cover 0 to 2400 seconds
    %for i=225*25/binSize:275*25/binSize % to cover 0 to 2400 seconds
    disp(floor(starttime/binSize) + "--" + i + "--" + floor(endtime/binSize))
    % Creating the tbin of 10 seconds
    tbin = [binSize*(i-1) binSize*i];


    % Exctracting the LFP
    % psdM are the LFP along ~4000Hz frequency band
    % FM is just the number of the frequency
    % allPowerVar is the variance for each frequecny - we don't use it
    [psdM, FM, allPowerVar] = psd_allchannels(myKsDir,'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(1,:));

    % Look at a specific freq like 1.5 to 4Hz here
    Freq=freqBands{1};
    psd_freq_1(counter,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

    Freq=freqBands{2};
    psd_freq_2(counter,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

    Freq=freqBands{3};
    psd_freq_3(counter,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

    Freq=freqBands{4};
    psd_freq_4(counter,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

    Freq=freqBands{5};
    psd_freq_5(counter,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

    %smooth psd profiles with gaussian-weighted moving average + convert to db
    psdplotM = smoothdata(10*log10(psdM),2,'gaussian','SmoothingFactor',.01);
    counter = counter + 1;

    % - - - - - - - -plot results for all channels
    % figure(1)
    % F = FM(1,:);
    % subplot(10,Column(i),i)
    % imagesc(0:maxtime/size(psd_freq_1,1):maxtime-1,1:size(psd_freq_1,2)*10,F,-yplot,psdplotM)
    % %caxis([-100 100])
    % colormap('jet')
    % %xlabel('F [Hz]')
    % %ylabel('depth on probe (ï¿½m)')
    % set(gca, 'ydir', 'normal', 'xlim', [0 100])
    % %title('Masters during adaptation')
    %
    % cbar.Label.String = 'PSD [db/Hz]';

end