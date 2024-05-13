%
%
% - - - - LFP ANALYSIS - - - - -
%
% Simon LAVAUD & Arnaud HALLEMANS

%close all
clear;clc;close all;
addpath(genpath("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Thesis_NERF\Mathias\LFP\simple_plot"));


if ~exist('init','var')
    clc
    disp("start")
    init = true;
    
    % - - - - Initialization of parameters about the NEUROPIXEL probe
    ycoords = repelem(20:20:3840,2);
    maxdepth = 1500;
    depth_binning = 50;
    evt_name = 'horridge';
    yplot = ycoords(ycoords<=maxdepth);
    
    maxch = find(ycoords<=maxdepth,1,'last');
    ybin = round(0:depth_binning:maxdepth);
    
    lfpFs = 2500;  % neuropixels phase3a
    nChansInFile = 385;  % neuropixels phase3a, from spikeGLX
    
    % - - - - - Initialization of the frequency bands of interest
    freqBands = {[0.5 4], [4 8], [8 12], [12 35], [35 200]};
    
    % Extract the recording of interest
    % - - - - - If you want to do a for loop you can use some commands like
    % - - - - - this
    
    %RecordingNumber=28;
    %spk = loadKSdir(sort_masters_horridge{RecordingNumber});
    
    % - - - - - - - To study a unique recording
    possiblePaths={'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220725\20220725_long_long_recording_g0'; ...
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220819\20220819_long_long_recording_g0'; ... %Recording seems to be corrupted or so, only 606 seconds long
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220831\20220831_long_long_recording_g0'; ...
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220921\20220921_long_long_recording_g0'; ...
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220923\20220923_long_long_recording_g0'; ...
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20221011\20221011_long_long_recording_g0'; ...
        '\\nerffs17\takeokalabwip2023\Jeremy\Recording\test\test_g0'};
   
    possiblePathsYokes =['Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220726\20220726_long_long_recording_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220905\20220905_long_long_recording_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220906\20220906_long_long_recording_g0';%No events.mat
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220922\20220922_long_long_recording_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220927\20220927_long_long_recording_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20221006\20221006_long_long_recording_g0'];  
     possibleBetween = [
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220207\20220207_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220214\20220214_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220221\20220221_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220310\20220310_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220314\20220314_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220324\20220324_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220419\20220419_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220428\20220428_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220505\20220505_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220525\20220525_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220530\20220530_between_g0';
         'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220602\20220602_between_g0';%no events.mat
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220616\20220616_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220622\20220622_between_g0';
        
        ];
    possibleBetweenYokes = [
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220208\20220208_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220301\20220301_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220311\20220311_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220323\20220323_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220413\20220413_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220427\20220427_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220429\20220429_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220607\20220607_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220614\20220614_between_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220623\20220623_between_g0';
        ];
    timeframes = [0 15 25 40 106 145 155 170;
            0 15 25 40 106 121 131 146;
            0 15 25 40 106 121 131 146;
            0 15 25 40 107 122 133 148;
            0 15 25 40 110 132 142 157;
            0 15 25 40 110 126 136 151;
            0 15 25 40 110 126 136 151];%times represent: start horridge after rest spont horridge after end
    
    
    timeframesYokes = [0 18 28 43 109 124 134 149;
            0 15 25 40 106 121 131 146;
            0 15 25 40 106 121 131 146;
            0 15 25 40 104 129 139 154;
            0 15 25 40 161 176 186 201;
            0 15 25 40 110 127 137 152];
    
   
    %len = [size(possibleBetween,1) size(possibleBetweenYokes,1)];
    len = size(possiblePaths,1);

    % Loops to loop through all the mice
    condpos = 4;%pick mouse from yoke list or not or between condition              %%%
    for k = 1                                                                       %%%
        for j = 1:len(k)
            cond = condpos(k);
            j;
            %close all
            if cond == 1%yoke
                myKsDir=possiblePathsYokes(j, :); 
                timeframe = timeframesYokes(j,:);
            elseif cond == 2%between
                myKsDir=possibleBetween(j, :); 
                d = dir(fullfile(myKsDir, '*.lf.bin'));
                nSamps = d.bytes/2/nChansInFile;
                timeframe = [0 0 0 0 nSamps/lfpFs/60];
            elseif cond == 3%between yoke
                myKsDir=possibleBetweenYokes(j, :); 
                d = dir(fullfile(myKsDir, '*.lf.bin'));
                nSamps = d.bytes/2/nChansInFile;
                timeframe = [0 0 0 0 nSamps/lfpFs/60];
            else
                myKsDir=possiblePaths{j};
                timeframe = timeframes(j,:);
            end
            disp(myKsDir)
            %going to inspect the resting period
            timeframe(end);
            % - - - - - Load the file that correspond to the events 
            % - - - - - One event correspond to an electrical shock during the horridge
            % - - - - - paradigm
            if cond<2%If mouse recording is not split, look at the spikes as well
                load(fullfile(myKsDir,'events.mat'),'events');
                evnames = {events.type};
                evtid = find(contains(evnames,'horridge'));
                if isempty(evtid)==1
                    evtid = find(contains(evnames,'opto'));
                end
            end
            % - - - - - Isolate the pathway and extract the lfp file
            lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
            lfpFilename = fullfile(myKsDir, lfpD(1).name);

            % - - - - - - Create variable for analysis
            Column=[1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 ...
                3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4];

            %% 
            % The analysis is working as a for loop on the recording.
            % It is extracting "clip" of 10seconds (you can change it in tbin) and analyse it one by one.
            % If you want to study the first minute of recording then you have to
            % loop from 0 to 6. Then tbin will return "clip" that goes :
            % 0 to 10
            % 10 to 20
            % 20 to 30 
            % ...
            %%

            %% 
            % Keep in mind that we record NORMALLY (it can change based on the recording) like this :
            % 0 to 900 seconds : rest time - beginning of the experiment, called "the before phase"
            % 900 to 1500 seconds : horridge paradigm - where shocks and learning happen
            % 1500 to 2400 seconds : rest time - called "the after phase"


            %%
            binSize = 10;%how many seconds you want each bin to be
            endrest = timeframe(5)*60;
            endtime = timeframe(end)*60;
            for i=1:floor(endtime)/binSize % to cover 0 to 2400 seconds
            %for i=225*25/binSize:275*25/binSize % to cover 0 to 2400 seconds

                % Creating the tbin of 10 seconds
                tbin = [binSize*(i-1) binSize*i];

                % - - - - - - - - - - - - - - -
                i; %indicate in which clip you are

                % Exctracting the LFP
                % psdM are the LFP along ~4000Hz frequency band
                % FM is just the number of the frequency
                % allPowerVar is the variance for each frequecny - we don't use it
                [psdM, FM, allPowerVar] = psd_allchannels(myKsDir,'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(1,:));

                % Look at a specific freq like 1.5 to 4Hz here
                Freq=freqBands{1};
                psd_freq_1(i,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

                Freq=freqBands{2};
                psd_freq_2(i,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

                Freq=freqBands{3};
                psd_freq_3(i,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

                Freq=freqBands{4};
                psd_freq_4(i,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

                Freq=freqBands{5};
                psd_freq_5(i,:)=mean(psdM(:,find(FM>=Freq(1),1):find(FM<Freq(2),1,'last')),2);

                %smooth psd profiles with gaussian-weighted moving average + convert to db
                psdplotM = smoothdata(10*log10(psdM),2,'gaussian','SmoothingFactor',.01);


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
            %% Normalise psds
            psd_freq_1 = normrobust(psd_freq_1);
            psd_freq_2 = normrobust(psd_freq_2);
            psd_freq_3 = normrobust(psd_freq_3);
            psd_freq_4 = normrobust(psd_freq_4);
            psd_freq_5 = normrobust(psd_freq_5);

            high_bound = 3.5;%log10(mu+2*sd);
            low_bound = -3.5;%log10(mu-2*sd);
            %%
            fig = figure(1);
            fig.WindowState = 'maximized';
            sgtitle("Power in corresponding frequency bands")
            subplot(5,1,1)
            if cond < 2
                psd_freq_no_training = psd_freq_1(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
                time_to_analyse = -timeframe(3)*60+timeframe(6)*60;
            else
                time_to_analyse = timeframe(end)*60;
                psd_freq_no_training = psd_freq_1;
            end
            imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');
            % xticks(1:endrest/5/binSize:endrest/binSize)
            % xticklabels(1:endrest/60/5:endrest/60)
            caxis([low_bound high_bound])
            a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
            title("0.5-4Hz")
            xlabel("time [s]")
            ylabel("Channel depth [µm]")

            subplot(5,1,2)
            if cond < 2
                psd_freq_no_training = psd_freq_2(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
            else
                psd_freq_no_training = psd_freq_2;
            end
            imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');
            caxis([low_bound high_bound])
            a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
            title("4-8Hz")
            xlabel("time [s]")
            ylabel("Channel depth [µm]")

            subplot(5,1,3)
            if cond < 2
                psd_freq_no_training = psd_freq_3(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
            else
                psd_freq_no_training = psd_freq_3;
            end
            imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');
            caxis([low_bound high_bound])
            a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
            title("8-12Hz")
            xlabel("time [s]")
            ylabel("Channel depth [µm]")

            subplot(5,1,4)
            if cond < 2
                psd_freq_no_training = psd_freq_4(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
            else
                psd_freq_no_training = psd_freq_4;
            end
            imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');
            caxis([low_bound high_bound])
            a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
            title("12-35Hz")
            xlabel("time [s]")
            ylabel("Channel depth [µm]")

            subplot(5,1,5)
            if cond < 2
                psd_freq_no_training = psd_freq_5(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
            else
                psd_freq_no_training = psd_freq_5;
            end
            imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');


            caxis([low_bound high_bound])
            a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
            title("35-200")
            xlabel("time [s]")
            ylabel("Channel depth [µm]")
            %% Save image

            % if cond == 2
            %     saveas(fig,"plots_Jeremy/psds_normalised/masters/"+lfpFilename(43:50)+".fig")
            %     saveas(fig,"plots_Jeremy/psds_normalised/masters/"+lfpFilename(43:50)+".emf")
            %     save("plots_Jeremy/psds_normalised/masters/"+lfpFilename(43:50)+".mat")
            % elseif cond == 3
            %     saveas(fig,"plots_Jeremy/psds_normalised/yokes/"+lfpFilename(43:50)+".fig")
            %     saveas(fig,"plots_Jeremy/psds_normalised/yokes/"+lfpFilename(43:50)+".emf")
            %     save("plots_Jeremy/psds_normalised/yokes/"+lfpFilename(43:50)+".mat")
            %end
        end
    end
end
%% resting or after
if ~exist('endrest','var')
    endrest = endtime;
end
rest = 1;

if rest == 0
    mintime = 1500;
    maxtime = 2400;
elseif cond >= 2
    maxtime = endtime;
else%You want the statistical analysis to start only after the training period clearly has stopped --> the stimulus should be constant --> look on figure where the stim stops and input it into prompt
    figure
    plot(events.onsets,1:length(events.onsets))
    xlabel('time[s]')
    ylabel('Accumulated stimulations')
    mintime = input('Where does resting period start start: ');%%%%check where actual resting starts by plotting events
    maxtime = endrest;
    endtime = timeframe(end)*60;
end


%% Display difference in frequency in between the non training period
Max=50;

clearvars Tab Tab_load
Tab=1:endtime*30;
Tab_load=zeros(endtime*30,1)';


Ev=find((events.onsets>0)|(events.onsets<1500 & events.onsets>900));
a=1;
for i=1:numel(Tab)

    if a>numel(Ev)
        break
    end

    if events.onsets(Ev(a))<i*(1/30)
    Tab_load(i-3:i+3)=1;
    a=a+1;
    else
    Tab_load(i)=0;
    end

end

%% Calculate histograms of activities 

freqBandsMat = cell2mat(freqBands);
figure
for i = 1:5
    eval(sprintf('%s = %s;','psd_band',strcat('psd_freq_',num2str(i))));
    subplot(5,1,i);
    h1=histogram(psd_band(1:900/binSize,:),50,'FaceColor','r');
    set(gca,'YScale','log')
    hold on
    histogram(psd_band(ceil(mintime/binSize):floor(maxtime/binSize),:),h1.BinEdges,'FaceColor','b');
    legend('Before','After')
end

freqBandsMat = cell2mat(freqBands);
figure
title("Difference in mean power over 150 channels")
for i = 1:5
    eval(sprintf('%s = %s;','psd_band',strcat('psd_freq_',num2str(i))));
    subplot(5,1,i);
    title(strcat(freqBandsMat(i*2-1),'-',freqBandsMat(i*2),'Hz'));
    histogram((mean(psd_band(ceil(mintime/binSize):floor(maxtime/binSize),:),1))-mean(psd_band(1:900/binSize,:),1),15);
    legend('After-Before')
end


%% Difference before / after all depth

% for j = 1:5
%     figure
%     
%     eval(sprintf('%s = %s;','psd_band',strcat('psd_freq_',num2str(j))));
%     
%     Start= mean(psd_band(2400/binSize:endtime/binSize,:),1);
%     End= mean(psd_band(1:900/binSize,:),1) ;
%     for i=1:150
%         plot(End(i),Start(i),'xb')
%         hold on
%     end
%     line([0 max(max(Start),max(End))],[0 max(max(Start),max(End))],'Color','k','LineWidth',1)
%     title(sprintf("mean %i-%iHz power before vs after",freqBandsMat(j*2-1),freqBandsMat(j*2)))
%     xlabel('Mean energy before')
%     ylabel('Mean energy after')
% 
% end
%% Subplots
figure

subplot(5,1,1)
End= mean(psd_freq_1(ceil(mintime/binSize):floor(maxtime/binSize),:),1);
Start= mean(psd_freq_1(1:900/binSize,:),1) ;
for i=1:150
    plot(Start(i),End(i),'xb')
    hold on
end
xlabel('Mean energy before')
ylabel('Mean energy after')
line([0 max(max(Start),max(End))],[0 max(max(Start),max(End))],'Color','k','LineWidth',1)

subplot(5,1,2)
End= mean(psd_freq_2(ceil(mintime/binSize):floor(maxtime/binSize),:),1);
Start= mean(psd_freq_2(1:900/binSize,:),1) ;
for i=1:150
    plot(Start(i),End(i),'xb')
    hold on
end
xlabel('Mean energy before')
ylabel('Mean energy after')
line([0 max(max(Start),max(End))],[0 max(max(Start),max(End))],'Color','k','LineWidth',1)

subplot(5,1,3)
End= mean(psd_freq_3(ceil(mintime/binSize):floor(maxtime/binSize),:),1);
Start= mean(psd_freq_3(1:900/binSize,:),1) ;
for i=1:150
    plot(Start(i),End(i),'xb')
    hold on
end
xlabel('Mean energy before')
ylabel('Mean energy after')
line([0 max(max(Start),max(End))],[0 max(max(Start),max(End))],'Color','k','LineWidth',1)

subplot(5,1,4)
End= mean(psd_freq_4(ceil(mintime/binSize):floor(maxtime/binSize),:),1);
Start= mean(psd_freq_4(1:900/binSize,:),1) ;
for i=1:150
    plot(Start(i),End(i),'xb')
    hold on
end
xlabel('Mean energy before')
ylabel('Mean energy after')
line([0 max(max(Start),max(End))],[0 max(max(Start),max(End))],'Color','k','LineWidth',1)

subplot(5,1,5)
End= mean(psd_freq_5(ceil(mintime/binSize):floor(maxtime/binSize),:),1);
Start= mean(psd_freq_5(1:900/binSize,:),1) ;
for i=1:150
    plot(Start(i),End(i),'xb')
    hold on
end
xlabel('Mean energy before')
ylabel('Mean energy after')
line([0 max(max(Start),max(End))],[0 max(max(Start),max(End))],'Color','k','LineWidth',1)





%%
figure
subplot(5,1,1)
Inc= mean(psd_freq_1(ceil(mintime/binSize):floor(maxtime/binSize),:),1) - mean(psd_freq_1(1:900/binSize,:),1) ;
[a, b]=ttest2(mean(psd_freq_1(ceil(mintime/binSize):floor(maxtime/binSize),:),1),mean(psd_freq_1(1:900/binSize,:),1));
if b<0.05
    C='g';
else
    C='r';
end
for i=1:150
    line([0 1],[0 Inc(i)]);
    hold on
end
line([0 1],[0 mean(Inc)],'LineWidth',2,'Color',C);

subplot(5,1,2)
Inc= mean(psd_freq_2(ceil(mintime/binSize):floor(maxtime/binSize),:),1) - mean(psd_freq_2(1:900/binSize,:),1) ;
[a, b]=ttest2(mean(psd_freq_2(ceil(mintime/binSize):floor(maxtime/binSize),:),1),mean(psd_freq_2(1:900/binSize,:),1));
if b<0.05
    C='g';
else
    C='r';
end
for i=1:150
    line([0 1],[0 Inc(i)]);
    hold on
end
line([0 1],[0 mean(Inc)],'LineWidth',2,'Color',C);

subplot(5,1,3)
Inc= mean(psd_freq_3(ceil(mintime/binSize):floor(maxtime/binSize),:),1) - mean(psd_freq_3(1:900/binSize,:),1) ;
[a, b]=ttest2(mean(psd_freq_3(ceil(mintime/binSize):floor(maxtime/binSize),:),1),mean(psd_freq_3(1:900/binSize,:),1));
if b<0.05
    C='g';
else
    C='r';
end
for i=1:150
    line([0 1],[0 Inc(i)])
    hold on
end
line([0 1],[0 mean(Inc)],'LineWidth',2,'Color',C)

subplot(5,1,4)
Inc= mean(psd_freq_4(ceil(mintime/binSize):floor(maxtime/binSize),:),1) - mean(psd_freq_4(1:900/binSize,:),1) ;
[a, b]=ttest2(mean(psd_freq_4(ceil(mintime/binSize):floor(maxtime/binSize),:),1),mean(psd_freq_4(1:900/binSize,:),1));
if b<0.05
    C='g';
else
    C='r';
end
for i=1:150
    line([0 1],[0 Inc(i)])
    hold on
end
line([0 1],[0 mean(Inc)],'LineWidth',2,'Color',C)


subplot(5,1,5)
Inc= mean(psd_freq_5(ceil(mintime/binSize):floor(maxtime/binSize),:),1) - mean(psd_freq_5(1:900/binSize,:),1) ;
[a, b]=ttest2(mean(psd_freq_5(ceil(mintime/binSize):floor(maxtime/binSize),:),1),mean(psd_freq_5(1:900/binSize,:),1));
if b<0.05
    C='g';
else
    C='r';
end
for i=1:150
    line([0 1],[0 Inc(i)])
    hold on
end
line([0 1],[0 mean(Inc)],'LineWidth',2,'Color',C)

