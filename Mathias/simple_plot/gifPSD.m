% Make a GIF of a wanted time period (timeband variable)
clear
close all

if ~exist('init','var')
    clc
    
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
    
    binSize = 10;
    % - - - - - Initialization of the time bands you're interested in the
    % PSD, a "binSize [s]" time frame after this time will be plotted. The different
    % rows are for different recordings
    timeBands = [5500 5800 6050 6180; 0 0 0 0; 3850 3970 4050 4243];%rest upactivation midactivation highactivity

    % Extract the recording of interest
    % - - - - - If you want to do a for loop you can use some commands like
    % - - - - - this
    
    %RecordingNumber=28;
    %spk = loadKSdir(sort_masters_horridge{RecordingNumber});
    
    % - - - - - - - To study a unique recording
    possiblePaths=['..\..\Mattia\Spinal_Cord\Recordings\20220725\20220725_long_long_recording_g0';
        '..\..\Mattia\Spinal_Cord\Recordings\20220819\20220819_long_long_recording_g0';%Recording seems to be corrupted or so, only 606 seconds long
        '..\..\Mattia\Spinal_Cord\Recordings\20220831\20220831_long_long_recording_g0';
        '..\..\Mattia\Spinal_Cord\Recordings\20220921\20220921_long_long_recording_g0';
        '..\..\Mattia\Spinal_Cord\Recordings\20220923\20220923_long_long_recording_g0';
        '..\..\Mattia\Spinal_Cord\Recordings\20221011\20221011_long_long_recording_g0'];
    possiblePathsYokes =['..\..\Mattia\Spinal_Cord\Recordings\20220726\20220726_long_long_recording_g0';
        '..\..\Mattia\Spinal_Cord\Recordings\20220905\20220905_long_long_recording_g0';
        '..\..\Mattia\Spinal_Cord\Recordings\20220906\20220906_long_long_recording_g0';%No events.mat
        '..\..\Mattia\Spinal_Cord\Recordings\20220922\20220922_long_long_recording_g0';
        '..\..\Mattia\Spinal_Cord\Recordings\20220927\20220927_long_long_recording_g0';
        '..\..\Mattia\Spinal_Cord\Recordings\20221006\20221006_long_long_recording_g0'];   
    timeframes = [0 15 25 40 106 145 155 170;
            0 15 25 40 106 121 131 146;
            0 15 25 40 106 121 131 146;
            0 15 25 40 107 122 133 148;
            0 15 25 40 110 132 142 157;
            0 15 25 40 110 126 136 151];%times represent: start horridge after rest spont horridge after end
        
    timeframesYokes = [0 18 28 43 109 124 134 149;
            0 15 25 40 106 121 131 146;
            0 15 25 40 106 121 131 146;
            0 15 25 40 104 129 139 154;
            0 15 25 40 161 176 186 201;
            0 15 25 40 110 127 137 152];
    i = 3;
    yokes = 0;%pick mouse from yoke list or not
    if yokes
        myKsDir=possiblePathsYokes(i, :); 
        timeframe = timeframesYokes(i,:);
    else
        myKsDir=possiblePaths(i, :);
        timeframe = timeframes(i,:);
        timeBand = timeBands(i,:);
        gifstart = timeBand(1);
        gifend = timeBand(end);
    end
    %going to inspect the resting period
    timeframe(end)
    % - - - - - Load the file that correspond to the events 
    % - - - - - One event correspond to an electrical shock during the horridge
    % - - - - - paradigm
    load(fullfile(myKsDir,'events.mat'),'events');
    evnames = {events.type};
    evtid = find(contains(evnames,'horridge'));
    if isempty(evtid)==1
        evtid = find(contains(evnames,'opto'));
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
    h = figure
    for i=gifstart:binSize:gifend-binSize % to cover 0 to 2400 seconds
    %for i=225*25/binSize:275*25/binSize % to cover 0 to 2400 seconds
          
        % Creating the tbin of 10 seconds
        tbin = [i i+binSize];
        
        % - - - - - - - - - - - - - - -
        i %indicate in which clip you are
        
        % Exctracting the LFP
        % psdM are the LFP along ~4000Hz frequency band
        % FM is just the number of the frequency
        % allPowerVar is the variance for each frequecny - we don't use it
        [psdM, FM, allPowerVar] = psd_allchannels(myKsDir,'probe','Neuropixels3A','chsel',1:maxch,'time_range',tbin(1,:));
            
        
        %smooth psd profiles with gaussian-weighted moving average + convert to db
        psdplotM = smoothdata(10*log10(psdM),2,'gaussian','SmoothingFactor',.01);
%         figure
%         maxF = find(FM < 25,1, 'last');
%         plot(FM(1:maxF), psdM(1,1:maxF))
%         hold on
%         plot(FM(1:maxF), psdM(75,1:maxF))
% 
%         plot(FM(1:maxF), psdM(150,1:maxF))
%         legend("Ventral","medial","Dorsal")
        
        % - - - - - - - -plot results for all channels in gif format, only
        % 1 plot
%         F = FM(1,:);
%         %subplot(10,Column(i),i)
%         maxtime = binSize;
%         if i == gifstart
%             h = imagesc(F,1:size(psdM,1)*10,flip(psdM,1));
%             ax = ancestor(h, 'axes');
%             title(ax,sprintf("time = %is",i))
%             caxis([0 10])
%             %colormap('jet')
%             cbar = colorbar
%             xlabel('F [Hz]')
%             ylabel('depth on probe (µm)')
%             set(gca, 'xlim', [0 20],'ylim', [1 size(psdM,1)*10])
%             cbar.Label.String = 'PSD';
%             gif("psd.gif",'DelayTime',1/3,"LoopCount",1) 
%         else
%             title(ax,sprintf("time = %is",i))
%             set(h,'xdata',F,'ydata',1:size(psdM,1)*10,'cdata',flip(psdM,1))
%             gif
%         end

         % - - - - - - - -plot results for all channels in gif format in 2
         % subplots (left & right)
        F = FM(1,:);
        %subplot(10,Column(i),i)
        maxtime = binSize;
        
        subplot(1,2,1)
        imagesc(F,1:size(psdM,1)*10,flip(psdM,1));
        caxis([0 50])
        %colormap('jet')
        cbar = colorbar
        xlabel('F [Hz]')
        ylabel('depth on probe (µm)')
        set(gca, 'xlim', [0 5],'ylim', [1 size(psdM,1)*10])
        %title('Masters during adaptation')        
        cbar.Label.String = 'PSD [db/Hz]';

        subplot(1,2,2)
        imagesc(F,1:size(psdM,1)*10,flip(psdM,1));
        caxis([0 3])
        %colormap('jet')
        cbar = colorbar
        xlabel('F [Hz]')
        ylabel('depth on probe (µm)')
        set(gca, 'xlim', [5 25],'ylim', [1 size(psdM,1)*10])
        %title('Masters during adaptation')        
        cbar.Label.String = 'PSD [db/Hz]';
        filename = "psd_subplots.gif";
        

        drawnow
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if i == gifstart
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        end 
        
        
        
        
        
    end
end