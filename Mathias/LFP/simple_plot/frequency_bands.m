%
% Plot different frequency bands

%close all


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
freqBands = {[0.5 4]};%, [4 8], [8 12], [12 35], [35 200]};

% Extract the recording of interest
% - - - - - If you want to do a for loop you can use some commands like
% - - - - - this

%RecordingNumber=28;
%spk = loadKSdir(sort_masters_horridge{RecordingNumber});

% - - - - - - - To study a unique recording
possiblePaths=['Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220725\20220725_long_long_recording_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220819\20220819_long_long_recording_g0';%Recording seems to be corrupted or so, only 606 seconds long
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220831\20220831_long_long_recording_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220921\20220921_long_long_recording_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20220923\20220923_long_long_recording_g0';
        'Y:\past lab members\Mattia\Spinal_Cord\Recordings\20221011\20221011_long_long_recording_g0';
        ];
timeframes = [0 15 25 40 106 145 155 170;0 0 0 0 0 0 0 0;
            0 15 25 40 106 121 131 146];%get them in LFP_analysis

for i = 3:3%size(possiblePaths,1) 
    myKsDir=possiblePaths(i, :);
    timeframe = timeframes(i,:);
    
    %% 
    % Keep in mind that we record NORMALLY (it can change based on the recording) like this :
    % 0 to 900 seconds : rest time - beginning of the experiment, called "the before phase"
    % 900 to 1500 seconds : horridge paradigm - where shocks and learning happen
    % 1500 to 2400 seconds : rest time - called "the after phase"
    %% Plot data to inspect
    timeStart = 1500;%At how many seconds do you want to start
    time = timeframe(6)*60-1500;%how long do you want to see 225
    nch = 150;%pick how many channels you would like to select
    Fs = 2500;
    chsel = [1,75,150];
    ch = [1,75,150];%Choose channel you want to be plotted
    sampStarts = Fs*timeStart;

    d = dir(fullfile(myKsDir, '*.lf.bin'));
    nSamps = d.bytes/2/nChansInFile;
    lfpFilename = fullfile(myKsDir, d(1).name);
    mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});
    tic
    thinning = 100;
    data = double(mmf.Data(1).x(chsel, (1:time*Fs)+sampStarts));
    toc
    %% Don't include the training period
    data_no_training = data;%[data(1:900*Fs),data(1500*Fs:end)];
    
    %% MatLab lowpass filer
    filtered_data = zeros(length(freqBands),size(data,2),size(data,1));
    for j = 1:length(freqBands)
        freqs = freqBands{j};
        filter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',freqs(1),'HalfPowerFrequency2',freqs(2),'SampleRate',2500,'DesignMethod','butter');
        filtered_data(j,:,:) = filtfilt(filter,data_no_training');
    end
    filtered_data = permute(filtered_data,[1 3 2]);
    

%   
end
%% plot raw vs delta
figure
sgtitle(sprintf("Raw signal vs delta waves for channels 1,75,150"))
for j = 1:length(freqBands)
    freqs = freqBands{j};
    subplot(2,1,2)
    for k = 1:size(filtered_data,2)
        plot(timeStart:1/Fs:time+timeStart-1/Fs,reshape(filtered_data(j,k,:),[],1))
        hold on
    end
    title(sprintf("%.1f-%.1fHz",freqs(1),freqs(2)))
    xlabel("time [s]")
    ylabel("Voltage [µV]")
    legend(sprintf("Channel %i",chsel(1)),sprintf("Channel %i",chsel(2)),sprintf("Channel %i",chsel(3)))
end
subplot(2,1,1)
for k = 1:size(filtered_data,2)
    plot(timeStart:1/Fs:time+timeStart-1/Fs,data(k,:))
    hold on
end
title(sprintf("Raw signal"))
xlabel("time [s]")
ylabel("Voltage [µV]")
legend(sprintf("Channel %i",chsel(1)),sprintf("Channel %i",chsel(2)),sprintf("Channel %i",chsel(3)))
%%



% %% Compare LFP data when using multiple channels
% 
% time = 100;%how long do you want to see
% nch = 385;
% Fs = 2500;
% chsel = 1:nch;
% ch = 1;%Choose channel you want to be plotted
% timeStart = 0;%At how many seconds do you want to start
% sampStarts = Fs*timeStart;
% 
% d = dir(fullfile(myKsDir, '*.lf.bin'));
% nSamps = floor(d.bytes/2/nch);
% lfpFilename = fullfile(myKsDir, d(1).name);
% mmf = memmapfile(lfpFilename, 'Format', {'int16', [nch nSamps], 'x'});
% data_no_training = double(mmf.Data(1).x(chsel, (1:time*Fs)+sampStarts));
% figure
% plot(0:1/Fs:time-1/Fs,data_no_training(1,:))
% hold on
% hhhh = data_no_training;
% 
% nch = 1;
% Fs = 2500;
% chsel = 1:nch;
% ch = 1;%Choose channel you want to be plotted
% timeStart = 0;%At how many seconds do you want to start
% sampStarts = Fs*timeStart;
% 
% d = dir(fullfile(myKsDir, '*.lf.bin'));
% nSamps = floor(d.bytes/2/nch);
% lfpFilename = fullfile(myKsDir, d(1).name);
% mmf = memmapfile(lfpFilename, 'Format', {'int16', [nch nSamps], 'x'});
% data_no_training = double(mmf.Data(1).x(chsel, (1:time*Fs)+sampStarts));
% plot(0:1/Fs:time-1/Fs,data_no_training(1,:))
% 
% legend('multiple','one')
%  