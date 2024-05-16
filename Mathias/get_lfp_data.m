clear;close all;clc;
maxdepth = 1500;
depth_binning = 50;


freqBands = {[0.5 4], [4 8], [8 12], [12 35], [35 200]};

nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

% - - - - - - - To study a unique recording

% LEARNER
% possiblePaths={'/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220725/20220725_long_long_recording_g0';
%         '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220819/20220819_long_long_recording2_g0';%Recording seems to be corrupted or so, only 606 seconds long
%         '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220831/20220831_long_long_recording_g0';
%         '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220921/20220921_long_long_recording_g0';
%         '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220923/20220923_long_long_recording_g0';
%         '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20221011/20221011_long_long_recording_g0';
%         '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20221116/20221116_long_long_recording_g0';
%         '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20221124/20221124_long_long_recording2_g0';
%         '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20221128/20221128_long_long_recording_g1';
%         };
 
% CONTROL
possiblePaths={'/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220705/20220705_long_long_recording_g0';
        '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220725/20220725_long_long_recording_g0';%Recording seems to be corrupted or so, only 606 seconds long
        '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220905/20220905_long_long_recording_g0';
        '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220922/20220922_long_long_recording_g0';
        '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220927/20220927_long_long_recording_g0';
        '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20221006/20221006_long_long_recording_g0';
        '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20221018/20221018_long_long_recording_g0';
        '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220906/20220906_long_long_recording2_g0';
        '/mnt/takeokalab/takeokalabwip2020/past lab members/Mattia/Spinal_Cord/Recordings/20220413/20220413_long_long_recording_g1';
        };

rest_start_times = [3665,3659,3655,3709,3921,3940,3040,3075,3207];
rest_end_times = [8700,7270,7270,7382,7924,7563,6902,7212,7213];

for i = [2,8,9]
    disp(i)
    myKsDir=possiblePaths{i};
    Fs = 2500;

    %% 
    % Keep in mind that we record NORMALLY (it can change based on the recording) like this :
    % 0 to 900 seconds : rest time - beginning of the experiment, called "the before phase"
    % 900 to 1500 seconds : horridge paradigm - where shocks and learning happen
    % 1500 to 2400 seconds : rest time - called "the after phase"
    %% Plot data to inspect


    %sampStarts = Fs*rest_start_times(i);
    %sampEnds = Fs*rest_end_times(i);
    sampStarts = rest_start_times(i)*Fs;
    sampEnds = rest_end_times(i)*Fs;

    d = dir(fullfile(myKsDir, '*.lf.bin'));
    nSamps = d.bytes/2/nChansInFile;
    lfpFilename = fullfile(myKsDir, d(1).name);
    mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});

    data = double(mmf.Data(1).x(:, (sampStarts:sampEnds)));
    save("/scratch/mathiass-takeokalab/01/LF_"+i+"y.mat", "data", "-v7.3");
    clearvars data
    % filtered_data = zeros(length(freqBands),size(data,2),size(data,1));
    % for j = 1:length(freqBands)
    %     freqs = freqBands{j};
    %     filter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',freqs(1),'HalfPowerFrequency2',freqs(2),'SampleRate',2500,'DesignMethod','butter');
    %     filtered_data(j,:,:) = filtfilt(filter,data');
    % end
    % filtered_data = permute(filtered_data,[1 3 2]);
    
end