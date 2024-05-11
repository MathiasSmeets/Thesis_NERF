clear;clc;close all;
Fs = 2500;
frequencies = [135,255];


for i = setdiff(1:9,[2,8,9])
    disp(i)
    % do it in multiple parts
    filter = designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',135,'HalfPowerFrequency2',255,'SampleRate',2500,'DesignMethod','butter');
    for j = 1:5
        load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_"+i+"m.mat")
        size_data = size(data,2);
        cur_data = data(:,1+(j-1)*size_data:j*size_data);
        clearvars data
        cur_filtered_data = filtfilt(filter,cur_data');
        save("/scratch/mathiass-takeokalab/01/LF_filtered_rec"+i+"_part"+j+"m.mat", "filtered_data", "-v7.3");
    end
    disp("baseline")

    load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_baseline_"+i+"m.mat")
    filtered_data = filtfilt(filter,data');
    save("/scratch/mathiass-takeokalab/01/LF_baseline_filtered_"+i+"m.mat", "filtered_data", "-v7.3");
end