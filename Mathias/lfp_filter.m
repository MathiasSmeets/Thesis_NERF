clear;clc;close all;
Fs = 2500;
frequencies = [135,255];


for i = setdiff(1:9,[2,8,9])
    disp(i)
    load("/mnt/takeokalab/takeokalabwip2020/Mathias/switch_data/LF_signals/LF_"+i+"m.mat")
    filter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',135,'HalfPowerFrequency2',255,'SampleRate',2500,'DesignMethod','butter');
    filtered_data = filtfilt(filter,data');
    save("/scratch/mathiass-takeokalab/01/LF_filtered_"+i+"m.mat", "filtered_data", "-v7.3");

    load("/mnt/takeokalab/takeokalabwip2020/Mathias/switch_data/LF_signals/LF_baseline_"+i+"m.mat")
    filtered_data = filtfilt(filter,data');
    save("/scratch/mathiass-takeokalab/01/LF_baseline_filtered_"+i+"m.mat", "filtered_data", "-v7.3");
end