clear;clc;close all;
Fs = 2500;
frequencies = [135,255];


for i = setdiff(1:9,[2,8,9])
    disp(i)
    % do it in multiple parts
    filter = designfilt('bandpassiir','FilterOrder',20,'PassbandFrequency1',11,'PassbandFrequency2',16,'StopbandAttenuation1',60,'PassbandRipple',1,'StopbandAttenuation2',60,'SampleRate',2500,'DesignMethod','ellip');
    for j = 1:5
        disp(i+": "+j)
        load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_"+i+"m.mat")
        size_data = floor(size(data,2)/5);
        cur_data = data(:,1+(j-1)*size_data:j*size_data);
        clearvars data
        cur_filtered_data = filtfilt(filter,cur_data');
        save("/scratch/mathiass-takeokalab/01/LF_filtered_rec"+i+"_part"+j+"m_spindle.mat", "cur_filtered_data", "-v7.3");
    end
    disp("baseline")

    load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_baseline_"+i+"m.mat")
    filtered_data = filtfilt(filter,data');
    save("/scratch/mathiass-takeokalab/01/LF_baseline_filtered_"+i+"m.mat", "filtered_data", "-v7.3");
end