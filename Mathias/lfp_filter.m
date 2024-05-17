clear;clc;close all;
Fs = 2500;
frequencies = [135,255];


for i = 8
    disp(i)
    % do it in multiple parts
    %filter = designfilt('bandpassiir','FilterOrder',20,'PassbandFrequency1',11,'PassbandFrequency2',16,'StopbandAttenuation1',60,'PassbandRipple',1,'StopbandAttenuation2',60,'SampleRate',2500,'DesignMethod','ellip');
    filter_ripple2 = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',100,'HalfPowerFrequency2',250,'SampleRate',2500,'DesignMethod','butter');
    %filter_spindle = designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',11,'HalfPowerFrequency2',20,'SampleRate',2500,'DesignMethod','butter');
    for j = 1:5%1:5
        disp(i+": "+j)
        load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_"+i+"y.mat")
        size_data = floor(size(data,2)/5);
        cur_data = data(:,1+(j-1)*size_data:j*size_data);
        clearvars data
        cur_filtered_data = filtfilt(filter_ripple2,cur_data');
        clearvars cur_data
        save("/scratch/mathiass-takeokalab/01/LF_filtered_rec"+i+"_part"+j+"y_ripple2.mat", "cur_filtered_data", "-v7.3");
    end
end