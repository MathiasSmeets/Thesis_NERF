clear;clc;close all;
Fs = 2500;

for i = 1:9
    disp(i)
    %load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_filtered_total"+i+"m.mat")
    data = [];
    for j = 1:5
        disp("prep: " + j)
        load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_filtered_rec"+i+"_part"+j+"m_ripple2.mat")
        data = [data;cur_filtered_data];
    end
    % calculate power
    cur_power = data.^2;
    cur_power = sum(cur_power,2);

    save("/scratch/mathiass-takeokalab/01/power_"+i+"m.mat", "power", "-v7.3");
end