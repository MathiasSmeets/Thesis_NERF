clear;clc;close all;

for i = setdiff(3:9,[2,8,9])
    disp(i)
    % do it in multiple parts
    data = [];
    for j = 1:5
        disp(j)
        load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_filtered_rec"+i+"_part"+j+"m_ripple2.mat")
        data = [data;cur_filtered_data];
    end
    save("/scratch/mathiass-takeokalab/01/LF_filtered_total"+i+"m_ripple2.mat", "data", "-v7.3");
end