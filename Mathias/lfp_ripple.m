clear;clc;close all;
Fs = 2500;
window = 5*Fs/25;

for i = setdiff(1:9,[2,8,9])
    disp(i)
    load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_filtered_total"+i+"m.mat")
    data = downsample(data,25);
    ripple = zeros(size(data,1),1);
    excluded = [];
    for j = 5*Fs/25+1:size(data,1)-1
        disp(j+"/"+size(data,1))
        cur_data = data(setdiff(j-window:j,excluded),:);
        cur_power = cur_data.^2;
        cur_power = sum(cur_power,2);
        cur_average = mean(cur_power);
        cur_mad = mad(cur_power);
        
        if 5*cur_mad + cur_average < (sum(data(j+1,:).^2))
            ripple(j+1) = 1;
            excluded = [excluded, j+1];
        end
    end
    save("/scratch/mathiass-takeokalab/01/ripple_"+i+".mat", "ripple", "-v7.3");
end