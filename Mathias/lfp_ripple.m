clear;clc;close all;
Fs = 2500;
window = 5*Fs;

for i = setdiff(1:9,[2,8,9])
    disp(i)
    load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_filtered_total"+i+"m.mat")
    ripple = zeros(size(data,2),1);
    for j = 5*Fs:size(data,1)-1
        cur_data = data(j-window:j,:);
        cur_power = cur_data.^2;
        cur_power = sum(cur_power,2);
        cur_average = mean(cur_power);
        cur_mad = mad(cur_power);
        
        if 5*cur_mad + cur_average < (sum(data(j+1,:).^2))
            ripple(j+1) = 1;
        end
    end
    save()
end