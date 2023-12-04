% this script appends available frM or frY tables

clear; clc; close all;
%datapath = "\\nerffs13\takeokalabwip2020\Mathias\";
datapath = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
savepath = "/scratch/mathiass-takeokalab/01/";

frM_7 = load(datapath + "frY_7.mat");
cur_table = struct2array(frM_7);
cur_table = cur_table(:,[1,7:9]);
clearvars frY_7

i = 8;
while i <= 34 
    if i ~= 20 && i ~= 32 && i~=33 && i ~=34
        new_table = load(datapath + "frY_" + i + ".mat");
        new_table = struct2array(new_table);

        % remove rows that contain zero in the recording column
        new_table(~new_table.Recording,:) = [];
        
        % append to total table
        cur_table = [cur_table; new_table(:,[1,7:9])];
    end

    i = i + 1;
    clearvars new_table
end

data = zeros(size(cur_table,1), 1 + size(cur_table.Fr_early{1,1},2)+size(cur_table.Fr_middle{1,1},2)+size(cur_table.Fr_late{1,1},2));
for i = 1:size(data,1)
    cur_data = [cur_table.Recording(i,1) ,cur_table.Fr_early{i,1}, cur_table.Fr_middle{i,1}, cur_table.Fr_late{i,1}];
    if size(cur_data,2) > size(data,2)
        data(i,:) = cur_data(1,1:size(data,2));
    elseif size(cur_data,2) < size(data,2)
        data(i,1:size(cur_data,2)) = cur_data;
    else
        data(i,:) = cur_data;
    end
end
save(savepath + "frY_peak_total.mat", "data", "-v7.3")