% this script appends available frM or frY tables

clear; clc; close all;
%datapath = "\\nerffs13\takeokalabwip2020\Mathias\";
%volume_base2 = '/mnt/takeokalab/';
datapath = "/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/tables";
savepath = "/scratch/mathiass-takeokalab/01/";

frM_1 = load(fullfile(datapath, "frM_switched_1.mat"));
cur_table = struct2array(frM_1);
cur_table = cur_table(:,[1,4,7:9,11,13,15]);
clearvars frY_1

% append tables
for i = 2:11
    new_table = load(fullfile(datapath, "frM_switched_" + i + ".mat"));
    new_table = struct2array(new_table);

    % remove rows that contain zero in the recording column
    new_table(~new_table.Recording,:) = [];

    % append to total table
    cur_table = [cur_table; new_table(:,[1,4,7:9,11,13,15])];

    clearvars new_table
end

%save(savepath + "frY_peak_total.mat", "data", "-v7.3")

% append 3 parts of horridge into one
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
%save(fullfile(savepath, "horridge_data.mat"), "data", "-v7.3")
clearvars data

% get switch data in array with recording nb in front
switch_data = zeros(size(cur_table.Fr_switch,1), size(cur_table.Fr_switch{1,1},2)+1);
for i = 1:size(switch_data,1)
    cur_switch_data = [cur_table.Recording(i,1), cur_table.Fr_switch{i,1}];
    switch_data(i,:) = cur_switch_data;
end

%save(fullfile(savepath, "switch_data.mat"), "switch_data", "-v7.3")
clearvars switch_data
clearvars cur_switch_data

% get after data in array with recording nb in front
after_data = zeros(size(cur_table.Fr_after,1), size(cur_table.Fr_after{1,1},2)+1);
for i = 1:size(after_data,1)
    cur_after_data = [cur_table.Recording(i,1), cur_table.Fr_after{i,1}];
    after_data(i,:) = cur_after_data;
end

%save(fullfile(savepath, "after_data.mat"), "after_data", "-v7.3")
clearvars after_data

% before_data
before_data = zeros(size(cur_table.Fr_spont,1), size(cur_table.Fr_spont{1,1},2)+1);
for i = 1:size(before_data,1)
    cur_before_data = [cur_table.Recording(i,1), cur_table.Fr_spont{i,1}];
    before_data(i,:) = cur_before_data;
end

%save(fullfile(savepath, "before_data_y.mat"), "before_data", "-v7.3")
clearvars after_data

% waiting_data
waiting_data = cell(11,1);
for i = 1:size(waiting_data,1)
    row_indices = cur_table.Recording == i;
    cur_waiting_data = [cur_table.Fr_waiting(row_indices,1)];
    cur_mouse = i;
    waiting_data{cur_mouse} = cur_waiting_data;
end
save(fullfile(savepath, "waiting_data_m.mat"), "waiting_data", "-v7.3")
