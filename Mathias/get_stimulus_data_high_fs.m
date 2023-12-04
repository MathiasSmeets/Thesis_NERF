clear; clc; close all;

%% load data
%path = "\\nerffs13\takeokalabwip2020\Mathias\10kfs\";
path = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
write_path = "/scratch/mathiass-takeokalab/01/";

data_m = load(path + "frM_peak_total.mat");
data_m = struct2array(data_m);
data_y = load(path + "frY_peak_total.mat");
data_y = struct2array(data_y);

new_data_m = zeros(size(data_m.zscore_early_full,1),1+size(data_m.zscore_early_full{1,1},2) + size(data_m.zscore_mid_full{1,1},2) + size(data_m.zscore_late_full{1,1},2));
new_data_y = zeros(size(data_y.zscore_early_full,1),1+size(data_y.zscore_early_full{1,1},2) + size(data_y.zscore_mid_full{1,1},2) + size(data_y.zscore_late_full{1,1},2));

data_m(~data_m.Recording,:) = [];
data_y(~data_y.Recording,:) = [];
for i = 1:15
    new_data_m(i,:) = [data_m.Recording(i,1) ,data_m.Fr_early{i,1}, data_m.Fr_middle{i,1}, data_m.Fr_late{i,1}];
    new_data_y(i,:) = [data_y.Recording(i,1) ,data_y.Fr_early{i,1}, data_y.Fr_middle{i,1}, data_y.Fr_late{i,1}];
end

data_m = new_data_m;
data_y = new_data_y;

event_main_M = struct();
event_main_M.event_7 = load(path + "\eventM_7.mat", "-mat", "events");
event_main_M.event_8 = load(path + "\eventM_8.mat", "-mat", "events");
event_main_M.event_9 = load(path + "\eventM_9.mat", "-mat", "events");
event_main_M.event_10 = load(path + "\eventM_10.mat", "-mat", "events");
event_main_M.event_11 = load(path + "\eventM_11.mat", "-mat", "events");
event_main_M.event_12 = load(path + "\eventM_12.mat", "-mat", "events");
event_main_M.event_13 = load(path + "\eventM_13.mat", "-mat", "events");
event_main_M.event_14 = load(path + "\eventM_14.mat", "-mat", "events");
event_main_M.event_15 = load(path + "\eventM_15.mat", "-mat", "events");
event_main_M.event_16 = load(path + "\eventM_16.mat", "-mat", "events");
event_main_M.event_17 = load(path + "\eventM_17.mat", "-mat", "events");
event_main_M.event_18 = load(path + "\eventM_18.mat", "-mat", "events");
event_main_M.event_19 = load(path + "\eventM_19.mat", "-mat", "events");
%event_20 = load(path + "\eventM_20.mat");
event_main_M.event_21 = load(path + "\eventM_21.mat", "-mat", "events");
event_main_M.event_22 = load(path + "\eventM_22.mat", "-mat", "events");
event_main_M.event_23 = load(path + "\eventM_23.mat", "-mat", "events");
event_main_M.event_24 = load(path + "\eventM_24.mat", "-mat", "events");
event_main_M.event_25 = load(path + "\eventM_25.mat", "-mat", "events");
event_main_M.event_26 = load(path + "\eventM_26.mat", "-mat", "events");
event_main_M.event_27 = load(path + "\eventM_27.mat", "-mat", "events");
event_main_M.event_28 = load(path + "\eventM_28.mat", "-mat", "events");
event_main_M.event_29 = load(path + "\eventM_29.mat", "-mat", "events");
event_main_M.event_30 = load(path + "\eventM_30.mat", "-mat", "events");
event_main_M.event_31 = load(path + "\eventM_31.mat", "-mat", "events");
%event_32 = load(path + "\eventM_32.mat");
%event_33 = load(path + "\eventM_33.mat");
%event_34 = load(path + "\eventM_34.mat");


event_main_Y = struct();
event_main_Y.event_7 = load(path + "\eventY_7.mat", "-mat", "events");
event_main_Y.event_8 = load(path + "\eventY_8.mat", "-mat", "events");
event_main_Y.event_9 = load(path + "\eventY_9.mat", "-mat", "events");
event_main_Y.event_10 = load(path + "\eventY_10.mat", "-mat", "events");
event_main_Y.event_11 = load(path + "\eventY_11.mat", "-mat", "events");
event_main_Y.event_12 = load(path + "\eventY_12.mat", "-mat", "events");
event_main_Y.event_13 = load(path + "\eventY_13.mat", "-mat", "events");
event_main_Y.event_14 = load(path + "\eventY_14.mat", "-mat", "events");
event_main_Y.event_15 = load(path + "\eventY_15.mat", "-mat", "events");
event_main_Y.event_16 = load(path + "\eventY_16.mat", "-mat", "events");
event_main_Y.event_17 = load(path + "\eventY_17.mat", "-mat", "events");
event_main_Y.event_18 = load(path + "\eventY_18.mat", "-mat", "events");
event_main_Y.event_19 = load(path + "\eventY_19.mat", "-mat", "events");
%event_20 = load(path + "\eventY_20.mat");
event_main_Y.event_21 = load(path + "\eventY_21.mat", "-mat", "events");
event_main_Y.event_22 = load(path + "\eventY_22.mat", "-mat", "events");
event_main_Y.event_23 = load(path + "\eventY_23.mat", "-mat", "events");
event_main_Y.event_24 = load(path + "\eventY_24.mat", "-mat", "events");
event_main_Y.event_25 = load(path + "\eventY_25.mat", "-mat", "events");
event_main_Y.event_26 = load(path + "\eventY_26.mat", "-mat", "events");
event_main_Y.event_27 = load(path + "\eventY_27.mat", "-mat", "events");
event_main_Y.event_28 = load(path + "\eventY_28.mat", "-mat", "events");
event_main_Y.event_29 = load(path + "\eventY_29.mat", "-mat", "events");
event_main_Y.event_30 = load(path + "\eventY_30.mat", "-mat", "events");
event_main_Y.event_31 = load(path + "\eventY_31.mat", "-mat", "events");
%event_32 = load(path + "\eventY_32.mat");
%event_33 = load(path + "\eventY_33.mat");
%event_34 = load(path + "\eventY_34.mat");

%% loop over each event file to get data after every stimulation

fields = fieldnames(event_main_M);
numbers = setdiff(7:34, [20, 32, 33, 34]);
%numbers = 10;

%after_stimulus_data_m = cell(length(numbers), 5217); % go to get_max_onsets.m to determine this 5217
%after_stimulus_data_y = cell(length(numbers), 5217);
%after_stimulus_data_m = cell(length(numbers), 157);
%after_stimulus_data_y = cell(length(numbers), 157);
after_stimulus_data_m = cell(1,3371);
after_stimulus_data_y = cell(1,2648);

start_neuron_m = 1;
start_neuron_y = 1;
[EventSizeM, EventSizeY] = get_start_end_events();
for i = 1:numel(fields)
    nb_mouse = numbers(i);

    cur_onsets_m = event_main_M.(fields{i}).events.onsets;
    cur_onsets_y = event_main_Y.(fields{i}).events.onsets;
    cur_onsets_m = cur_onsets_m(EventSizeM(nb_mouse,1):EventSizeM(nb_mouse,2));
    cur_onsets_y = cur_onsets_y(EventSizeY(nb_mouse,1):EventSizeY(nb_mouse,2));
    cur_onsets_m = cur_onsets_m - cur_onsets_m(1,1);
    cur_onsets_y = cur_onsets_y - cur_onsets_y(1,1);

    % Make array with 10ms before and 60ms after each stimulus
    cur_nb_neurons_m = sum(data_m(:,1) == nb_mouse);
    cur_nb_neurons_y = sum(data_y(:,1) == nb_mouse);

    for j = 2:length(cur_onsets_m)
        if round(cur_onsets_m(j)*10000)+599 > 6000000
            after_stimulus_data_m{i,j-1} = data_m(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1,round(cur_onsets_m(j)*1000)*10-100+2:end);
        else
            after_stimulus_data_m{i,j-1} = data_m(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1,round(cur_onsets_m(j)*1000)*10-100+2:round(cur_onsets_m(j)*1000)*10+599+2);
        end
    end
    for j = 2:length(cur_onsets_y)
        if round(cur_onsets_y(j)*10000)+599 > 6000000
            after_stimulus_data_y{i,j-1} = data_y(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1,round(cur_onsets_y(j)*1000)*10-100+2:end);
        else
            after_stimulus_data_y{i,j-1} = data_y(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1,round(cur_onsets_y(j)*1000)*10-100+2:round(cur_onsets_y(j)*1000)*10+599+2);
        end
    end
    start_neuron_m = start_neuron_m + cur_nb_neurons_m;
    start_neuron_y = start_neuron_y + cur_nb_neurons_y;
end


save(write_path + "after_stimulus_data_m.mat", "after_stimulus_data_m", "-v7.3")
save(write_path + "after_stimulus_data_y.mat", "after_stimulus_data_y", "-v7.3")








