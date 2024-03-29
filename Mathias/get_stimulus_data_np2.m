clear; clc; close all;

if ispc
    volume_base = 'X:\';
    volume_base2 = 'X:\';
    data_path = fullfile(volume_base,"Mathias\switch_data\tabled_data");
    stimulus_path = fullfile(volume_base, "Mathias\switch_data\tables");
    save_path = data_path;
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/takeokalabwip2023';
    volume_base2 = '/mnt/takeokalab/takeokalabwip2023';
    data_path = fullfile(volume_base,"Mathias/switch_data/tabled_data");
    stimulus_path = fullfile(volume_base, "Mathias/switch_data/tables");
    save_path = "/scratch/mathiass-takeokalab/01";
end


%% load data

%path = "/mnt/takeokalab/takeokalabwip2023/Mathias/";



% data_m = load(fullfile(data_path,"horridge_data_m.mat"));
% data_m = struct2array(data_m);
data_y = load(fullfile(data_path, "horridge_data_np2.mat"));
data_y = struct2array(data_y);
%data_y = uint16(data_y);


% event_main_M = struct();
% event_main_M.event_1 = load(fullfile(stimulus_path, "frM_stim_switched_1.mat"), "-mat", "frM_stim");
% event_main_M.event_2 = load(fullfile(stimulus_path, "frM_stim_switched_2.mat"), "-mat", "frM_stim");
% event_main_M.event_3 = load(fullfile(stimulus_path, "frM_stim_switched_3.mat"), "-mat", "frM_stim");
% event_main_M.event_4 = load(fullfile(stimulus_path, "frM_stim_switched_4.mat"), "-mat", "frM_stim");
% event_main_M.event_5 = load(fullfile(stimulus_path, "frM_stim_switched_5.mat"), "-mat", "frM_stim");
% event_main_M.event_6 = load(fullfile(stimulus_path, "frM_stim_switched_6.mat"), "-mat", "frM_stim");
% event_main_M.event_7 = load(fullfile(stimulus_path, "frM_stim_switched_7.mat"), "-mat", "frM_stim");
% event_main_M.event_8 = load(fullfile(stimulus_path, "frM_stim_switched_8.mat"), "-mat", "frM_stim");
% event_main_M.event_9 = load(fullfile(stimulus_path, "frM_stim_switched_9.mat"), "-mat", "frM_stim");
% event_main_M.event_10 = load(fullfile(stimulus_path, "frM_stim_switched_10.mat"), "-mat", "frM_stim");
% event_main_M.event_11 = load(fullfile(stimulus_path, "frM_stim_switched_11.mat"), "-mat", "frM_stim");

event_main_Y = struct();
event_main_Y.event_1 = load(fullfile(stimulus_path, "frM_np2_stim_1.mat"), "-mat", "frM_stim");
event_main_Y.event_2 = load(fullfile(stimulus_path, "frM_np2_stim_2.mat"), "-mat", "frM_stim");
event_main_Y.event_3 = load(fullfile(stimulus_path, "frM_np2_stim_3.mat"), "-mat", "frM_stim");
% event_main_Y.event_4 = load(fullfile(stimulus_path, "frY_stim_switched_4.mat"), "-mat", "frM_stim");
% event_main_Y.event_5 = load(fullfile(stimulus_path, "frY_stim_switched_5.mat"), "-mat", "frM_stim");
% event_main_Y.event_6 = load(fullfile(stimulus_path, "frY_stim_switched_6.mat"), "-mat", "frM_stim");
% event_main_Y.event_7 = load(fullfile(stimulus_path, "frY_stim_switched_7.mat"), "-mat", "frM_stim");
% event_main_Y.event_8 = load(fullfile(stimulus_path, "frY_stim_switched_8.mat"), "-mat", "frM_stim");
% event_main_Y.event_9 = load(fullfile(stimulus_path, "frY_stim_switched_9.mat"), "-mat", "frM_stim");
% event_main_Y.event_10 = load(fullfile(stimulus_path, "frY_stim_switched_10.mat"), "-mat", "frM_stim");
% event_main_Y.event_11 = load(fullfile(stimulus_path, "frY_stim_switched_11.mat"), "-mat", "frM_stim");

% event_main_M.event_7 = load(path + "\eventM_7.mat", "-mat", "events");
% event_main_M.event_8 = load(path + "\eventM_8.mat", "-mat", "events");
% event_main_M.event_9 = load(path + "\eventM_9.mat", "-mat", "events");
% event_main_M.event_10 = load(path + "\eventM_10.mat", "-mat", "events");
% event_main_M.event_11 = load(path + "\eventM_11.mat", "-mat", "events");
% event_main_M.event_12 = load(path + "\eventM_12.mat", "-mat", "events");
% event_main_M.event_13 = load(path + "\eventM_13.mat", "-mat", "events");
% event_main_M.event_14 = load(path + "\eventM_14.mat", "-mat", "events");
% event_main_M.event_15 = load(path + "\eventM_15.mat", "-mat", "events");
% event_main_M.event_16 = load(path + "\eventM_16.mat", "-mat", "events");
% event_main_M.event_17 = load(path + "\eventM_17.mat", "-mat", "events");
% event_main_M.event_18 = load(path + "\eventM_18.mat", "-mat", "events");
% event_main_M.event_19 = load(path + "\eventM_19.mat", "-mat", "events");
% %event_20 = load(path + "\eventM_20.mat");
% event_main_M.event_21 = load(path + "\eventM_21.mat", "-mat", "events");
% event_main_M.event_22 = load(path + "\eventM_22.mat", "-mat", "events");
% event_main_M.event_23 = load(path + "\eventM_23.mat", "-mat", "events");
% event_main_M.event_24 = load(path + "\eventM_24.mat", "-mat", "events");
% event_main_M.event_25 = load(path + "\eventM_25.mat", "-mat", "events");
% event_main_M.event_26 = load(path + "\eventM_26.mat", "-mat", "events");
% event_main_M.event_27 = load(path + "\eventM_27.mat", "-mat", "events");
% event_main_M.event_28 = load(path + "\eventM_28.mat", "-mat", "events");
% event_main_M.event_29 = load(path + "\eventM_29.mat", "-mat", "events");
% event_main_M.event_30 = load(path + "\eventM_30.mat", "-mat", "events");
% event_main_M.event_31 = load(path + "\eventM_31.mat", "-mat", "events");
% %event_32 = load(path + "\eventM_32.mat");
% %event_33 = load(path + "\eventM_33.mat");
% %event_34 = load(path + "\eventM_34.mat");
% 
% 
% event_main_Y = struct();
% event_main_Y.event_7 = load(path + "\eventY_7.mat", "-mat", "events");
% event_main_Y.event_8 = load(path + "\eventY_8.mat", "-mat", "events");
% event_main_Y.event_9 = load(path + "\eventY_9.mat", "-mat", "events");
% event_main_Y.event_10 = load(path + "\eventY_10.mat", "-mat", "events");
% event_main_Y.event_11 = load(path + "\eventY_11.mat", "-mat", "events");
% event_main_Y.event_12 = load(path + "\eventY_12.mat", "-mat", "events");
% event_main_Y.event_13 = load(path + "\eventY_13.mat", "-mat", "events");
% event_main_Y.event_14 = load(path + "\eventY_14.mat", "-mat", "events");
% event_main_Y.event_15 = load(path + "\eventY_15.mat", "-mat", "events");
% event_main_Y.event_16 = load(path + "\eventY_16.mat", "-mat", "events");
% event_main_Y.event_17 = load(path + "\eventY_17.mat", "-mat", "events");
% event_main_Y.event_18 = load(path + "\eventY_18.mat", "-mat", "events");
% event_main_Y.event_19 = load(path + "\eventY_19.mat", "-mat", "events");
% %event_20 = load(path + "\eventY_20.mat");
% event_main_Y.event_21 = load(path + "\eventY_21.mat", "-mat", "events");
% event_main_Y.event_22 = load(path + "\eventY_22.mat", "-mat", "events");
% event_main_Y.event_23 = load(path + "\eventY_23.mat", "-mat", "events");
% event_main_Y.event_24 = load(path + "\eventY_24.mat", "-mat", "events");
% event_main_Y.event_25 = load(path + "\eventY_25.mat", "-mat", "events");
% event_main_Y.event_26 = load(path + "\eventY_26.mat", "-mat", "events");
% event_main_Y.event_27 = load(path + "\eventY_27.mat", "-mat", "events");
% event_main_Y.event_28 = load(path + "\eventY_28.mat", "-mat", "events");
% event_main_Y.event_29 = load(path + "\eventY_29.mat", "-mat", "events");
% event_main_Y.event_30 = load(path + "\eventY_30.mat", "-mat", "events");
% event_main_Y.event_31 = load(path + "\eventY_31.mat", "-mat", "events");
%event_32 = load(path + "\eventY_32.mat");
%event_33 = load(path + "\eventY_33.mat");
%event_34 = load(path + "\eventY_34.mat");

%% loop over each event file to get data after every stimulation


fs = 1; % ms
multiplier = 1/fs;

fields = fieldnames(event_main_Y);
%numbers = setdiff(7:34, [20, 32, 33, 34]);
numbers = 1:3;


%after_stimulus_data_m = cell(length(numbers), 157);
%after_stimulus_data_y = cell(length(numbers), 157);

%
start_neuron_y = 1;
%[EventSizeM, EventSizeY] = get_start_end_events();
%[~, EventSizeY, ~, maxy, ~, waitingeventY, ~, switcheventY] = get_start_end_events_np2(event_main_Y);
% after_stimulus_data_m = cell(length(numbers), maxm); % go to get_max_onsets.m to determine this 5217
maxy=1;
after_stimulus_data_y = cell(length(numbers), maxy);

for i = 1:numel(fields)
    
    nb_mouse = numbers(i);

    %cur_onsets_m = event_main_M.(fields{i}).events.onsets;
    %cur_onsets_y = event_main_Y.(fields{i}).events.onsets;
    % cur_onsets_m = event_main_M.(fields{i}).frM_stim.horridge{1,1};
    cur_onsets_y = event_main_Y.(fields{i}).frM_stim.horridge{1,1};
    % cur_onsets_m = cur_onsets_m;%(EventSizeM(nb_mouse,1):EventSizeM(nb_mouse,2));
    % cur_onsets_y = cur_onsets_y;%(EventSizeY(nb_mouse,1):EventSizeY(nb_mouse,2));
    % cur_onsets_m = cur_onsets_m - cur_onsets_m(1,1);
    cur_onsets_y = cur_onsets_y - cur_onsets_y(1,1);

    % Make array with 10ms before and 60ms after each stimulus
    % cur_nb_neurons_m = sum(data_m(:,1) == nb_mouse);
    cur_nb_neurons_y = sum(data_y(:,1) == nb_mouse);

    % for j = 2:length(cur_onsets_m)
    %     if round(cur_onsets_m(j)*1000)+59 > 600000
    %         after_stimulus_data_m{i,j-1} = data_m(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1,round(cur_onsets_m(j)*1000)-10+2:end);
    %         %after_stimulus_data_m{i,j-1}(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1, round(cur_onsets_m(j)*1000)+1:round(cur_onsets_m(j)*1000)+1) = 0;
    %     else
    %         after_stimulus_data_m{i,j-1} = data_m(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1,round(cur_onsets_m(j)*1000)-10+2:round(cur_onsets_m(j)*1000)+59+2);
    %         %after_stimulus_data_m{i,j-1}(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1, round(cur_onsets_m(j)*1000)+1:round(cur_onsets_m(j)*1000)+1) = 0;
    %     end
    % end
    for j = 2:length(cur_onsets_y)
        if round(cur_onsets_y(j)*1000)+59 > 600000
            after_stimulus_data_y{i,j-1} = data_y(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1,round(cur_onsets_y(j)*1000)-10+2:end);
            %after_stimulus_data_y{i,j-1}(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1, round(cur_onsets_y(j)*1000)+1:round(cur_onsets_y(j)*1000)+1) = 0;
        else
            after_stimulus_data_y{i,j-1} = data_y(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1,round(cur_onsets_y(j)*1000)-10+2:round(cur_onsets_y(j)*1000)+59+2);
            %after_stimulus_data_y{i,j-1}(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1, round(cur_onsets_y(j)*1000)+1:round(cur_onsets_y(j)*1000)+1) = 0;
        end
    end
    %data_y([start_neuron_y:start_neuron_y+cur_nb_neurons_y-1],:) = [];
    % start_neuron_m = start_neuron_m + cur_nb_neurons_m;
    start_neuron_y = start_neuron_y + cur_nb_neurons_y;
end
%save(fullfile(save_path, "after_stimulus_data_m_horridge.mat"), "after_stimulus_data_m", "-v7.3")
save(fullfile(save_path, "after_stimulus_data_np2_horridge.mat"), "after_stimulus_data_y", "-v7.3")
clearvars after_stimulus_data_m_horridge data_m after_stimulus_data_y_horridge data_y

%% now do the same for switch data
disp("second part")
%switch_m = load(fullfile(data_path, "switch_data_m.mat"));
%switch_m = struct2array(switch_m);
switch_y = load(fullfile(data_path, "switch_data_y.mat"));
switch_y = struct2array(switch_y);
%switch_y = uint16(switch_y);

start_neuron_m = 1;
start_neuron_y = 1;
%after_stimulus_switch_m = cell(length(numbers), 1); %% to do
after_stimulus_switch_y = cell(length(numbers), maxy);


for i = 1:numel(fields)
    nb_mouse = numbers(i);

    %cur_onsets_m = event_main_M.(fields{i}).events.onsets;
    %cur_onsets_y = event_main_Y.(fields{i}).events.onsets;
    %cur_onsets_m = event_main_M.(fields{i}).frM_stim.switch{1,1};
    cur_onsets_y = event_main_Y.(fields{i}).frM_stim.switch{1,1};
    %cur_onsets_m = cur_onsets_m;%(EventSizeM(nb_mouse,1):EventSizeM(nb_mouse,2));
    cur_onsets_y = cur_onsets_y;%(EventSizeY(nb_mouse,1):EventSizeY(nb_mouse,2));
    %cur_onsets_m = cur_onsets_m - cur_onsets_m(1,1);
    cur_onsets_y = cur_onsets_y - cur_onsets_y(1,1);

    % Make array with 10ms before and 60ms after each stimulus
    %cur_nb_neurons_m = sum(switch_m(:,1) == nb_mouse);
    cur_nb_neurons_y = sum(switch_y(:,1) == nb_mouse);

    % for j = 2:length(cur_onsets_m)
    %     if round(cur_onsets_m(j)*1000)+59 > 600000
    %         after_stimulus_switch_m{i,j-1} = switch_m(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1,round(cur_onsets_m(j)*1000)-10+2:end);
    %         after_stimulus_data_m{i,j-1}(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1, round(cur_onsets_m(j)*1000)+1:round(cur_onsets_m(j)*1000)+1) = 0;
    %     else
    %         after_stimulus_switch_m{i,j-1} = switch_m(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1,round(cur_onsets_m(j)*1000)-10+2:round(cur_onsets_m(j)*1000)+59+2);
    %         after_stimulus_data_m{i,j-1}(start_neuron_m:start_neuron_m + cur_nb_neurons_m - 1, round(cur_onsets_m(j)*1000)+1:round(cur_onsets_m(j)*1000)+1) = 0;
    %     end
    % end
    for j = 2:length(cur_onsets_y)
        if round(cur_onsets_y(j)*1000)+59 > 600000
            after_stimulus_switch_y{i,j-1} = switch_y(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1,round(cur_onsets_y(j)*1000)-10+2:end);
            %after_stimulus_switch_y{i,j-1}(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1, round(cur_onsets_y(j)*1000)+1:round(cur_onsets_y(j)*1000)+1) = 0;
        else
            after_stimulus_switch_y{i,j-1} = switch_y(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1,round(cur_onsets_y(j)*1000)-10+2:round(cur_onsets_y(j)*1000)+59+2);
            %after_stimulus_switch_y{i,j-1}(start_neuron_y:start_neuron_y + cur_nb_neurons_y - 1, round(cur_onsets_y(j)*1000)+1:round(cur_onsets_y(j)*1000)+1) = 0;
        end
    end
    %start_neuron_m = start_neuron_m + cur_nb_neurons_m;
    start_neuron_y = start_neuron_y + cur_nb_neurons_y;
end


%save(fullfile(save_path, "after_stimulus_data_m_switch.mat"), "after_stimulus_switch_m", "-v7.3")
save(fullfile(save_path, "after_stimulus_data_np2_switch.mat"), "after_stimulus_switch_y", "-v7.3")



