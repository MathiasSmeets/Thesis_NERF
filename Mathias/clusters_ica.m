clear; clc; close all;

%% get data

%path_to_code = "\\nerffs13\takeokalabwip2020\Mathias\data\";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

stimulus_data_m = load(path_to_code + "data_after_stimulus_m.mat");
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_y = load(path_to_code + "data_after_stimulus_y.mat");
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

%% put this in matrix format
mouse_1_m = stimulus_data_m{1,1};
for i = 1:size(stimulus_data_m,2)
    if ~isempty(stimulus_data_m{1,i})
        mouse_1_m = [mouse_1_m, stimulus_data_m{1,i}];
    end
end

%% arnoud function
[predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = ica_assembly_detection(mouse_1_m,1);