clear; clc; close all;

if ispc
    volume_base = '\\nerffs13';
    volume_base2 = '\\nerffs17';
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/';
    volume_base2 = '/mnt/takeokalab/';
end

path_to_code = "takeokalabwip2023/Mathias/data";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

stimulus_data_m = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_m.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_y = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_y.mat"));
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

neurons_of_interest_m = load(fullfile(volume_base2, path_to_code, "output_m.mat"));
neurons_of_interest_m = neurons_of_interest_m.output_m;
inhibited_neurons_m = load(fullfile(volume_base2, path_to_code, "inhibited_m.mat"));
inhibited_neurons_m = inhibited_neurons_m.inhibited_m;

folder = fileparts(which("clusters_cpd.m"));
addpath(genpath(folder))

%% load in ica data

ica_activity = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\activity.mat"); ica_activity = ica_activity.total_activity;
ica_nb_neurons = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\nb_neurons.mat"); ica_nb_neurons = ica_nb_neurons.total_nb_neurons;
ica_assemblies = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\assemblies.mat"); ica_assemblies = ica_assemblies.total_assemblies;
ica_nb_assemblies = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\nb_assemblies.mat"); ica_nb_assemblies = ica_nb_assemblies.total_nb_assemblies;
ica_data = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\data.mat"); ica_data = ica_data.total_data;
ica_neurons_of_interest = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\neurons_of_interest.mat"); ica_neurons_of_interest = ica_neurons_of_interest.total_neurons_of_interest;

%% initialize data
interval_step = 10; % depends on how data was created

%% get cluster shape for each assembly

templates = cell(size(ica_activity));
% loop over each mouse
for k = size(templates,1)
    % loop over each interval
    for i = size(templates,2)
        cur_assembly = ica_assemblies{k,i};
        if isempty(cur_assembly)
            templates{k,i} = [];
        elseif size(cur_assembly,1) == 1
            ...
        else
            ...
        end
    end
end




