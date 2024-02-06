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

%% variables initialization

interval_size = size(stimulus_data_m{1,1},2);
wanted_bin_size = 10;
interval_step = 10;
neuron_counter = 1;
indices = ceil((1:interval_size)/wanted_bin_size);
total_nb_assemblies = cell(size(stimulus_data_m));
total_nb_neurons = cell(size(stimulus_data_m));
total_assemblies = cell(size(stimulus_data_m));
total_activity = cell(size(stimulus_data_m));
total_neurons_of_interest = cell(size(stimulus_data_m));
total_data = cell(size(stimulus_data_m));

%% create 3D tensor with (dimension neurons of interest) x (interval_size/wanted_bin_size) x (interval_step)

for k = 1:size(stimulus_data_m, 1)
    index_counter = 1;
    cur_neurons_of_interest = get_neurons_of_interest(stimulus_data_m{k,1}, neurons_of_interest_m, inhibited_neurons_m, neuron_counter);
    cur_tensor_10ms_neurons_oi = zeros(length(cur_neurons_of_interest), size(stimulus_data_m{k,1},2)/wanted_bin_size);
    for i = 1:interval_step:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{k,i})
            for j = i:i+wanted_bin_size-1
                if ~isempty(stimulus_data_m{k,j})

                    cur_mouse = stimulus_data_m{k,j}(cur_neurons_of_interest,:);
                    cur_mouse_fs_adjusted = zeros(size(cur_mouse,1),size(cur_mouse,2)/wanted_bin_size);
                    for jj = 1:size(cur_mouse,1)
                        cur_mouse_fs_adjusted(jj,:) = accumarray(indices',cur_mouse(jj,:)',[],@sum)';
                    end

                    if j == i
                        cur_tensor_10ms_neurons_oi(:,:,1) = cur_mouse_fs_adjusted;
                    else
                        cur_tensor_10ms_neurons_oi(:,:,j-i+1) = cur_mouse_fs_adjusted;
                    end
                end
            end
            % remove neurons that are not active
            for jj = size(cur_tensor_10ms_neurons_oi,1):-1:1
                if all(cur_tensor_10ms_neurons_oi(jj,:,:)==cur_tensor_10ms_neurons_oi(jj,1,1))
                    cur_tensor_10ms_neurons_oi(jj,:,:) = [];
                    cur_neurons_of_interest(jj) = [];
                end
            end

            [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = cpd_assembly_detection(cur_tensor_10ms_neurons_oi, ica_assemblies{k,index_counter});
            if predicted_nbr_assemblies ~= 0
                total_neurons_of_interest{k,index_counter} = cur_neurons_of_interest;
                total_nb_assemblies{k,index_counter} = predicted_nbr_assemblies;
                total_nb_neurons{k,index_counter} = predicted_nbr_neurons;
                total_assemblies{k,index_counter} = assemblies;
                total_activity{k,index_counter} = activity;
                total_data{k,index_counter} = cur_tensor_10ms_neurons_oi;
            end
            index_counter = index_counter + 1;
        end
    end
    neuron_counter = neuron_counter + size(stimulus_data_m{k,1},1);
    disp(k)
end









