clear; clc; close all;

%% get data

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

sos_results_m = load(fullfile(volume_base2, path_to_code, "sos_results_m.mat"));
sos_results_m = sos_results_m.sos_results_m;

folder = fileparts(which("clusters_ica.m"));
addpath(genpath(folder))

%% variables initialization

interval_size = size(stimulus_data_m{1,1},2);
wanted_bin_size = 10;
interval_step = 20;
create_plots = false;
neuron_counter = 1;
index_counter = 1;
indices = ceil((1:interval_size)/wanted_bin_size);
    
total_nb_assemblies = cell(size(stimulus_data_m));
total_nb_neurons = cell(size(stimulus_data_m));
total_assemblies = cell(size(stimulus_data_m));
total_activity = cell(size(stimulus_data_m));
total_neurons_of_interest = cell(size(stimulus_data_m));
total_data = cell(size(stimulus_data_m));

%% calculate assembly for each interval

% loop over each mouse
for k = 1:size(stimulus_data_m,1)
    % loop over each interval of this mouse
    cur_neurons_of_interest = get_neurons_of_interest(stimulus_data_m{k,1}, neurons_of_interest_m, inhibited_neurons_m, neuron_counter);
    cur_neurons_of_interest = 1:size(stimulus_data_m{k,1});
    index_counter = 1;
    for i = 1:interval_step:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{k,i})
            cur_total_mouse = [];
            for ii = i:i+interval_step-1
                if ~isempty(stimulus_data_m{k,ii})
                    cur_mouse = stimulus_data_m{k,ii}(cur_neurons_of_interest,:);

                    % transform to 10ms bins
                    cur_mouse_fs_adjusted = zeros(size(cur_mouse,1),size(cur_mouse,2)/wanted_bin_size);
                    for j = 1:size(cur_mouse,1)
                        cur_mouse_fs_adjusted(j,:) = accumarray(indices',cur_mouse(j,:)',[],@sum)';
                    end
                    cur_total_mouse = [cur_total_mouse, cur_mouse_fs_adjusted];
                end
            end

            % calculate zscores
            cur_std = sos_results_m{k,2};
            cur_std(cur_std==0) = 0.01;
            cur_total_mouse_zscore = (cur_total_mouse - sos_results_m{k,1}(cur_neurons_of_interest,:)) ./ cur_std(cur_neurons_of_interest,:);
            
            % remove neurons that are not active
            for jj = size(cur_total_mouse_zscore,1):-1:1
                if all(cur_total_mouse_zscore(jj,:)==cur_total_mouse_zscore(jj,1))
                    cur_total_mouse_zscore(jj,:) = [];
                    cur_neurons_of_interest(jj) = [];
                end
            end
            if k == 7 && index_counter == 9
                disp("break")
            end
            [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = ica_assembly_detection(cur_total_mouse_zscore', create_plots);
            if predicted_nbr_assemblies ~= 0
                total_neurons_of_interest{k,index_counter} = cur_neurons_of_interest;
                total_nb_assemblies{k,index_counter} = predicted_nbr_assemblies;
                total_nb_neurons{k,index_counter} = predicted_nbr_neurons;
                total_assemblies{k,index_counter} = assemblies;
                total_activity{k,index_counter} = activity;
                total_data{k,index_counter} = cur_total_mouse_zscore;
            end
            index_counter = index_counter + 1;
        end
    end
    neuron_counter = neuron_counter + size(stimulus_data_m{k,1},1);
    disp(k)
end


% activity = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\activity.mat"); activity = activity.total_activity;
% nb_neurons = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\nb_neurons.mat"); nb_neurons = nb_neurons.total_nb_neurons;
% assemblies = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\assemblies.mat"); assemblies = assemblies.total_assemblies;
% nb_assemblies = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\nb_assemblies.mat"); nb_assemblies = nb_assemblies.total_nb_assemblies;

