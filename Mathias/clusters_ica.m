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

folder = fileparts(which("clusters_ica.m"));
addpath(genpath(folder))

%% variables initialization

interval_size = size(stimulus_data_m{1,1},2);
wanted_bin_size = 10;
create_plots = false;
neuron_counter = 1;
indices = ceil((1:interval_size)/10);
total_nb_assemblies = cell(size(stimulus_data_m));
total_nb_neurons = cell(size(stimulus_data_m));
total_assemblies = cell(size(stimulus_data_m));
total_activity = cell(size(stimulus_data_m));

%% calculate assembly for each interval

% loop over each mouse
for k = 1:size(stimulus_data_m,1)
    % loop over each interval of this mouse
    for i = 1:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{k,i})
            
            % get neurons of interest
            cur_neurons_of_interest = [];
            for ii = 1:size(neurons_of_interest_m,2)-1 % do not use this for inhibited neurons
                for jj = 1:size(neurons_of_interest_m{1,ii}{1,1},1)
                    if neurons_of_interest_m{1,ii}{1,1}{jj,1} >= neuron_counter && neurons_of_interest_m{1,ii}{1,1}{jj,1} < neuron_counter + size(stimulus_data_m{k,1},1)
                        cur_neurons_of_interest = [cur_neurons_of_interest, neurons_of_interest_m{1,ii}{1,1}{jj,1}];
                    end
                end
            end
            inhibited_neurons_of_interest = inhibited_neurons_m(inhibited_neurons_m >= neuron_counter & inhibited_neurons_m < neuron_counter + size(stimulus_data_m{k,1},1));
            cur_neurons_of_interest = [cur_neurons_of_interest, inhibited_neurons_of_interest];
            cur_neurons_of_interest = unique(cur_neurons_of_interest);
            cur_neurons_of_interest = cur_neurons_of_interest - neuron_counter + 1;

            cur_mouse = stimulus_data_m{k,i}(cur_neurons_of_interest,:);

            % transform to 10ms bins
            cur_mouse_fs_adjusted = zeros(size(cur_mouse,1),size(cur_mouse,2)/wanted_bin_size);
            for j = 1:size(cur_mouse,1)
                cur_mouse_fs_adjusted(j,:) = accumarray(indices',cur_mouse(j,:)',[],@sum)';
            end

            % calculate zscores
            cur_mean_mouse = mean(cur_mouse_fs_adjusted,2);
            cur_std_mouse = std(cur_mouse_fs_adjusted,[],2);
            cur_std_mouse(cur_std_mouse == 0) = 0.1;
            cur_mouse_zscore = (cur_mouse_fs_adjusted - cur_mean_mouse) ./ cur_std_mouse;

            [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = ica_assembly_detection(cur_mouse_zscore,create_plots);
            if predicted_nbr_assemblies ~= 0
                total_nb_assemblies{k,i} = predicted_nbr_assemblies;
                total_nb_neurons{k,i} = predicted_nbr_neurons;
                total_assemblies{k,i} = assemblies;
                total_activity{k,i} = activity;
            end
        end
    end
    neuron_counter = neuron_counter + size(stimulus_data_m{k,1},1);
end




