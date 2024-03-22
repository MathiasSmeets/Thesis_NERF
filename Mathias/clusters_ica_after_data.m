clear; clc; close all;
% use this code for waiting
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

path_to_data = "takeokalabwip2023/Mathias/switch_data/data_after_stimulus";
path_to_neurons_of_interest = "takeokalabwip2023/Mathias/switch_data/neurons_of_interest";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

stimulus_data_m = load(fullfile(volume_base2, path_to_data,"waiting_data_np2.mat"));
stimulus_data_m = stimulus_data_m.waiting_data;
%stimulus_data_y = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_y.mat"));
%stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

neurons_of_interest_m = load(fullfile(volume_base2, path_to_neurons_of_interest, "neurons_of_interest_horridge_np2.mat"));
neurons_of_interest_m = neurons_of_interest_m.output_m;
inhibited_neurons_m = load(fullfile(volume_base2, path_to_neurons_of_interest, "inhibited_horridge_np2.mat"));
inhibited_neurons_m = inhibited_neurons_m.inhibited_m;

sos_results_m = load(fullfile(volume_base2, "takeokalabwip2023","Mathias", "switch_data","sos_data", "sos_results_np2.mat"));
sos_results_m = sos_results_m.sos_results_m;

folder = fileparts(which("clusters_ica.m"));
addpath(genpath(folder))

%% variables initialization

wanted_bin_size = 15;
interval_step = 1;
interval_size = 30*70;
create_plots = false;
neuron_counter = 1;
index_counter = 1;
indices = ceil((1:interval_size)/wanted_bin_size);

largest_data_mouse = 0;
for i = 1:size(stimulus_data_m,1)
    cur_data = max(largest_data_mouse,size(stimulus_data_m{i},1));
end
total_nb_assemblies = cell(size(stimulus_data_m,1),ceil(largest_data_mouse/interval_size));
total_nb_neurons = cell(size(stimulus_data_m,1),ceil(largest_data_mouse/interval_size));
total_assemblies = cell(size(stimulus_data_m,1),ceil(largest_data_mouse/interval_size));
total_activity = cell(size(stimulus_data_m,1),ceil(largest_data_mouse/interval_size));
total_neurons_of_interest = cell(size(stimulus_data_m,1),ceil(largest_data_mouse/interval_size));
total_data = cell(size(stimulus_data_m,1),ceil(largest_data_mouse/interval_size));

%% calculate assembly for each interval

% loop over each mouse
for k = 1:size(stimulus_data_m,1)
    cur_after_data = stimulus_data_m{k};
    cur_after_data = cell2mat(cur_after_data);
    %cur_neurons_of_interest = get_neurons_of_interest(cur_after_data, neurons_of_interest_m, inhibited_neurons_m, neuron_counter);
    cur_neurons_of_interest = 1:size(cur_after_data,1);
    reduced_cur_after_data = cur_after_data(cur_neurons_of_interest,:);
    index_counter = 1;
    for i = 1:interval_size:size(reduced_cur_after_data,2)-interval_size+1
        cur_data = reduced_cur_after_data(:,i:i+interval_size-1);
        % transform to larger bins
        cur_mouse_fs_adjusted = zeros(size(cur_data,1),ceil(interval_size/wanted_bin_size));
        for j = 1:size(cur_data,1)
            cur_mouse_fs_adjusted(j,:) = accumarray(indices',cur_data(j,:)',[],@sum)';
        end
        % calculate zscores
        cur_std = sos_results_m{k,2};
        cur_std(cur_std==0) = 0.01;
        %cur_neurons_of_interest = get_neurons_of_interest(cur_after_data, neurons_of_interest_m, inhibited_neurons_m, neuron_counter);
        cur_neurons_of_interest = 1:size(cur_after_data,1);
        cur_data_zscore = (cur_mouse_fs_adjusted - sos_results_m{k,1}(cur_neurons_of_interest,:)) ./ cur_std(cur_neurons_of_interest,:);

        % remove neurons that are not active
        for jj = size(cur_data_zscore,1):-1:1
            if all(cur_data_zscore(jj,:)==cur_data_zscore(jj,1))
                cur_data_zscore(jj,:) = [];
                cur_neurons_of_interest(jj) = [];
            end
        end

        [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity,~] = ica_assembly_detection(cur_data_zscore', create_plots);
        if predicted_nbr_assemblies ~= 0
            total_neurons_of_interest{k,index_counter} = cur_neurons_of_interest;
            total_nb_assemblies{k,index_counter} = predicted_nbr_assemblies;
            total_nb_neurons{k,index_counter} = predicted_nbr_neurons;
            total_assemblies{k,index_counter} = assemblies;
            total_activity{k,index_counter} = activity;
            total_data{k,index_counter} = cur_data_zscore;
        end
        index_counter = index_counter + 1;
    end
    neuron_counter = neuron_counter + size(cur_after_data,1);
    disp(k)
end

savepath = "/scratch/mathiass-takeokalab/01/";
save(fullfile(savepath, "neurons_of_interest_after_np2.mat"), "total_neurons_of_interest", "-v7.3")
save(fullfile(savepath, "nb_assemblies_after_np2.mat"), "total_nb_assemblies", "-v7.3")
save(fullfile(savepath, "nb_neurons_after_np2.mat"), "total_nb_neurons", "-v7.3")
save(fullfile(savepath, "assemblies_after_np2.mat"), "total_assemblies", "-v7.3")
save(fullfile(savepath, "activity_after_np2.mat"), "total_activity", "-v7.3")
save(fullfile(savepath, "data_after_np2.mat"), "total_data", "-v7.3")


% activity = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\activity.mat"); activity = activity.total_activity;
% nb_neurons = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\nb_neurons.mat"); nb_neurons = nb_neurons.total_nb_neurons;
% assemblies = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\assemblies.mat"); assemblies = assemblies.total_assemblies;
% nb_assemblies = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\nb_assemblies.mat"); nb_assemblies = nb_assemblies.total_nb_assemblies;

%% create figure on how this evolves

intervals_to_combine = 3;
population_similarities_start = zeros(2+1,size(total_assemblies,1));
population_similarities_end = zeros(2+1,size(total_assemblies,1));
all_cluster_matrices = cell(1,size(stimulus_data_m,1));
% loop over all mice
for i = 1:size(stimulus_data_m,1)
    % get last interval
    cur_after_data = stimulus_data_m{i};
    cur_after_data = cell2mat(cur_after_data);
    last_interval_index = floor(size(cur_after_data,2)/interval_size);

    % create clustermatrix
    cluster_matrix = zeros(size(cur_after_data,1), last_interval_index);
    for k = 1:size(total_assemblies,2)
        for l = 1:length(total_assemblies{i,k})
            cluster_matrix(total_neurons_of_interest{i,k}(total_assemblies{i,k}{l}),k) = 1;
        end
    end
    figure
    heatmap(cluster_matrix,'CellLabelColor','none')
    xlabel("Time Bins (" + interval_size + " together)")
    ylabel("Neurons")
    title("Neurons in a cluster")
    all_cluster_matrices{i} = cluster_matrix;

    % row_indices = [3,14,19,30];
    % for column_index = 1:size(cluster_matrix, 2)
    %     for cur_row = 1:numel(row_indices)
    %         if cluster_matrix(row_indices(cur_row), column_index) == 1
    %             % If all are 1, transform those elements to 2
    %             cluster_matrix(row_indices(cur_row), column_index) = 2;
    %         end
    %     end
    % end
    

    % calculate population similarity for this matrix
    % similarity at the start
    % counter = 0;
    % for k = 1:intervals_to_combine-1
    %     for l = k+1:intervals_to_combine
    %         if ~all(cluster_matrix(:,k)==0) && ~all(cluster_matrix(:,l)==0)
    %             counter = counter + 1;
    %             cur_x = cluster_matrix(:,k);
    %             cur_y = cluster_matrix(:,l);
    %             population_similarities_start(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
    %         end
    %     end
    % end
    % % similarity at the end
    % counter = 0;
    % last_interval = size(cluster_matrix,2);
    % for k = last_interval:-1:last_interval-intervals_to_combine+2
    %     for l = k-1:-1:last_interval-intervals_to_combine+1
    %         if ~all(cluster_matrix(:,k)==0) && ~all(cluster_matrix(:,l)==0)
    %             counter = counter + 1;
    %             cur_x = cluster_matrix(:,k);
    %             cur_y = cluster_matrix(:,l);
    %             population_similarities_end(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
    %         end
    %     end
    % end
    % % similarity between start and end
    % counter = 0;
    % last_interval = size(cluster_matrix,2);
    % for k = 1:intervals_to_combine
    %     for l = last_interval:-1:last_interval-intervals_to_combine+1
    %         if ~all(cluster_matrix(:,k)==0) && ~all(cluster_matrix(:,l)==0)
    %             counter = counter + 1;
    %             cur_x = cluster_matrix(:,k);
    %             cur_y = cluster_matrix(:,l);
    %             population_similarities_start_end(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
    %         end
    %     end
    % end
end
save(fullfile(savepath, "cluster_matrices_after_np2.mat"), "all_cluster_matrices", "-v7.3")

