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

path_to_data = "takeokalabwip2023/Mathias/switch_data/data_after_stimulus";
path_to_neurons_of_interest = "takeokalabwip2023/Mathias/switch_data/neurons_of_interest";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

stimulus_data_m = load(fullfile(volume_base2, path_to_data,"after_stimulus_data_m_horridge.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
%stimulus_data_y = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_y.mat"));
%stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

neurons_of_interest_m = load(fullfile(volume_base2, path_to_neurons_of_interest, "neurons_of_interest_horridge_m.mat"));
neurons_of_interest_m = neurons_of_interest_m.output_m;
inhibited_neurons_m = load(fullfile(volume_base2, path_to_neurons_of_interest, "inhibited_horridge_m.mat"));
inhibited_neurons_m = inhibited_neurons_m.inhibited_m;

sos_results_m = load(fullfile(volume_base2, "takeokalabwip2023","Mathias", "switch_data","sos_data", "sos_results_m.mat"));
sos_results_m = sos_results_m.sos_results_m;

folder = fileparts(which("clusters_ica.m"));
addpath(genpath(folder))

%% variables initialization

interval_size = 60;
wanted_bin_size = 15;
interval_step = 30;
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
total_vector = cell(size(stimulus_data_m));

%% calculate assembly for each interval

% loop over each mouse
for k = 1:size(stimulus_data_m,1)
    % loop over each interval of this mouse
    cur_neurons_of_interest = get_neurons_of_interest(stimulus_data_m{k,1}, neurons_of_interest_m, inhibited_neurons_m, neuron_counter);
    %cur_neurons_of_interest = 1:size(stimulus_data_m{k,1});
    index_counter = 1;
    for i = 1:interval_step:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{k,i})
            cur_total_mouse = [];
            for ii = i:i+interval_step-1
                if ii <= size(stimulus_data_m,2)
                    if 70 <= size(stimulus_data_m{k,ii},2)
                        if ~isempty(stimulus_data_m{k,ii})
                            cur_mouse = stimulus_data_m{k,ii}(cur_neurons_of_interest,11:end);

                            % transform to 15ms bins
                            cur_mouse_fs_adjusted = zeros(size(cur_mouse,1),ceil(size(cur_mouse,2)/wanted_bin_size));
                            for j = 1:size(cur_mouse,1)
                                cur_mouse_fs_adjusted(j,:) = accumarray(indices',cur_mouse(j,:)',[],@sum)';
                            end
                            cur_total_mouse = [cur_total_mouse, cur_mouse_fs_adjusted];
                        end
                    end
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
            % up until now: simply create input for ica method
            % now we perform the ica method
            [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity,vector] = ica_assembly_detection(cur_total_mouse_zscore', create_plots);
            if predicted_nbr_assemblies ~= 0
                total_neurons_of_interest{k,index_counter} = cur_neurons_of_interest;
                total_nb_assemblies{k,index_counter} = predicted_nbr_assemblies;
                total_nb_neurons{k,index_counter} = predicted_nbr_neurons;
                total_assemblies{k,index_counter} = assemblies;
                total_activity{k,index_counter} = activity;
                total_data{k,index_counter} = cur_total_mouse_zscore;
                total_vector{k, index_counter} = vector;
            end
            index_counter = index_counter + 1;
        end
    end
    neuron_counter = neuron_counter + size(stimulus_data_m{k,1},1);
    disp(k)
end

% in the second mouse of the control group, there is no cluster detected

savepath = "X:\Mathias\switch_data\clusters";
%savepath = "/scratch/mathiass-takeokalab/01/";
save(fullfile(savepath, "neurons_of_interest_horridge_m.mat"), "total_neurons_of_interest", "-v7.3")
save(fullfile(savepath, "nb_assemblies_horridge_m.mat"), "total_nb_assemblies", "-v7.3")
save(fullfile(savepath, "nb_neurons_horridge_m.mat"), "total_nb_neurons", "-v7.3")
save(fullfile(savepath, "assemblies_horridge_m.mat"), "total_assemblies", "-v7.3")
save(fullfile(savepath, "activity_horridge_m.mat"), "total_activity", "-v7.3")
save(fullfile(savepath, "data_horridge_m.mat"), "total_data", "-v7.3")
save(fullfile(savepath, "ica_vector_horridge_m.mat"), "total_vector", "-v7.3")


% activity = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\activity.mat"); activity = activity.total_activity;
% nb_neurons = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\nb_neurons.mat"); nb_neurons = nb_neurons.total_nb_neurons;
% assemblies = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\assemblies.mat"); assemblies = assemblies.total_assemblies;
% nb_assemblies = load("X:\Mathias\cluster_output\bin_10ms_neurons_oi\nb_assemblies.mat"); nb_assemblies = nb_assemblies.total_nb_assemblies;

%% create figure on how this evolves
% ignore this part, old

% intervals_to_combine = 3;
% population_similarities_start = zeros(2+1,size(total_assemblies,1));
% population_similarities_end = zeros(2+1,size(total_assemblies,1));
% all_cluster_matrices = cell(1,size(stimulus_data_m,1));
% % loop over all mice
% for i = 1:3
%     % get last interval
%     j = size(total_data,2);
%     while isempty(total_data{i,j}) && j > 1
%         j = j-1;
%     end
%     last_interval_index = j;
% 
%     % create clustermatrix
%     cluster_matrix = zeros(size(stimulus_data_m{i,1},1), last_interval_index);
%     for k = 1:last_interval_index
%         for l = 1:length(total_assemblies{i,k})
%             cluster_matrix(total_neurons_of_interest{i,k}(total_assemblies{i,k}{l}),k) = 1;
%         end
%     end
%     figure
%     imagesc(cluster_matrix)
%     % heatmap(cluster_matrix,'CellLabelColor','none')
%     % xlabel("Intervals (" + interval_step + " together)")
%     % ylabel("Neurons")
%     % title("Neurons in a cluster")
%     all_cluster_matrices{i} = cluster_matrix;
%     % row_indices = [3,14,19,30];
%     % for column_index = 1:size(cluster_matrix, 2)
%     %     for cur_row = 1:numel(row_indices)
%     %         if cluster_matrix(row_indices(cur_row), column_index) == 1
%     %             % If all are 1, transform those elements to 2
%     %             cluster_matrix(row_indices(cur_row), column_index) = 2;
%     %         end
%     %     end
%     % end
% 
% 
%     % calculate population similarity for this matrix
%     % similarity at the start
%     % % % % % counter = 0;
%     % % % % % for k = 1:intervals_to_combine-1
%     % % % % %     for l = k+1:intervals_to_combine
%     % % % % %         if ~all(cluster_matrix(:,k)==0) && ~all(cluster_matrix(:,l)==0)
%     % % % % %             counter = counter + 1;
%     % % % % %             cur_x = cluster_matrix(:,k);
%     % % % % %             cur_y = cluster_matrix(:,l);
%     % % % % %             population_similarities_start(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
%     % % % % %         end
%     % % % % %     end
%     % % % % % end
%     % % % % % % similarity at the end
%     % % % % % counter = 0;
%     % % % % % last_interval = size(cluster_matrix,2);
%     % % % % % for k = last_interval:-1:last_interval-intervals_to_combine+2
%     % % % % %     for l = k-1:-1:last_interval-intervals_to_combine+1
%     % % % % %         if ~all(cluster_matrix(:,k)==0) && ~all(cluster_matrix(:,l)==0)
%     % % % % %             counter = counter + 1;
%     % % % % %             cur_x = cluster_matrix(:,k);
%     % % % % %             cur_y = cluster_matrix(:,l);
%     % % % % %             population_similarities_end(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
%     % % % % %         end
%     % % % % %     end
%     % % % % % end
%     % % % % % % similarity between start and end
%     % % % % % counter = 0;
%     % % % % % last_interval = size(cluster_matrix,2);
%     % % % % % for k = 1:intervals_to_combine
%     % % % % %     for l = last_interval:-1:last_interval-intervals_to_combine+1
%     % % % % %         if ~all(cluster_matrix(:,k)==0) && ~all(cluster_matrix(:,l)==0)
%     % % % % %             counter = counter + 1;
%     % % % % %             cur_x = cluster_matrix(:,k);
%     % % % % %             cur_y = cluster_matrix(:,l);
%     % % % % %             population_similarities_start_end(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
%     % % % % %         end
%     % % % % %     end
%     % % % % % end
% end
% save(fullfile(savepath, "cluster_matrices_between_m.mat"), "all_cluster_matrices", "-v7.3")
