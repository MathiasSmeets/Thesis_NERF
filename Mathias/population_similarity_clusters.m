%clear; clc; close all;

if ispc
    volume_base = '\\nerffs13';
    volume_base2 = '\\nerffs17';
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/';
    volume_base2 = '/mnt/takeokalab/';
end

path_to_code = "takeokalabwip2023/Mathias/switch_data/clusters";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

clusters_after = load(fullfile(volume_base2, path_to_code,"cluster_matrices_after_m.mat"));
clusters_after = clusters_after.all_cluster_matrices;
clusters_before = load(fullfile(volume_base2, path_to_code,"cluster_matrices_before_m.mat"));
clusters_before = clusters_before.all_cluster_matrices;
clusters_between = load(fullfile(volume_base2, path_to_code,"cluster_matrices_between_m.mat"));
clusters_between = clusters_between.all_cluster_matrices;
% stimulus_data_m = load(fullfile(volume_base2, path_to_code,"after_stimulus_data_m_switch.mat"));
% stimulus_data_m = stimulus_data_m.after_stimulus_switch_m;
% stimulus_data_y = load(fullfile(volume_base2, path_to_code,"after_stimulus_data_y_switch.mat"));
% stimulus_data_y = stimulus_data_y.after_stimulus_switch_y;

folder = fileparts(which("clusters_cpd.m"));
addpath(genpath(folder))


%% calculate similarities (over intervals_to_combine)
nb_combinations_before = size(cur_cluster_before,2)*size(cur_cluster_between,2);
nb_combinations_after = size(cur_cluster_after,2)*size(cur_cluster_between,2);
population_similarities_before = zeros(nb_combinations_before, size(clusters_before,1));
population_similarities_after = zeros(nb_combinations_after, size(clusters_after,1));
% loop over each mouse
for i = 1:numel(clusters_after)
    
    cur_cluster_before = clusters_before{i};
    cur_cluster_after = clusters_after{i};
    cur_cluster_between = clusters_between{i};

    % loop over all columns of both matrices
    counter = 1;
    for k = 1:size(cur_cluster_before,2)
        for l = 1:size(cur_cluster_between,2)
            cur_x = cur_cluster_before(:,k);
            cur_y = cur_cluster_between(:,l);
            if ~all(cur_x==0) && ~all(cur_y==0)
                population_similarities_before(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
            end
            counter = counter + 1;
        end
    end

    counter = 1;
    for k = 1:size(cur_cluster_after,2)
        for l = 1:size(cur_cluster_between,2)
            cur_x = cur_cluster_after(:,k);
            cur_y = cur_cluster_between(:,l);
            if ~all(cur_x==0) && ~all(cur_y==0)
                population_similarities_after(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
            end
            counter = counter + 1;
        end
    end
end

average_similarities_start = mean(population_similarities_before)';
average_similarities_end = mean(population_similarities_after)';

figure;hold on;
for i = 1:size(clusters_after,1)
    subplot(6,4,i)
    plot(nonzeros(population_similarities(i,:)))
    title("Population Similarities Horridge")
    xlabel("Intervals")
    ylabel("Cosine Similarity")
    ylim([0 1])
end
%save("X:\Mathias\switch_data\population_similarities\original_large_dataset_y.mat", "population_similarities")

%% create plots

% population_similarities_m = load("X:\Mathias\switch_data\population_similarities\horridge_m.mat"); population_similarities_m = population_similarities_m.population_similarities;
% population_similarities_y = load("X:\Mathias\switch_data\population_similarities\horridge_y.mat"); population_similarities_y = population_similarities_y.population_similarities;
% % first intervals
% figure;boxplot([population_similarities_m(:,1) population_similarities_y(:,1)]);title("Horridge; Learner vs Control; First 5 intervals");
% % last intervals
% for i = 1:11
% population_similarities_m_figure(i) = population_similarities_m(i,find(population_similarities_m(i,:),1,'last'));
% population_similarities_y_figure(i) = population_similarities_y(i,find(population_similarities_y(i,:),1,'last'));
% end
% figure;boxplot([population_similarities_m_figure' population_similarities_y_figure']);title("Horridge; Learner vs Control; Last 5 intervals")




