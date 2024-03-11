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

path_to_code = "takeokalabwip2023/Mathias/data";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

stimulus_data_m = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_y.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_y;
% stimulus_data_m = load(fullfile(volume_base2, path_to_code,"after_stimulus_data_m_switch.mat"));
% stimulus_data_m = stimulus_data_m.after_stimulus_switch_m;
% stimulus_data_y = load(fullfile(volume_base2, path_to_code,"after_stimulus_data_y_switch.mat"));
% stimulus_data_y = stimulus_data_y.after_stimulus_switch_y;

folder = fileparts(which("clusters_cpd.m"));
addpath(genpath(folder))

%% hyperparameters

intervals_to_combine = 5;


%% calculate all possible connections
nb_combinations = 0;
for i = 1:intervals_to_combine-1
    nb_combinations = nb_combinations+i;
end

%% calculate mean activity after stimulus (12ms-70ms)

all_mean_activity = cell(size(stimulus_data_m));
% loop over all mice
for i = 1:size(stimulus_data_m,1)
    % loop over all intervals
    for j = 1:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{i,j})
            all_mean_activity{i,j} = mean(stimulus_data_m{i,j}(:,12:end),2);
        end
    end
end

%% calculate similarities (over intervals_to_combine)

population_similarities_start = zeros(nb_combinations, size(stimulus_data_m,1));
population_similarities_end = zeros(nb_combinations, size(stimulus_data_m,1));
% loop over each mouse
for i = 1:size(stimulus_data_m,1)
    % calculate last interval
    j = size(stimulus_data_m,2);
    while isempty(stimulus_data_m{i,j})
        j = j-1;
    end
    last_interval = j;

    % loop over first intervals_to_combine intervals   
    counter = 0;
    for k = 1:intervals_to_combine-1
        for l = k+1:intervals_to_combine
            if ~isempty(stimulus_data_m{i,k}) && ~isempty(stimulus_data_m{i,l})
                counter = counter + 1;
                cur_x = all_mean_activity{i,k};
                cur_y = all_mean_activity{i,l};
                population_similarities_start(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
            end
        end
    end

    % loop over last intervals_to_combine intervals
    counter = 0;
    for k = last_interval:-1:last_interval-intervals_to_combine+2
        for l = k-1:-1:last_interval-intervals_to_combine+1
            if ~isempty(stimulus_data_m{i,k}) && ~isempty(stimulus_data_m{i,l})
                counter = counter + 1;
                cur_x = all_mean_activity{i,k};
                cur_y = all_mean_activity{i,l};
                population_similarities_end(counter,i) = dot(cur_x,cur_y) / sqrt(dot(cur_x, cur_x) * dot(cur_y, cur_y));
            end
        end
    end
end

average_similarities_start = mean(population_similarities_start)';
average_similarities_end = mean(population_similarities_end)';

figure;hold on;
for i = 1:size(stimulus_data_m,1)
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






