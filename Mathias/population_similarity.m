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

%% calculate mean activity after stimulus (12ms-70ms)

all_mean_activity = cell(size(stimulus_data_m));
% loop over all mice
for i = 1:size(stimulus_data_m,1)
    % loop over all intervals
    for j = 1:size(stimulus_data_m,2)
        %all_mean_activity{i,j} = sum(stimulus_data_m{i,j}(:,12:end),2);
        all_mean_activity{i,j} = mean(stimulus_data_m{i,j}(:,12:end),2);
    end
end

%% calculate similarities (over intervals_to_combine)

population_similarities = zeros(size(stimulus_data_m,1), ceil(size(stimulus_data_m,2) / intervals_to_combine));
% loop over each mouse
for i = 1:size(stimulus_data_m,1)
    % loop over intervals
    counter = 0;
    for j = 1:intervals_to_combine:size(stimulus_data_m,2)-intervals_to_combine+1
        if size(stimulus_data_m,2) > j+intervals_to_combine-1
            if ~isempty(stimulus_data_m{i,j+intervals_to_combine-1})
                counter = counter + 1;
                cur_numerator = 0;
                cur_x_denum = 0;
                cur_y_denum = 0;
                for k = 0:intervals_to_combine-2
                    for l = k+1:intervals_to_combine-1
                        cur_x = all_mean_activity{i,j+k};
                        cur_y = all_mean_activity{i,j+l};
                        cur_numerator = cur_numerator + dot(cur_x,cur_y);
                        cur_x_denum = cur_x_denum + dot(cur_x, cur_x);
                        cur_y_denum = cur_y_denum + dot(cur_y, cur_y);
                    end
                end
                population_similarities(i,counter) = cur_numerator / sqrt(cur_x_denum*cur_y_denum);
            end
        end
    end
end
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






