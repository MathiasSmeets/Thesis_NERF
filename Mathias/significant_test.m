clear;clc;close all;
if ispc
    volume_base = '\\nerffs13';
    volume_base2 = '\\nerffs17';
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/';
    volume_base2 = '/mnt/takeokalab/';
end

path_to_correlations = "takeokalabwip2023/Mathias/switch_data/correlations";

raw_before = load(fullfile(volume_base2, path_to_correlations, "raw_before_m.mat")); raw_before = raw_before.cur_correlation_before;
raw_between = load(fullfile(volume_base2, path_to_correlations, "raw_between_m.mat")); raw_between = raw_between.cur_correlation_between;
raw_after = load(fullfile(volume_base2, path_to_correlations, "raw_after_m.mat")); raw_after = raw_after.cur_correlation_after;
raw_horridge = load(fullfile(volume_base2, path_to_correlations, "raw_horridge_m.mat")); raw_horridge = raw_horridge.cur_correlation_horridge;

correlation_distribution_before = load(fullfile(volume_base2, path_to_correlations, "correlation_distribution_before_m.mat")); correlation_distribution_before = correlation_distribution_before.correlation_distribution_before;
correlation_distribution_between = load(fullfile(volume_base2, path_to_correlations, "correlation_distribution_between_m.mat")); correlation_distribution_between = correlation_distribution_between.correlation_distribution_between;
correlation_distribution_after = load(fullfile(volume_base2, path_to_correlations, "correlation_distribution_after_m.mat")); correlation_distribution_after = correlation_distribution_after.correlation_distribution_after;
correlation_distribution_horridge = load(fullfile(volume_base2, path_to_correlations, "correlation_distribution_horridge_m.mat")); correlation_distribution_horridge = correlation_distribution_horridge.correlation_distribution_horridge;

mouse_to_exclude = 0;
%mouse_to_exclude = 2;
%mouse_to_exclude = 4:9;

%%

threshold_before = zeros(numel(setdiff(1:size(raw_after,1),mouse_to_exclude)),1);
threshold_between = zeros(numel(setdiff(1:size(raw_after,1),mouse_to_exclude)),1);
threshold_after = zeros(numel(setdiff(1:size(raw_after,1),mouse_to_exclude)),1);
threshold_horridge = zeros(numel(setdiff(1:size(raw_after,1),mouse_to_exclude)),1);
thresholded_before = zeros(size(raw_before));
thresholded_between = cell(size(raw_between));
thresholded_after = cell(size(raw_after));
thresholded_horridge = cell(size(raw_horridge));
normalized_averages = zeros(numel(setdiff(1:size(raw_after,1),mouse_to_exclude)), 4);
numbers_above_threshold = zeros(numel(setdiff(1:size(raw_after,1),mouse_to_exclude)), 4);
for i = setdiff(1:size(raw_after,1),mouse_to_exclude)
    % Sort the distribution by values
    sorted_distribution_before = sortrows(correlation_distribution_before{i}, 1);
    sorted_distribution_between = sortrows(correlation_distribution_between{i}, 1);
    sorted_distribution_after = sortrows(correlation_distribution_after{i}, 1);
    sorted_distribution_horridge = sortrows(correlation_distribution_horridge{i}, 1);

    % Compute cumulative sum of frequencies
    cumulative_sum_before = cumsum(sorted_distribution_before(:, 2));
    cumulative_sum_between = cumsum(sorted_distribution_between(:, 2));
    cumulative_sum_after = cumsum(sorted_distribution_after(:, 2));
    cumulative_sum_horridge = cumsum(sorted_distribution_horridge(:, 2));

    % Find the index where cumulative sum exceeds 95% of total frequency
    total_frequency_before = sum(sorted_distribution_before(:, 2));
    total_frequency_between = sum(sorted_distribution_between(:, 2));
    total_frequency_after = sum(sorted_distribution_after(:, 2));
    total_frequency_horridge = sum(sorted_distribution_horridge(:, 2));
    
    percentage = 0.99;
    percentile_index_before = find(cumulative_sum_before >= percentage * total_frequency_before, 1);
    percentile_index_between = find(cumulative_sum_between >= percentage * total_frequency_between, 1);
    percentile_index_after = find(cumulative_sum_after >= percentage * total_frequency_after, 1);
    percentile_index_horridge = find(cumulative_sum_horridge >= percentage * total_frequency_horridge, 1);

    % Extract the value at the 99.99th percentile
    threshold_before(i) = sorted_distribution_before(percentile_index_before, 1);
    threshold_between(i) = sorted_distribution_between(percentile_index_between, 1);
    threshold_after(i) = sorted_distribution_after(percentile_index_after, 1);
    threshold_horridge(i) = sorted_distribution_horridge(percentile_index_horridge, 1);

    % apply thresholds
    cur_before = raw_before(i,:);
    cur_before(cur_before < threshold_before(i)) = 0;
    thresholded_before(i,:) = cur_before;

    thresholded_between{i} = raw_between{i};
    thresholded_between{i}(thresholded_between{i} < threshold_between(i)) = 0;
    
    thresholded_after{i} = raw_after{i};
    thresholded_after{i}(thresholded_after{i} < threshold_after(i)) = 0;

    thresholded_horridge{i} = raw_horridge{i};
    thresholded_horridge{i}(thresholded_horridge{i} < threshold_horridge(i)) = 0;
   

    % calculate averages
    normalized_averages(i,1) = (sum(thresholded_before(i,:))) / numel(thresholded_before(i,:));
    normalized_averages(i,2) = (sum(thresholded_between{i})) / numel(thresholded_between{i});
    normalized_averages(i,3) = (sum(thresholded_after{i})) / numel(thresholded_after{i});
    normalized_averages(i,4) = (sum(thresholded_horridge{i})) / numel(thresholded_horridge{i});

    numbers_above_threshold(i,1) = sum(thresholded_before(i,:)>0)/size(thresholded_before(i,:),2);
    numbers_above_threshold(i,2) = sum(thresholded_between{i}>0)/size(thresholded_between{i},2);
    numbers_above_threshold(i,3) = sum(thresholded_after{i}>0)/size(thresholded_after{i},2);
    numbers_above_threshold(i,4) = sum(thresholded_horridge{i}>0)/size(thresholded_horridge{i},2);
end

%% figures

figure
boxplot([normalized_averages(:,1),normalized_averages(:,2),normalized_averages(:,3), normalized_averages(:,4)], 'Labels', {'Baseline', 'Experiment', 'Rest', 'Horridge'})
hold on
scatter(ones(size(normalized_averages(:,1),1)),normalized_averages(:,1), 'filled', 'blue')
scatter(ones(size(normalized_averages(:,2),1))*2,normalized_averages(:,2), 'filled', 'blue')
scatter(ones(size(normalized_averages(:,3),1))*3,normalized_averages(:,3), 'filled', 'blue')
scatter(ones(size(normalized_averages(:,4),1))*4,normalized_averages(:,4), 'filled', 'blue')
line([ones(size(normalized_averages(:,1))), ones(size(normalized_averages(:,2)))*2]',[normalized_averages(:,1), normalized_averages(:,2)]','Color','green')
line([ones(size(normalized_averages(:,2)))*2, ones(size(normalized_averages(:,3)))*3]',[normalized_averages(:,2), normalized_averages(:,3)]','Color','green')
line([ones(size(normalized_averages(:,3)))*3, ones(size(normalized_averages(:,4)))*4]',[normalized_averages(:,3), normalized_averages(:,4)]','Color','green')


figure
boxplot([normalized_averages(:,1), normalized_averages(:,3)], 'Labels', {'Baseline', 'Rest'})
hold on
scatter(ones(size(normalized_averages(:,1),1)),normalized_averages(:,1), 'filled', 'blue')
scatter(ones(size(normalized_averages(:,3),1))*2,normalized_averages(:,3), 'filled', 'blue')
line([ones(size(normalized_averages(:,1))), ones(size(normalized_averages(:,3)))*2]',[normalized_averages(:,1), normalized_averages(:,3)]','Color','green')




