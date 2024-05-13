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
path_to_stimulus_data = "takeokalabwip2023/Mathias/switch_data/data_after_stimulus";

raw_before = load(fullfile(volume_base2, path_to_correlations, "raw_before_m.mat")); raw_before = raw_before.cur_correlation_before;
raw_between = load(fullfile(volume_base2, path_to_correlations, "raw_between_m.mat")); raw_between = raw_between.cur_correlation_between;
raw_after = load(fullfile(volume_base2, path_to_correlations, "raw_after_m.mat")); raw_after = raw_after.cur_correlation_after;
raw_horridge = load(fullfile(volume_base2, path_to_correlations, "raw_horridge_m.mat")); raw_horridge = raw_horridge.cur_correlation_horridge;

correlation_distribution_before = load(fullfile(volume_base2, path_to_correlations, "correlation_distribution_before_m.mat")); correlation_distribution_before = correlation_distribution_before.correlation_distribution_before;
correlation_distribution_between = load(fullfile(volume_base2, path_to_correlations, "correlation_distribution_between_m.mat")); correlation_distribution_between = correlation_distribution_between.correlation_distribution_between;
correlation_distribution_after = load(fullfile(volume_base2, path_to_correlations, "correlation_distribution_after_m.mat")); correlation_distribution_after = correlation_distribution_after.correlation_distribution_after;
correlation_distribution_horridge = load(fullfile(volume_base2, path_to_correlations, "correlation_distribution_horridge_m.mat")); correlation_distribution_horridge = correlation_distribution_horridge.correlation_distribution_horridge;

%last_interval_data = load(fullfile(volume_base2, path_to_stimulus_data, "last_interval_data.mat")); last_interval_data = last_interval_data.last_interval_data;

mouse_to_exclude = 0; % m
%mouse_to_exclude = 2; % y
%mouse_to_exclude = 4:9; % np2

stimulations = cell(1,9);
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_1.mat");stimulations{1} = cur_stim.frM_stim.horridge;
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_2.mat");stimulations{2} = cur_stim.frM_stim.horridge;
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_3.mat");stimulations{3} = cur_stim.frM_stim.horridge;
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_4.mat");stimulations{4} = cur_stim.frM_stim.horridge;
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_5.mat");stimulations{5} = cur_stim.frM_stim.horridge;
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_6.mat");stimulations{6} = cur_stim.frM_stim.horridge;
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_7.mat");stimulations{7} = cur_stim.frM_stim.horridge;
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_8.mat");stimulations{8} = cur_stim.frM_stim.horridge;
cur_stim = load("X:\Mathias\switch_data\tables\frM_stim_switched_9.mat");stimulations{9} = cur_stim.frM_stim.horridge;
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

    
    % figure
    % subplot(1,2,1)
    % lines = 1000*(stimulations{i}{1}-stimulations{i}{1}(1));
    % xline(lines)
    % hold on
    % plot(movmean(thresholded_between{i},10000))
    % title("Horridge")
    % subplot(1,2,2)
    % plot(movmean(thresholded_after{i},100000))
    % title("Rest")
end


%%

% average_correlation_between = zeros(size(thresholded_between{1}));
% std_correlation_between = zeros(size(thresholded_between{1}));
% average_correlation_after = zeros(size(thresholded_after{2}));
% std_correlation_after = zeros(size(thresholded_after{2}));
% 
% for i = 1:numel(average_correlation_between)
%     average_correlation_between(i) = mean([thresholded_between{1}(i),thresholded_between{2}(i),thresholded_between{3}(i),thresholded_between{4}(i),thresholded_between{5}(i),thresholded_between{6}(i),thresholded_between{7}(i),thresholded_between{8}(i),thresholded_between{9}(i)]);
%     std_correlation_between(i) = std([thresholded_between{1}(i),thresholded_between{2}(i),thresholded_between{3}(i),thresholded_between{4}(i),thresholded_between{5}(i),thresholded_between{6}(i),thresholded_between{7}(i),thresholded_between{8}(i),thresholded_between{9}(i)]);
% end
% for i = 1:numel(average_correlation_after)
%     average_correlation_after(i) = mean([thresholded_after{1}(i),thresholded_after{2}(i),thresholded_after{3}(i),thresholded_after{4}(i),thresholded_after{5}(i),thresholded_after{6}(i),thresholded_after{7}(i),thresholded_after{8}(i),thresholded_after{9}(i)]);
%     std_correlation_after(i) = std([thresholded_after{1}(i),thresholded_after{2}(i),thresholded_after{3}(i),thresholded_after{4}(i),thresholded_after{5}(i),thresholded_after{6}(i),thresholded_after{7}(i),thresholded_after{8}(i),thresholded_after{9}(i)]);
% end
% figure
% plot(movmean(average_correlation_between,10000))
% hold on
% plot(movmean(average_correlation_between+std_correlation_between,10000))
% plot(movmean(average_correlation_between-std_correlation_between,10000))
% 
% figure
% plot(movmean(average_correlation_after,100000))
% hold on
% plot(movmean(average_correlation_after+std_correlation_after,100000))
% plot(movmean(average_correlation_after-std_correlation_after,100000))

%%
differences = normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3)-normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1);

%% figures

figure
boxplot([normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1),normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),2),normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3), normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),4)], 'Labels', {'Baseline', 'Horridge', 'Rest', 'Switch'})
hold on
scatter(ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1),1)),normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1), 'filled', 'blue')
scatter(ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),2),1))*2,normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),2), 'filled', 'blue')
scatter(ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3),1))*3,normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3), 'filled', 'blue')
scatter(ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),4),1))*4,normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),4), 'filled', 'blue')
line([ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1))), ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),2)))*2]',[normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1), normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),2)]','Color','green')
line([ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),2)))*2, ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3)))*3]',[normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),2), normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3)]','Color','green')
line([ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3)))*3, ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),4)))*4]',[normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3), normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),4)]','Color','green')

figure
boxplot([normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1), normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3)], 'Labels', {'Baseline', 'Rest'})
hold on
scatter(ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1),1)),normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1), 'filled', 'blue')
scatter(ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3),1))*2,normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3), 'filled', 'blue')
line([ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1))), ones(size(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3)))*2]',[normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1), normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3)]','Color','green')


p_value = signrank(normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),1),normalized_averages(setdiff(1:size(raw_after,1),mouse_to_exclude),3));
%disp("p-value baseline-rest: "+p_value)

save(fullfile(volume_base2, path_to_correlations, "results_m.mat"), "normalized_averages")


%% results without previous

results_m = load("X:\Mathias\switch_data\correlations\results_m.mat");results_m = results_m.normalized_averages;
results_np2 = load("X:\Mathias\switch_data\correlations\results_np2.mat");results_np2 = results_np2.normalized_averages;
results_y = load("X:\Mathias\switch_data\correlations\results_y.mat");results_y = results_y.normalized_averages;
results_m = [results_m;results_np2];
figure
data = [results_m(:,1); results_m(:,3); results_y([1,3:9],1); results_y([1,3:9],3)];
g = [zeros(length(results_m(:,1)),1); ones(length(results_m(:,3)),1); 2*ones(length(results_y([1,3:9],1)),1); 3*ones(length(results_y([1,3:9],3)),1)];
boxplot(data, g, 'Labels', {'Baseline', 'Rest', 'Baseline', 'Rest'})
hold on
scatter(ones(size(results_m(:,1),1)),results_m(:,1), 'filled', 'blue')
scatter(ones(size(results_m(:,3),1))*2,results_m(:,3), 'filled', 'blue')
line([ones(size(results_m(:,1),1),1), ones(size(results_m(:,3),1),1)*2]',[results_m(:,1), results_m(:,3)]','Color','green')
%figure
%boxplot([results_y([1,3:9],1), results_y([1,3:9],3)], 'Labels', {'Baseline', 'Rest'})
%hold on
scatter(ones(size(results_y([1,3:9],1),1))*3,results_y([1,3:9],1), 'filled', 'blue')
scatter(ones(size(results_y([1,3:9],3),1))*4,results_y([1,3:9],3), 'filled', 'blue')
line([ones(size(results_y([1,3:9],1),1),1)*3, ones(size(results_y([1,3:9],3),1),1)*4]',[results_y([1,3:9],1), results_y([1,3:9],3)]','Color','green')

p_m = signrank(results_m(:,1), results_m(:,3));
p_y = signrank(results_y([1,3:9],1), results_y([1,3:9],3));

disp("P-value Learner: "+ p_m)
disp("P-value Control: "+ p_y)


figure
data2 = [results_m(:,1); results_m(:,2);results_m(:,3);results_y([1,3:9],1);results_y([1,3:9],2);results_y([1,3:9],3)];
g2 = [zeros(length(results_m(:,1)),1); ones(length(results_m(:,2)),1); 2*ones(length(results_m(:,3)),1); 3*ones(length(results_y([1,3:9],1)),1); 4*ones(length(results_y([1,3:9],2)),1); 5*ones(length(results_y([1,3:9],3)),1)];
boxplot(data2, g2, 'Labels',{'Baseline', 'Horridge','Rest','Baseline', 'Horridge','Rest'})
hold on
scatter(ones(size(results_m(:,1),1)), results_m(:,1), 'filled', 'blue')
scatter(ones(size(results_m(:,2),1))*2, results_m(:,2), 'filled', 'blue')
scatter(ones(size(results_m(:,3),1))*3, results_m(:,3), 'filled', 'blue')
line([ones(size(results_m(:,1),1),1), ones(size(results_m(:,2),1),1)*2]',[results_m(:,1), results_m(:,2)]','Color','green')
line([ones(size(results_m(:,2),1),1)*2, ones(size(results_m(:,3),1),1)*3]',[results_m(:,2), results_m(:,3)]','Color','green')

scatter(ones(size(results_y([1,3:9],1),1))*4, results_y([1,3:9],1), 'filled', 'blue')
scatter(ones(size(results_y([1,3:9],2),1))*5, results_y([1,3:9],2), 'filled', 'blue')
scatter(ones(size(results_y([1,3:9],3),1))*6, results_y([1,3:9],3), 'filled', 'blue')
line([ones(size(results_y([1,3:9],1),1),1)*4, ones(size(results_y([1,3:9],2),1),1)*5]',[results_y([1,3:9],1), results_y([1,3:9],2)]','Color','green')
line([ones(size(results_y([1,3:9],2),1),1)*5, ones(size(results_y([1,3:9],1),1),1)*6]',[results_y([1,3:9],2), results_y([1,3:9],3)]','Color','green')

