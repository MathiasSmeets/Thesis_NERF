clear;clc;close all;

correlation_path = "X:\Mathias\switch_data\correlations";
correlation_before_m = load(fullfile(correlation_path, "thresholded_before_m.mat"));correlation_before_m = correlation_before_m.thresholded_before;
correlation_between_m = load(fullfile(correlation_path, "thresholded_between_m.mat"));correlation_between_m = correlation_between_m.thresholded_between;
correlation_after_m = load(fullfile(correlation_path, "thresholded_after_m.mat"));correlation_after_m = correlation_after_m.thresholded_after;
correlation_before_y = load(fullfile(correlation_path, "thresholded_before_y.mat"));correlation_before_y = correlation_before_y.thresholded_before;
correlation_between_y = load(fullfile(correlation_path, "thresholded_between_y.mat"));correlation_between_y = correlation_between_y.thresholded_between;
correlation_after_y = load(fullfile(correlation_path, "thresholded_after_y.mat"));correlation_after_y = correlation_after_y.thresholded_after;

%% learner mice

for i=1:9
    figure
    cur_cor = movmean(correlation_after_m{i},100*46);
    plot(cur_cor)
end

foot_position2 = load("X:\Mathias\switch_data\foot_position\switch\pair6\ZmasterD1P6.mat");foot_position2 = foot_position2.ZMaster(:,1); % 220819 --> 2
foot_position3 = load("X:\Mathias\switch_data\foot_position\switch\pair7\ZmasterD1P7.mat");foot_position3 = foot_position3.ZMaster(:,1); % 220831 --> 3
foot_position4 = load("X:\Mathias\switch_data\foot_position\switch\pair8\ZmasterD1P8.mat");foot_position4 = foot_position4.ZMaster(:,1); % 220821 --> 4

threshold = 3;
percentage2 = sum(foot_position2>threshold,"all")/numel(foot_position2);
percentage3 = sum(foot_position3>threshold,"all")/numel(foot_position3);
percentage4 = sum(foot_position4>threshold,"all")/numel(foot_position4);

correlation_improvement2 = (mean(correlation_after_m{2}) - mean(correlation_before_m(2))) / mean(correlation_between_m{2});
correlation_improvement3 = (mean(correlation_after_m{3}) - mean(correlation_before_m(3))) / mean(correlation_between_m{3});
correlation_improvement4 = (mean(correlation_after_m{4}) - mean(correlation_before_m(4))) / mean(correlation_between_m{4});

%% testing
load("X:\Mathias\01\01\template_m.mat");
load("X:\Mathias\01\01\threshold_m.mat");
load("X:\Mathias\01\01\thresholded_after_m.mat");
load("X:\Mathias\01\01\raw_after_m.mat")
figure
threshold_max = zeros(9,1);
for i = 1:9
subplot(3,3,i)
plot(cur_correlation_after{i})
cur_max = maxk(cur_correlation_after{i},5);
threshold_max(i) = cur_max(end)*0.6;
ylim([0 1.5])
hold on
plot(ones(size(cur_correlation_after{i}))*threshold(i),'Linewidth',3)
plot(ones(size(cur_correlation_after{i}))*mean(threshold),'Linewidth',3)
plot(ones(size(cur_correlation_after{i}))*threshold_max(i),'Linewidth',3)

end
%%
alt_after = cur_correlation_after;
percentage_in_cluster = zeros(9,2);
alt_before = cell(9,1);
for i = 1:9
    alt_after{i}(alt_after{i}<threshold(i)) = [];
    alt_before{i} = cur_correlation_before(i,:);
    alt_before{i}(alt_before{i}<threshold(i)) = [];
    percentage_in_cluster(i,1) = numel(alt_before{i})/numel(cur_correlation_before(i,:));
    percentage_in_cluster(i,2) = numel(alt_after{i})/numel(cur_correlation_after{i});
end

figure
boxplot([percentage_in_cluster(:,1), percentage_in_cluster(:,2)], 'Labels', {'Baseline', 'Rest'})
hold on
scatter(ones(size(percentage_in_cluster(:,1),1)),percentage_in_cluster(:,1), 'filled', 'blue')
scatter(ones(size(percentage_in_cluster(:,2),1))*2,percentage_in_cluster(:,2), 'filled', 'blue')
line([ones(size(percentage_in_cluster(:,1))), ones(size(percentage_in_cluster(:,2)))*2]',[percentage_in_cluster(:,1), percentage_in_cluster(:,2)]','Color','green')




