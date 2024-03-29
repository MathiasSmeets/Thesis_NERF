clear;clc;close all;

%% learner mice

foot_position2 = load("X:\Mathias\switch_data\foot_position\switch\pair6\ZmasterD1P6.mat");foot_position2 = foot_position2.ZMaster(:,1); % 220819 --> 2
foot_position3 = load("X:\Mathias\switch_data\foot_position\switch\pair7\ZmasterD1P7.mat");foot_position3 = foot_position3.ZMaster(:,1); % 220831 --> 3
foot_position4 = load("X:\Mathias\switch_data\foot_position\switch\pair8\ZmasterD1P8.mat");foot_position4 = foot_position4.ZMaster(:,1); % 220821 --> 4

threshold = 3;
percentage2 = sum(foot_position2>threshold,"all")/numel(foot_position2);
percentage3 = sum(foot_position3>threshold,"all")/numel(foot_position3);
percentage4 = sum(foot_position4>threshold,"all")/numel(foot_position4);

% correlation_improvement2 = (mean(correlation_after_m{2}) - mean(correlation_before_m(2))) / mean(correlation_between_m{2});
% correlation_improvement3 = (mean(correlation_after_m{3}) - mean(correlation_before_m(3))) / mean(correlation_between_m{3});
% correlation_improvement4 = (mean(correlation_after_m{4}) - mean(correlation_before_m(4))) / mean(correlation_between_m{4});

%% testing
switch_wanted = 0;
if switch_wanted
    added_string = "_switchbased";
else
    added_string = "";
end

load("X:\Mathias\01\01\template_m"+added_string+".mat");
load("X:\Mathias\01\01\threshold_m"+added_string+".mat");
load("X:\Mathias\01\01\thresholded_after_m"+added_string+".mat");
load("X:\Mathias\01\01\raw_before_m"+added_string+".mat")
load("X:\Mathias\01\01\raw_between_m"+added_string+".mat")
load("X:\Mathias\01\01\raw_after_m"+added_string+".mat")
load("X:\Mathias\01\01\raw_horridge_m"+added_string+".mat")

figure
threshold_max = zeros(9,1);
for i = 1:9
subplot(3,3,i)
plot(cur_correlation_after{i})
cur_max = maxk(cur_correlation_after{i},5);
threshold_max(i) = cur_max(end)*0.6;
ylim([0 1.3])
hold on
plot(ones(size(cur_correlation_after{i}))*threshold(i),'Linewidth',3)
%plot(ones(size(cur_correlation_after{i}))*mean(threshold),'Linewidth',3)
%plot(ones(size(cur_correlation_after{i}))*threshold_max(i),'Linewidth',3)

end
%%
alt_after = cur_correlation_after;
alt_horridge = cur_correlation_horridge;
alt_between =  cur_correlation_between;
percentage_in_cluster = zeros(9,2);
alt_before = cell(9,1);
for i = 1:9
    alt_between{i}(alt_between{i}<threshold(i)) = [];
    alt_horridge{i}(alt_horridge{i}<threshold(i)) = [];
    alt_after{i}(alt_after{i}<threshold(i)) = [];
    alt_before{i} = cur_correlation_before(i,:);
    alt_before{i}(alt_before{i}<threshold(i)) = [];
    percentage_in_cluster(i,1) = numel(alt_before{i})/numel(cur_correlation_before(i,:));
    percentage_in_cluster(i,2) = numel(alt_between{i})/numel(cur_correlation_between{i});
    percentage_in_cluster(i,3) = numel(alt_after{i})/numel(cur_correlation_after{i});
    percentage_in_cluster(i,4) = numel(alt_horridge{i})/numel(cur_correlation_horridge{i});
end

figure
boxplot([percentage_in_cluster(:,1), percentage_in_cluster(:,3)], 'Labels', {'Baseline', 'Rest'})
hold on
scatter(ones(size(percentage_in_cluster(:,1),1)),percentage_in_cluster(:,1), 'filled', 'blue')
scatter(ones(size(percentage_in_cluster(:,3),1))*2,percentage_in_cluster(:,3), 'filled', 'blue')
line([ones(size(percentage_in_cluster(:,1))), ones(size(percentage_in_cluster(:,3)))*2]',[percentage_in_cluster(:,1), percentage_in_cluster(:,3)]','Color','green')

figure
boxplot([percentage_in_cluster(:,2), percentage_in_cluster(:,4)], 'Labels', {'Conditioning', 'Horridge'})
hold on
scatter(ones(size(percentage_in_cluster(:,2),1)),percentage_in_cluster(:,2), 'filled', 'blue')
scatter(ones(size(percentage_in_cluster(:,4),1))*2,percentage_in_cluster(:,4), 'filled', 'blue')
line([ones(size(percentage_in_cluster(:,2))), ones(size(percentage_in_cluster(:,4)))*2]',[percentage_in_cluster(:,2), percentage_in_cluster(:,4)]','Color','green')




