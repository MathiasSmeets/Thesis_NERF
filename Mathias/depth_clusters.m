clear;clc;close all;

load("X:\Mathias\switch_data\correlations\template_cluster_m.mat")
load("X:\Mathias\switch_data\correlations\template_smoothed_3_m.mat")
load("X:\Mathias\switch_data\depth\depth_m.mat")

depths = cell(1,9);
for i = 1:9
    cur_clusters = template_cluster{i};
    depths{i} = depth_m{i}(cur_clusters);
end
depths_learner_horridge = depths;
figure
hold on
for i = 1:9
    scatter(ones(size(depths{i}))*i, depths{i}, 'blue', 'filled')
end
scatter(1, 640, "green", "filled")
scatter(4,1000, "green", "filled")
set(gca, 'XAxisLocation','top', 'YDir', 'reverse')
xlim([0.5 9.5])
ylim([0 2000])
% all_depths_array = sort(cell2mat(depth_m));
% depths_array = sort(cell2mat(depths));
% figure
% violin(depths_array')
% %xlim([0 max(all_depths_array)])clear
% 
% 
% figure
% violin(all_depths_array')
% %xlim([0 max(all_depths_array)])
% 
% peak_spike = cell(1,9);
% for i = 1:9
%     [~,peak_spike{i}] = max(template_smoothed{i},[],2);
% end
% 
% figure
% hold on
% for i = 1:9
%     subplot(3,3,i)
%     scatter(peak_spike{i}, depths{i},'blue', 'filled')
%     hold on
% 
%     % % draw lines in order of max spike
%     % for j = 1:i-1
%     %     line([peak_spike{j}, peak_spike{j}]',[results_y([1,3:9],2), results_y([1,3:9],3)]','Color','green')
%     % end
% end

load("C:\Users\Mathi\Downloads\frM_20220802_11corrected(4).mat")
recording_1 = find(frM.Recording==28 & frM.optotagged == 1) - find(frM.Recording==28,1);
recording_2 = find(frM.Recording==29 & frM.optotagged == 1) - find(frM.Recording==29,1);
recording_4 = find(frM.Recording==30 & frM.optotagged == 1) - find(frM.Recording==30,1);
recording_5 = find(frM.Recording==31 & frM.optotagged == 1) - find(frM.Recording==31,1);
recording_6 = find(frM.Recording==32 & frM.optotagged == 1) - find(frM.Recording==32,1);
recording_7 = find(frM.Recording==34 & frM.optotagged == 1) - find(frM.Recording==34,1);

% the only ptf1alpha neuron is 72 in recording 6
depths_master = depths;
%% control
load("X:\Mathias\switch_data\correlations\template_cluster_y.mat")
load("X:\Mathias\switch_data\correlations\template_smoothed_3_y.mat")
load("X:\Mathias\switch_data\depth\depth_y.mat")

depths = cell(1,9);
for i = 1:9
    cur_clusters = template_cluster{i};
    depths{i} = depth_y{i}(cur_clusters);
end
depths_control_horridge = depths;
figure
hold on
for i = 1:9
    scatter(ones(size(depths{i}))*i, depths{i}, 'blue', 'filled')
end
set(gca, 'XAxisLocation','top', 'YDir', 'reverse')
xlim([0.5 9.5])
ylim([0 2000])
depths_yoke = depths;

%% compare
array_m = [];
array_y = [];

for i = 1:size(depths_master,2)
    array_m = [array_m, depths_master{i}];
    array_y = [array_y, depths_yoke{i}];
end

figure
boxplot([array_m', array_y'])
ylim([0 2000])
set(gca, 'YDir', 'reverse')

%% rest
% 
% 
% path1 = "X:\Mathias\switch_data\clusters";
% path2 = "X:\Mathias\switch_data\data_after_stimulus";
% total_assemblies = load(fullfile(path1,"assemblies_switch_m.mat"));total_assemblies = total_assemblies.total_assemblies;
% total_neurons_of_interest = load(fullfile(path1, "neurons_of_interest_switch_m.mat"));total_neurons_of_interest = total_neurons_of_interest.total_neurons_of_interest;
% total_data = load(fullfile(path1, "data_switch_m.mat"));total_data = total_data.total_data;
% load("X:\Mathias\switch_data\correlations\template_cluster_m.mat")
% load("X:\Mathias\switch_data\correlations\template_smoothed_3_m.mat")
% load("X:\Mathias\switch_data\depth\depth_m.mat")
% 
% clusters_to_check = cell(1,9);
% for i = 1:9
%     % get last interval
%     j = size(total_data,2);
%     while isempty(total_data{i,j}) && j > 1
%         j = j-1;
%     end
%     last_interval_index = j;
% 
%     all_assemblies = {};
%     all_assemblies_count = [];
%     for k = 1:last_interval_index
%         for l = 1:length(total_assemblies{i,k})
%             cur_assembly = total_neurons_of_interest{i,k}(total_assemblies{i,k}{l});
%             idx = find(cellfun(@(x) isequal(x, cur_assembly), all_assemblies));
%             if ~isempty(idx)
%                 all_assemblies_count(idx) = all_assemblies_count(idx) + 1;
%             else
%                 all_assemblies{end+1} = cur_assembly;
%                 all_assemblies_count = [all_assemblies_count, 1];
%             end
%         end
%     end
%     counter = 1;
%     for x = 1:numel(all_assemblies)
%         cur_assembly = all_assemblies{x};
%         cur_count = all_assemblies_count(x);
%         cur_clusters = template_cluster{i};
%         common_elements = intersect(cur_assembly, cur_clusters);
%         if cur_count > 0.03*last_interval_index && numel(common_elements)>1
%             clusters_to_check{i}{counter} = cur_assembly;
%             counter = counter + 1;
%         end
%     end
% end
% 
% 
% depths_rest = cell(1,9);
% for i = 1:9
%     cur_depths = [];
%     for j = 1:numel(clusters_to_check{i})
%         cur_clusters = clusters_to_check{i}{j};
%         cur_depths = [cur_depths, depth_m{i}(cur_clusters)];
% 
%     end
%     depths_rest{i} = unique(cur_depths);
% end
% 
% figure
% hold on
% for i = 1:9
%     scatter(ones(size(depths_rest{i}))*i, depths_rest{i}, 'blue', 'filled')
% end

% set(gca, 'XAxisLocation','top', 'YDir', 'reverse')
% xlim([0.5 9.5])
% ylim([0 2000])
% title('Learner Rest')



%% control rest
% 
% path1 = "X:\Mathias\switch_data\clusters";
% path2 = "X:\Mathias\switch_data\data_after_stimulus";
% total_assemblies = load(fullfile(path1,"assemblies_switch_y.mat"));total_assemblies = total_assemblies.total_assemblies;
% total_neurons_of_interest = load(fullfile(path1, "neurons_of_interest_switch_y.mat"));total_neurons_of_interest = total_neurons_of_interest.total_neurons_of_interest;
% total_data = load(fullfile(path1, "data_switch_y.mat"));total_data = total_data.total_data;
% load("X:\Mathias\switch_data\correlations\template_cluster_y.mat")
% load("X:\Mathias\switch_data\correlations\template_smoothed_3_y.mat")
% load("X:\Mathias\switch_data\depth\depth_y.mat")
% 
% clusters_to_check = cell(1,9);
% for i = 1:9
%     % get last interval
%     j = size(total_data,2);
%     while isempty(total_data{i,j}) && j > 1
%         j = j-1;
%     end
%     last_interval_index = j;
% 
%     all_assemblies = {};
%     all_assemblies_count = [];
%     for k = 1:last_interval_index
%         for l = 1:length(total_assemblies{i,k})
%             cur_assembly = total_neurons_of_interest{i,k}(total_assemblies{i,k}{l});
%             idx = find(cellfun(@(x) isequal(x, cur_assembly), all_assemblies));
%             if ~isempty(idx)
%                 all_assemblies_count(idx) = all_assemblies_count(idx) + 1;
%             else
%                 all_assemblies{end+1} = cur_assembly;
%                 all_assemblies_count = [all_assemblies_count, 1];
%             end
%         end
%     end
%     counter = 1;
%     for x = 1:numel(all_assemblies)
%         cur_assembly = all_assemblies{x};
%         cur_count = all_assemblies_count(x);
%         cur_clusters = template_cluster{i};
%         common_elements = intersect(cur_assembly, cur_clusters);
%         if cur_count > 0.1*last_interval_index% && numel(common_elements)>1
%             clusters_to_check{i}{counter} = cur_assembly;
%             counter = counter + 1;
%         end
%     end
% end
% 
% 
% depths_rest_y = cell(1,9);
% for i = 1:9
%     cur_depths = [];
%     for j = 1:numel(clusters_to_check{i})
%         cur_clusters = clusters_to_check{i}{j};
%         cur_depths = [cur_depths, depth_y{i}(cur_clusters)];
% 
%     end
%     depths_rest_y{i} = unique(cur_depths);
% end
% 
% figure
% hold on
% for i = 1:9
%     scatter(ones(size(depths_rest_y{i}))*i, depths_rest_y{i}, 'blue', 'filled')
% end
% title('Control Rest')
% set(gca, 'XAxisLocation','top', 'YDir', 'reverse')
% xlim([0.5 9.5])
% ylim([0 2000])
% 
% array_m = [];
% array_y = [];
% 
% for i = 1:size(depths_master,2)
%     array_m = [array_m, depths_rest{i}];
%     array_y = [array_y, depths_rest_y{i}];
% end
% 
% 
% hold on
% figure
% boxplot([array_m'; array_y'], [zeros(length(array_m'),1); ones(length(array_y'),1)], 'Labels', {'Learner', 'Control'})
% ylim([0 2000])
% set(gca, 'YDir', 'reverse')


%%

array = cell2mat(depths_learner_horridge);
% Define the bin edges
binEdges = 0:250:2000;
% Count the number of values in each interval
counts = histcounts(array, binEdges);
% Define the y-axis values (midpoints of the intervals)
y = binEdges(1:end-1) + diff(binEdges)/2;
% Plot the data with switched axes
figure;
plot(counts/numel(array), y);
ylabel('Value Ranges');
xlabel('Number of Values');
title('Number of Values in Each Interval');
yticks(binEdges); % This sets the y-axis ticks to the bin edges
hold on
array = cell2mat(depths_control_horridge);
% Define the bin edges
binEdges = 0:250:2000;
% Count the number of values in each interval
counts = histcounts(array, binEdges);
% Define the y-axis values (midpoints of the intervals)
y = binEdges(1:end-1) + diff(binEdges)/2;
% Plot the data with switched axes
plot(counts/numel(array), y);
legend("Learner", "Control")

set(gca, 'YDir', 'reverse')






