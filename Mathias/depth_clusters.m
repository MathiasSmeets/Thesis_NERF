clear;clc;close all;

load("X:\Mathias\switch_data\correlations\template_cluster_m.mat")
load("X:\Mathias\switch_data\correlations\template_smoothed_3_m.mat")
load("X:\Mathias\switch_data\depth\depth_m.mat")

depths = cell(1,9);
for i = 1:9
    cur_clusters = template_cluster{i};
    depths{i} = depth_m{i}(cur_clusters);
end

figure
hold on
for i = 1:9
    scatter(ones(size(depths{i}))*i, depths{i}, 'blue', 'filled')
end
scatter(1, 640, "green", "filled")
scatter(4,1000, "green", "filled")
set(gca, 'XAxisLocation','top', 'YDir', 'reverse')
xlim([0.5 9.5])
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
