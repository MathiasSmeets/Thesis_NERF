clear;clc;close all
if ispc
    volume_base = '\\nerffs13';
    volume_base2 = '\\nerffs17';
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/';
    volume_base2 = '/mnt/takeokalab/';
end
path_to_raw = "takeokalabwip2023/Mathias/switch_data/tabled_data";
path_to_data = "takeokalabwip2023/Mathias/switch_data/data_after_stimulus";
path_to_clusters = "takeokalabwip2023/Mathias/switch_data/clusters";
path_to_noi = "takeokalabwip2023/Mathias/switch_data/neurons_of_interest";
path_to_correlations = "takeokalabwip2023/Mathias/switch_data/correlations";
stimulus_data_m = load(fullfile(volume_base2, path_to_data, "after_stimulus_data_m_horridge.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_m = stimulus_data_m(1:9,:);

load("X:\Mathias\switch_data\neurons_of_interest\neurons_of_interest_horridge_m.mat")
load("X:\Mathias\switch_data\neurons_of_interest\inhibited_horridge_m.mat")

sos_results_m = load("X:\Mathias\switch_data\sos_data\sos_results_m.mat");
sos_results_m = sos_results_m.sos_results_m;
last_interval_data = load("X:\Mathias\switch_data\data_after_stimulus\last_interval_data.mat");
last_interval_data = last_interval_data.last_interval_data;
mouse = 3;
%%

%noi = get_neurons_of_interest(stimulus_data_m{mouse,1}, output_m, inhibited_m, 1);

% cur_sum = zeros(size(stimulus_data_m{1,1}));
% not_norm = zeros(size(stimulus_data_m{1,1}));
% for i = 1:last_interval_data(mouse)
%     cur_std = sos_results_m{mouse,2};
%     cur_std(cur_std==0) = 0.01;
%     cur_sum = cur_sum + ((stimulus_data_m{mouse,i} - sos_results_m{mouse,1}) ./ cur_std);
%     not_norm = not_norm + stimulus_data_m{mouse,i};
% end
% cur_sum = cur_sum./last_interval_data(mouse);
% not_norm = not_norm./last_interval_data(mouse);
% 
% %cur_sum = cur_sum(noi,:);
% figure
% heatmap(cur_sum);grid('off')
% figure
% heatmap(not_norm);grid('off')
% 
% load("X:\Mathias\switch_data\correlations\template_cluster_m.mat")

%%
sos_results_m = load("X:\Mathias\switch_data\sos_data\sos_results_m.mat");
sos_results_m = sos_results_m.sos_results_m;
cur_std = sos_results_m{3,2};
cur_std(cur_std==0) = 0.01;

figure
heatmap((stimulus_data_m{3,200}-sos_results_m{3,1}) ./ cur_std)
heatmap(stimulus_data_m{3,200});grid('off')