clear;clc;close all;

clusters_m = load("X:\Mathias\switch_data\correlations\template_cluster_m.mat"); clusters_m = clusters_m.template_cluster;
clusters_y = load("X:\Mathias\switch_data\correlations\template_cluster_y.mat"); clusters_y = clusters_y.template_cluster;
load("X:\Mathias\switch_data\data_after_stimulus\after_stimulus_data_m_horridge.mat")
load("X:\Mathias\switch_data\data_after_stimulus\after_stimulus_data_y_horridge.mat")
output_m = load("X:\Mathias\switch_data\neurons_of_interest\neurons_of_interest_horridge_m.mat"); output_m = output_m.output_m;
inhibited_m = load("X:\Mathias\switch_data\neurons_of_interest\inhibited_horridge_m.mat"); inhibited_m = inhibited_m.inhibited_m;
output_y = load("X:\Mathias\switch_data\neurons_of_interest\neurons_of_interest_horridge_y.mat"); output_y = output_y.output_m;
inhibited_y = load("X:\Mathias\switch_data\neurons_of_interest\inhibited_horridge_y.mat"); inhibited_y = inhibited_y.inhibited_m;

secondary_m = cell2mat(output_m.(1){1,1});
secondary_m = unique(secondary_m(:,1));
other_m = cell2mat(output_m.(2){1,1});
other_m = unique(other_m(:,1));
neuron_counter = 0;
for i = 1:9
    cur_neurons = clusters_m{i} + neuron_counter;
    cur_secondary = numel(intersect(cur_neurons,secondary_m));
    cur_others = numel(intersect(cur_neurons,other_m));
    %cur_inhibited = numel(intersect(cur_neurons,inhibited_m));

    disp("secondary: "+ cur_secondary)
    disp("others: " + cur_others)

    neuron_counter = neuron_counter + size(after_stimulus_data_m{i,1},1);
end


secondary_y = cell2mat(output_y.(1){1,1});
secondary_y = unique(secondary_y(:,1));
othery_y = cell2mat(output_y.(2){1,1});
other_y = unique(other_y(:,1));