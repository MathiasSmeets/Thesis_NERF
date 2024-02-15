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

stimulus_data_m = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_m.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_y = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_y.mat"));
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

neurons_of_interest_m = load(fullfile(volume_base2, path_to_code, "output_m.mat"));
neurons_of_interest_m = neurons_of_interest_m.output_m;
inhibited_neurons_m = load(fullfile(volume_base2, path_to_code, "inhibited_m.mat"));
inhibited_neurons_m = inhibited_neurons_m.inhibited_m;

folder = fileparts(which("clusters_cpd.m"));
addpath(genpath(folder))

%% load in ica data

ica_activity = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\activity.mat"); ica_activity = ica_activity.total_activity;
ica_nb_neurons = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\nb_neurons.mat"); ica_nb_neurons = ica_nb_neurons.total_nb_neurons;
ica_assemblies = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\assemblies.mat"); ica_assemblies = ica_assemblies.total_assemblies;
ica_nb_assemblies = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\nb_assemblies.mat"); ica_nb_assemblies = ica_nb_assemblies.total_nb_assemblies;
ica_data = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\data.mat"); ica_data = ica_data.total_data;
ica_neurons_of_interest = load("X:\Mathias\cluster_output\ica_10ms_bins_10_intervals\neurons_of_interest.mat"); ica_neurons_of_interest = ica_neurons_of_interest.total_neurons_of_interest;

%% initialize data
interval_step = 10; % depends on how data was created

%% get cluster shape for each assembly

orders = cell(size(ica_activity));
% loop over each mouse
for k = 1:size(orders,1)
    % loop over each interval
    for i = 1:size(orders,2)
        cur_assembly = ica_assemblies{k,i};
        if isempty(cur_assembly)
            orders{k,i} = [];
        elseif size(cur_assembly,1) == 1
            cur_assembly = cur_assembly{1,1};
            cur_assembly = ica_neurons_of_interest{k,i}(cur_assembly);
            cur_data = ica_data{k,i};
            [pks, locs] = findpeaks(ica_activity{k,i},"NPeaks",10,"MinPeakHeight",2*mean(ica_activity{k,i}));
            candidate_templates = cell(1,length(locs));
            TOTAL = [];
            for j = 1:length(locs)
                raw_data_index = (i-1)*interval_step + ceil((locs(j))/7);
                position_in_data = mod(locs(j)-1,7)*10+1;
                cur_raw_data = stimulus_data_m{k,raw_data_index};
                cur_assembly_data = cur_raw_data(cur_assembly, position_in_data:position_in_data+9);

                % find the indices of the first and last non-zero elements in each row
                first_nonzero = find(any(cur_assembly_data, 1), 1, 'first');
                last_nonzero = find(any(cur_assembly_data, 1), 1, 'last');
                candidate_templates{1,j} = cur_assembly_data(:, first_nonzero:last_nonzero);

            end
            % find the order of the neurons firing
            orders{k,i} = find_order_neurons(candidate_templates);
            disp("k")
            % go to each peak and get template there (use median to get total template?)
            % maybe first calculate median length
            % stack time in data in 3D on top of each other
            % template = median(stacked_templates,3);
        else
            % do the same but first need a step where you detect which
            % assembly is which peak based on ica_assemblies abs size
        end
    end
end

[1;
 3;
 2;
 4];
