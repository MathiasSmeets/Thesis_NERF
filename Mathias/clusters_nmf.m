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

folder = fileparts(which("clusters_nmf.m"));
addpath(genpath(folder))

%% variables initialization

interval_size = size(stimulus_data_m{1,1},2);
wanted_bin_size = 1;
create_plots = false;
neuron_counter = 1;
indices = ceil((1:interval_size)/wanted_bin_size);
total_nb_assemblies = cell(size(stimulus_data_m));
total_nb_neurons = cell(size(stimulus_data_m));
total_assemblies = cell(size(stimulus_data_m));
total_activity = cell(size(stimulus_data_m));

%% nmf (NOT on zscores, only positive data)


% loop over each mouse
for k = 1:size(stimulus_data_m,1)
    % loop over each interval of this mouse
    cur_neurons_of_interest = get_neurons_of_interest(stimulus_data_m{k,1}, neurons_of_interest_m, inhibited_neurons_m, neuron_counter);
    for i = 1:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{k,i})
            
            cur_mouse = stimulus_data_m{k,i}(cur_neurons_of_interest,:);

            % transform to wanted_bin_size ms bins
            cur_mouse_fs_adjusted = zeros(size(cur_mouse,1),size(cur_mouse,2)/wanted_bin_size);
            for j = 1:size(cur_mouse,1)
                cur_mouse_fs_adjusted(j,:) = accumarray(indices',cur_mouse(j,:)',[],@sum)';
            end
            
            D=zeros(1,10);
            for iii = 1:23
                [W,H, D(iii)] = nnmf(cur_mouse_fs_adjusted, iii);
            end

        end
    end
    neuron_counter = neuron_counter + size(stimulus_data_m{k,1},1);
end