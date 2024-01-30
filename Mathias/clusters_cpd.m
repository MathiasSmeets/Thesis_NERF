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

%% variables initialization

interval_size = size(stimulus_data_m{1,1},2);

%% create 3th order tensors

%% 1ms bins, all neurons --> does not really work

% loop over each mouse
for k = 1:size(stimulus_data_m, 1)
    cur_tensor_1ms_all_neurons = stimulus_data_m{k,1};
    % loop over each interval of this mouse
    for i = 2:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{k,i})
            cur_tensor_1ms_all_neurons(:,:,i) = stimulus_data_m{k,i};
        end
    end
end

%% 1ms bins, neurons of interest

% loop over each mouse
for k = 1:size(stimulus_data_m, 1)
    cur_neurons_of_interest = get_neurons_of_interest(stimulus_data_m{k,1}, neurons_of_interest_m, inhibited_neurons_m, neuron_counter);
    cur_tensor_1ms_neurons_oi = stimulus_data_m{k,1}(cur_neurons_of_interest,:);
    for i = 2:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{k,i})
            cur_mouse = stimulus_data_m{k,i}(cur_neurons_of_interest,:);
            cur_tensor_1ms_neurons_oi(:,:,i) = stimulus_data_m{};
        end
    end
end












