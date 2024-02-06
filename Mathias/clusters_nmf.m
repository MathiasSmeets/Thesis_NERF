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

sos_results_m = load(fullfile(volume_base2, path_to_code, "sos_results_m.mat"));
sos_results_m = sos_results_m.sos_results_m;

%% variables initialization

interval_size = size(stimulus_data_m{1,1},2);
wanted_bin_size = 10;
create_plots = false;
neuron_counter = 1;
interval_step = 10;
indices = ceil((1:interval_size)/wanted_bin_size);
total_nb_assemblies = cell(size(stimulus_data_m));
total_nb_neurons = cell(size(stimulus_data_m));
total_assemblies = cell(size(stimulus_data_m));
total_activity = cell(size(stimulus_data_m));

%% nmf (NOT on zscores, only positive data)

% the rank r of the factorization is generally chosen so that (n+m)r < nm  (https://www.cs.columbia.edu/~blei/fogm/2020F/readings/LeeSeung1999.pdf) also nice source regarding constraints different methods

% loop over each mouse
for k = 1:size(stimulus_data_m,1)
    % loop over each interval of this mouse
    cur_neurons_of_interest = get_neurons_of_interest(stimulus_data_m{k,1}, neurons_of_interest_m, inhibited_neurons_m, neuron_counter);
    cur_total_mouse = [];
    for i = 1:size(stimulus_data_m,2)
        if ~isempty(stimulus_data_m{k,i})
            for ii = i:i+interval_step-1

                cur_mouse = stimulus_data_m{k,ii}(cur_neurons_of_interest,:);

                % transform to wanted_bin_size ms bins
                cur_mouse_fs_adjusted = zeros(size(cur_mouse,1),size(cur_mouse,2)/wanted_bin_size);
                for j = 1:size(cur_mouse,1)
                    cur_mouse_fs_adjusted(j,:) = accumarray(indices',cur_mouse(j,:)',[],@sum)';
                end

                cur_total_mouse = [cur_total_mouse, cur_mouse_fs_adjusted];
            end

            % remove neurons that are not active
            for jj = size(cur_total_mouse, 1):-1:1
                if all(cur_total_mouse(jj,:) == 0)
                    cur_total_mouse(jj,:) = [];
                    cur_neurons_of_interest(jj) = [];
                end
            end

            % calculate zscores
            cur_std = sos_results_m{k,2};
            cur_std(cur_std==0) = 0.01;
            cur_total_mouse_zscore = (cur_total_mouse - sos_results_m{k,1}(cur_neurons_of_interest,:)) ./ cur_std(cur_neurons_of_interest,:);
            cur_total_mouse_zscore(cur_total_mouse_zscore<0) = 0;
            cur_total_mouse = cur_total_mouse_zscore;

            rank = determine_rank(cur_total_mouse);
            [W,H] = nnmf(cur_total_mouse, rank);
            %[W_seq, H_seq] = seqNMF(cur_total_mouse,'K',3, 'L', 10,'lambda', 0,'showplot',1);
        end

    end
    neuron_counter = neuron_counter + size(stimulus_data_m{k,1},1);
    disp(k)
end


