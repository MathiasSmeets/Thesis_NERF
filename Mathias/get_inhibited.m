clear; clc; close all; % make sure to reference paper lab if i end up using this code

path_to_code = "X:\Mathias\switch_data\data_after_stimulus\";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
%other_path = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

% stimulus_data_m = load(path_to_code + "data_after_stimulus_m.mat");
% stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
% stimulus_data_y = load(path_to_code + "data_after_stimulus_y.mat");
% stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

stimulus_data_m = load(path_to_code + "after_stimulus_data_y_horridge.mat");
stimulus_data_m = stimulus_data_m.after_stimulus_data_y;
stimulus_data_m = cellfun(@double, stimulus_data_m, 'UniformOutput', false);
stimulus_data_y = load(path_to_code + "after_stimulus_data_y_horridge.mat");
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;
stimulus_data_y = cellfun(@double, stimulus_data_y, 'UniformOutput', false);


interval_size = 70;

%% get total amount of spikes + amount of spikes of each neuron
total_m = zeros(1,interval_size);
total_y = zeros(1,interval_size);

total_nb_neurons_m = 0;
total_nb_neurons_y = 0;
for i = 1:size(stimulus_data_m,1)
    total_nb_neurons_m = total_nb_neurons_m + size(stimulus_data_m{i,1},1);
end
for i = 1:size(stimulus_data_y,1)
    total_nb_neurons_y = total_nb_neurons_y + size(stimulus_data_y{i,1},1);
end
neuron_spikes_m = zeros(total_nb_neurons_m, interval_size);
neuron_spikes_y = zeros(total_nb_neurons_y, interval_size);

cur_neuron_nb_m = 1;
cur_neuron_nb_y = 1;
for i = 1:size(stimulus_data_m,1)
    for j = 1:size(stimulus_data_m,2)

        if ~isempty(stimulus_data_m{i,j})
            if size(stimulus_data_m{i,j},2) < size(total_m,2)
                total_m(1:size(stimulus_data_m{i,j},2)) = total_m(1:size(stimulus_data_m{i,j},2)) + sum(stimulus_data_m{i,j});
                neuron_spikes_m(cur_neuron_nb_m:cur_neuron_nb_m + size(stimulus_data_m{i,j},1) - 1, 1:size(stimulus_data_m{i,j},2)) = neuron_spikes_m(cur_neuron_nb_m:cur_neuron_nb_m + size(stimulus_data_m{i,j}, 1) - 1 ,1:size(stimulus_data_m{i,j},2)) + stimulus_data_m{i,j};
            else
                total_m = total_m + sum(stimulus_data_m{i,j});
                neuron_spikes_m(cur_neuron_nb_m:cur_neuron_nb_m + size(stimulus_data_m{i,j},1) - 1, :) = neuron_spikes_m(cur_neuron_nb_m:cur_neuron_nb_m + size(stimulus_data_m{i,j}, 1) - 1, :) + stimulus_data_m{i,j};
            end
        end
    end
    cur_neuron_nb_m = cur_neuron_nb_m + size(stimulus_data_m{i,1},1);
end

for i = 1:size(stimulus_data_y,1)
    for j = 1:size(stimulus_data_y,2)
            if size(stimulus_data_y{i,j},2) < size(total_y,2)
                total_y(1:size(stimulus_data_y{i,j},2)) = total_y(1:size(stimulus_data_y{i,j},2)) + sum(stimulus_data_y{i,j});
                neuron_spikes_y(cur_neuron_nb_y:cur_neuron_nb_y + size(stimulus_data_y{i,j},1) - 1, 1:size(stimulus_data_y{i,j},2)) = neuron_spikes_y(cur_neuron_nb_y:cur_neuron_nb_y + size(stimulus_data_y{i,j}, 1) - 1 ,1:size(stimulus_data_y{i,j},2)) + stimulus_data_y{i,j};
            else
                total_y = total_y + sum(stimulus_data_y{i,j});
                neuron_spikes_y(cur_neuron_nb_y:cur_neuron_nb_y + size(stimulus_data_y{i,j}, 1) - 1, :) = neuron_spikes_y(cur_neuron_nb_y:cur_neuron_nb_y + size(stimulus_data_y{i,j}, 1) - 1, :) + stimulus_data_y{i,j};
            end
    end
    cur_neuron_nb_y = cur_neuron_nb_y + size(stimulus_data_y{i,1},1);
end

%neuron_spikes_m(:,9:11) = 0;
%neuron_spikes_y(:,9:11) = 0;


%% detect inhibited neurons using cusum                                                                                     %% todo: enkel kijken naar periode na de stimulus

%decreased_FR_idx_m = find(norm_mean_after_m < norm_mean_before_m - 5*norm_std_before_m);
%decreased_FR_idx_y = find(norm_mean_after_y < norm_mean_before_y - 5*norm_std_before_y);

[peaks_m_max, loc_peaks_m_max] = max(neuron_spikes_m, [], 2);
[peaks_y_max, loc_peaks_y_max] = max(neuron_spikes_y, [], 2);
peak_interval_m = loc_peaks_m_max - 1 : loc_peaks_m_max + 1;

mean_before_m = mean(neuron_spikes_m(:,1:9), 2);
mean_after_m = mean(neuron_spikes_m(:,12:70), 2);

std_before_m = std(neuron_spikes_m(:,1:9), [], 2);
std_after_m = std(neuron_spikes_m(:,12:70), [] , 2);

% compute zscores of psth
std_before_m(std_before_m == 0) = 0.1;
neuron_spikes_m_zscore = (neuron_spikes_m - mean_before_m) ./ std_before_m;

inhibited_m = [];
ilower_m = cell(size(neuron_spikes_m_zscore,1),1);
for i = 1:size(neuron_spikes_m_zscore,1)
    %[~, ilower_m{i,1}] = cusum_(neuron_spikes_m_zscore(i,:),5,1,0,2.5);
    %if std(neuron_spikes_m_zscore(i,1:10)) == 0
    %    [~, ilower_m{i,1}] = cusum_(neuron_spikes_m_zscore(i,11:end),5,1,mean(neuron_spikes_m_zscore(i,1:10)),0.1, "all");
    %else
    [~, ilower_m{i,1}] = cusum_(neuron_spikes_m_zscore(i,12:end),5,1,0,1, "all");
    %end

    if ~isempty(ilower_m{i,1}) && 40 <= sum(neuron_spikes_m(i,1:10))
        inhibited_m = [inhibited_m, i];
    end
end

%% validation

% close all;
% for i = 1:40
%    figure;subplot(2,1,1);plot(neuron_spikes_m(inhibited_m(1,i),:));subplot(2,1,2);cusum_(neuron_spikes_m_zscore(inhibited_m(1,i),:),5,1,mean(neuron_spikes_m_zscore(inhibited_m(1,i),1:10)),1,"all");
% end


