% calculate threshold based on random permutated data for each recording
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


flat_stimulus_data_m = load(fullfile(volume_base2, path_to_raw, "horridge_data_m.mat"));
flat_stimulus_data_m = flat_stimulus_data_m.data;
flat_hor_data_m = load(fullfile(volume_base2, path_to_raw, "switch_data_m.mat"));
flat_hor_data_m = flat_hor_data_m.switch_data;

stimulus_data_m = load(fullfile(volume_base2, path_to_data, "after_stimulus_data_m_horridge.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_m = stimulus_data_m(1:9,:);
%stimulus_data_m = stimulus_data_m(1:9,:);

horridge_data_m = load(fullfile(volume_base2, path_to_data, "after_stimulus_data_m_switch.mat"));
horridge_data_m = horridge_data_m.after_stimulus_switch_m;
horridge_data_m = horridge_data_m(1:9,:);

output_m = load(fullfile(volume_base2, path_to_noi, "neurons_of_interest_horridge_m.mat"));
output_m = output_m.output_m;

inhibited_m = load(fullfile(volume_base2, path_to_noi, "inhibited_horridge_m.mat"));
inhibited_m = inhibited_m.inhibited_m;

after_data_m = load(fullfile(volume_base2, path_to_data, "waiting_data_m.mat"));
after_data_m = after_data_m.waiting_data;
after_data_m = after_data_m(1:9,:);

before_data_m = load(fullfile(volume_base2, path_to_data, "before_data_m.mat"));
before_data_m = before_data_m.before_data;
before_data_m(before_data_m(:,1)>=10,:) = [];

ica_assemblies = load(fullfile(volume_base2,path_to_clusters, "assemblies_horridge_m.mat")); ica_assemblies = ica_assemblies.total_assemblies;
ica_data = load(fullfile(volume_base2,path_to_clusters, "data_horridge_m.mat")); ica_data = ica_data.total_data;
ica_neurons_of_interest = load(fullfile(volume_base2,path_to_clusters, "neurons_of_interest_horridge_m.mat")); ica_neurons_of_interest = ica_neurons_of_interest.total_neurons_of_interest;
ica_activity = load(fullfile(volume_base2,path_to_clusters, "activity_horridge_m.mat")); ica_activity = ica_activity.total_activity;
ica_vector = load(fullfile(volume_base2,path_to_clusters, "ica_vector_horridge_m.mat")); ica_vector = ica_vector.total_vector;

ica_neurons_of_interest_before = load(fullfile(volume_base2, path_to_clusters, "neurons_of_interest_before_m.mat"));ica_neurons_of_interest_before = ica_neurons_of_interest_before.total_neurons_of_interest;
ica_assemblies_before = load(fullfile(volume_base2, path_to_clusters, "assemblies_before_m.mat")); ica_assemblies_before = ica_assemblies_before.total_assemblies;
cluster_matrices_between_m = load(fullfile(volume_base2, path_to_clusters, "cluster_matrices_between_m.mat"));cluster_matrices_between_m = cluster_matrices_between_m.all_cluster_matrices;

mouse_to_exclude = 1:8;      % m

load(fullfile(volume_base2,path_to_correlations,"template_smoothed_3_m_switchbased.mat"));
load(fullfile(volume_base2,path_to_correlations,"template_cluster_m_switchbased.mat"));

%% check for replay in before and after data
last_interval_data = zeros(1,size(stimulus_data_m,1));
last_interval_data_horridge = zeros(1,size(stimulus_data_m,1));
correlation_distribution_before = cell(1,9);
correlation_distribution_after = cell(1,9);
correlation_distribution_between = cell(1,9);
correlation_distribution_horridge = cell(1,9);
for i = setdiff(1:size(stimulus_data_m,1),mouse_to_exclude)
    cur_before_data = before_data_m(before_data_m(:,1) == i,:);
    cur_before_data = cur_before_data(:,2:end);
    %cur_after_data = after_data_m(after_data_m(:,1) == i,:);
    %cur_after_data = cur_after_data(:,2:end);
    cur_after_data = after_data_m{i};
    cur_after_data = cell2mat(cur_after_data);
    cur_stim_data = flat_stimulus_data_m(flat_stimulus_data_m(:,1)==i,:);
    cur_stim_data = cur_stim_data(:,2:end);
    cur_hor_data = flat_hor_data_m(flat_hor_data_m(:,1)==i,:);
    cur_hor_data = cur_hor_data(:,2:end);

    counter = size(stimulus_data_m,2);
    last_interval_data(i) = counter;
    while isempty(stimulus_data_m{i,counter})
        counter = counter - 1;
        last_interval_data(i) = counter;
    end
    counter = size(horridge_data_m,2);
    last_interval_data_horridge(i) = counter;
    while isempty(horridge_data_m{i,counter})
        counter = counter - 1;
        last_interval_data_horridge(i) = counter;
    end

    cur_template = template_smoothed{i};
    %cur_template = template{i};
    cur_cluster = template_cluster{i};

    cur_correlation_before = zeros(1,size(cur_before_data,2)-size(cur_template,2)+1);
    cur_correlation_after = zeros(1,size(cur_after_data,2)-size(cur_template,2)+1);
    cur_correlation_between = zeros(1,size(cur_stim_data,2)-size(cur_template,2)+1);
    cur_correlation_horridge = zeros(1,size(cur_hor_data,2)-size(cur_template,2)+1);

    correlation_distribution_before{i} = zeros(0,2);
    correlation_distribution_after{i} = zeros(0,2);
    correlation_distribution_between{i} = zeros(0,2);
    correlation_distribution_horridge{i} = zeros(0,2);
    for random_iteration = 1:1000
        disp(i + ": " +random_iteration)
        % calculate random data
        random_before = zeros(numel(cur_cluster),size(cur_before_data,2));
        random_after = zeros(numel(cur_cluster),size(cur_after_data,2));
        random_between = zeros(numel(cur_cluster),size(cur_stim_data,2));
        random_horridge = zeros(numel(cur_cluster),size(cur_hor_data,2));
        for row_in_cluster = 1:numel(cur_cluster)
            rand_int_before = randperm(size(cur_before_data,2),1);
            rand_int_after = randperm(size(cur_after_data,2),1);
            rand_int_between = randperm(size(cur_stim_data,2),1);
            rand_int_horridge = randperm(size(cur_hor_data,2),1);
            random_before(row_in_cluster,:) = cur_before_data(cur_cluster(row_in_cluster),[rand_int_before:end,1:rand_int_before-1]);
            random_after(row_in_cluster,:) = cur_after_data(cur_cluster(row_in_cluster),[rand_int_after:end,1:rand_int_after-1]);
            random_between(row_in_cluster,:) = cur_stim_data(cur_cluster(row_in_cluster),[rand_int_between:end,1:rand_int_between-1]);
            random_horridge(row_in_cluster,:) = cur_hor_data(cur_cluster(row_in_cluster),[rand_int_horridge:end,1:rand_int_horridge-1]);
        end
        %error("test random numels vs raw numels")
        % calculate current correlation
        for j = 1:size(random_before,2)-size(cur_template,2)+1
            cur_correlation_before(j) = sum(cur_template.*random_before(:,j:j+size(cur_template,2)-1),'all');% / (size(cur_template,1) * size(cur_template,2));
        end
        for j = 1:size(random_after,2)-size(cur_template,2)+1
            cur_correlation_after(j) = sum(cur_template.*random_after(:,j:j+size(cur_template,2)-1),'all');% / (size(cur_template,1) * size(cur_template,2));
        end
        for j = 1:size(random_between,2)-size(cur_template,2)+1
            cur_correlation_between(j) = sum(cur_template.*random_between(:,j:j+size(cur_template,2)-1),'all');% / (size(cur_template,1) * size(cur_template,2));
        end
        for j = 1:size(random_horridge,2)-size(cur_template,2)+1
            cur_correlation_horridge(j) = sum(cur_template.*random_horridge(:,j:j+size(cur_template,2)-1),'all');% / (size(cur_template,1) * size(cur_template,2));
        end
        % put all unique values in array, along with its frequency counts
        [unique_values_before, ~, freq_idx_before] = unique(cur_correlation_before);
        [unique_values_after, ~, freq_idx_after] = unique(cur_correlation_after);
        [unique_values_between, ~, freq_idx_between] = unique(cur_correlation_between);
        [unique_values_horridge, ~, freq_idx_horridge] = unique(cur_correlation_horridge);
        counts_before = accumarray(freq_idx_before, 1);
        counts_after = accumarray(freq_idx_after, 1);
        counts_between = accumarray(freq_idx_between, 1);
        counts_horridge = accumarray(freq_idx_horridge, 1);

        % put them in total correlation distributions
        for j = 1:numel(unique_values_before)
            idx = find(correlation_distribution_before{i}(:, 1) == unique_values_before(j), 1);
            if isempty(idx)
                correlation_distribution_before{i}(end+1, :) = [unique_values_before(j), counts_before(j)];
            else
                correlation_distribution_before{i}(idx, 2) = correlation_distribution_before{i}(idx, 2) + counts_before(j);
            end
        end
        for j = 1:numel(unique_values_after)
            idx = find(correlation_distribution_after{i}(:, 1) == unique_values_after(j), 1);
            if isempty(idx)
                correlation_distribution_after{i}(end+1, :) = [unique_values_after(j), counts_after(j)];
            else
                correlation_distribution_after{i}(idx, 2) = correlation_distribution_after{i}(idx, 2) + counts_after(j);
            end
        end
        for j = 1:numel(unique_values_between)
            idx = find(correlation_distribution_between{i}(:, 1) == unique_values_between(j), 1);
            if isempty(idx)
                correlation_distribution_between{i}(end+1, :) = [unique_values_between(j), counts_between(j)];
            else
                correlation_distribution_between{i}(idx, 2) = correlation_distribution_between{i}(idx, 2) + counts_between(j);
            end
        end
        for j = 1:numel(unique_values_horridge)
            idx = find(correlation_distribution_horridge{i}(:, 1) == unique_values_horridge(j), 1);
            if isempty(idx)
                correlation_distribution_horridge{i}(end+1, :) = [unique_values_horridge(j), counts_horridge(j)];
            else
                correlation_distribution_horridge{i}(idx, 2) = correlation_distribution_horridge{i}(idx, 2) + counts_horridge(j);
            end
        end
    end
    save("/scratch/mathiass-takeokalab/01/correlation_distribution_before_m_switchbased9.mat","correlation_distribution_before")
    save("/scratch/mathiass-takeokalab/01/correlation_distribution_after_m_switchbased9.mat","correlation_distribution_after")
    save("/scratch/mathiass-takeokalab/01/correlation_distribution_between_m_switchbased9.mat","correlation_distribution_between")
    save("/scratch/mathiass-takeokalab/01/correlation_distribution_horridge_m_switchbased9.mat","correlation_distribution_horridge")
end

save("/scratch/mathiass-takeokalab/01/correlation_distribution_before_m_switchbased9.mat","correlation_distribution_before")
save("/scratch/mathiass-takeokalab/01/correlation_distribution_after_m_switchbased9.mat","correlation_distribution_after")
save("/scratch/mathiass-takeokalab/01/correlation_distribution_between_m_switchbased9.mat","correlation_distribution_between")
save("/scratch/mathiass-takeokalab/01/correlation_distribution_horridge_m_switchbased9.mat","correlation_distribution_horridge")


