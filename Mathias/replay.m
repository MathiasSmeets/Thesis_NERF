if ispc
    volume_base = '\\nerffs13';
    volume_base2 = '\\nerffs17';
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/';
    volume_base2 = '/mnt/takeokalab/';
end

path_to_data = "takeokalabwip2023/Mathias/switch_data/data_after_stimulus";
path_to_clusters = "takeokalabwip2023/Mathias/switch_data/clusters";
path_to_noi = "takeokalabwip2023/Mathias/switch_data/neurons_of_interest";

stimulus_data_m = load(fullfile(volume_base2, path_to_data, "after_stimulus_data_m_horridge.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
%stimulus_data_m = stimulus_data_m(1:9,:);
stimulus_data_m = stimulus_data_m(1:9,:);

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

mouse_to_exclude = 0;
%mouse_to_exclude = 2;
%% get template

interval_size = 60;
interval_size_before = 70;
intervals_together = 30;
bins_together = 15;
intervals_together_before = intervals_together*interval_size_before;
last_interval_data = zeros(1,size(stimulus_data_m,1));
template = cell(1,size(stimulus_data_m,1));
template_smoothed = cell(1,size(stimulus_data_m,1));
template_cluster = cell(1,size(stimulus_data_m,1));
template_vector = cell(1,size(stimulus_data_m,1));
template_cluster_count = cell(1,size(stimulus_data_m,1));
for i = setdiff(1:size(stimulus_data_m,1),mouse_to_exclude)
    % get last interval
    % get last interval
    counter = size(stimulus_data_m,2);
    last_interval_data(i) = counter;
    while isempty(stimulus_data_m{i,counter})
        counter = counter - 1;
        last_interval_data(i) = counter;
    end
    cur_avg = zeros(size(stimulus_data_m{i,last_interval_data(i)}));
    %% template: avergae neurons in cluster that is most common
    % find cluster that is most common
    cur_last_interval = ceil(last_interval_data(i)/intervals_together);
    all_assemblies = {};
    all_assemblies_count = [];
    all_vectors = {};
    for j = 1:cur_last_interval
        for k = 1:numel(ica_assemblies{i,j})
            cur_assembly = ica_neurons_of_interest{i,j}(ica_assemblies{i,j}{k});
            cur_vector = ica_vector{i,j}{k};
            idx = find(cellfun(@(x) isequal(x, cur_assembly), all_assemblies));
            if ~isempty(idx)
                all_assemblies_count(idx) = all_assemblies_count(idx) + 1;
                all_vectors{idx} = [all_vectors{idx}, cur_vector];
            else
                all_assemblies{end+1} = cur_assembly;
                all_assemblies_count = [all_assemblies_count, 1];
                all_vectors{end + 1} = cur_vector;
            end
        end
    end
    % check if cluster that is found is not also most occurring in before data, if this is the case, take next one that does not contain any of these neurons
    cur_before_data = before_data_m(before_data_m(:,1) == i,:);
    cur_before_data = cur_before_data(:,2:end);
    last_interval_before = ceil(size(cur_before_data,2)/intervals_together_before);
    all_assemblies_before = {};
    all_assemblies_count_before = [];
    for j = 1:last_interval_before
        for k = 1:numel(ica_assemblies_before{i,j})
            cur_assembly_before = ica_neurons_of_interest_before{i,j}(ica_assemblies_before{i,j}{k});
            idx_before = find(cellfun(@(x) isequal(x, cur_assembly_before), all_assemblies_before));
            if ~isempty(idx_before)
                all_assemblies_count_before(idx_before) = all_assemblies_count_before(idx_before) + 1;
            else
                all_assemblies_before{end+1} = cur_assembly_before;
                all_assemblies_count_before = [all_assemblies_count_before, 1];
            end
        end
    end
    occurrences = zeros(1,size(cur_before_data,1));
    % get occurence in cluster for each neuron
    for j = 1:numel(all_assemblies_before)
        occurrences = occurrences + histcounts(all_assemblies_before{j}, [1:size(cur_before_data,1)+1]-0.5) * all_assemblies_count_before(j);
    end

    % if neurons are active in clusters similar to between, avoid these
    % similar = at least activity during experiment - 10%
    neuron_between_activities = sum(cluster_matrices_between_m{i},2)/size(cluster_matrices_between_m{i},2); % in percentage
    variable_threshold = 0.5*neuron_between_activities*last_interval_before;
    neurons_to_avoid = find(occurrences>variable_threshold');

    %active_assemblies_before = all_assemblies_before(all_assemblies_count_before>threshold_10percent);
    %all_active_assemblies_before{i} = active_assemblies_before;
    
    % 
    maxvalue = sort(all_assemblies_count,'descend');
    value_counter = 1;
    cur_value = maxvalue(value_counter);
    most_common_cluster = all_assemblies(all_assemblies_count==cur_value);
    actual_most_common_cluster = most_common_cluster{end};
    internal_counter = 0;
    % contains neuron to avoid
    while numel(find(ismember(actual_most_common_cluster, neurons_to_avoid))) >= 1
        % if another has equal count
        if ~isequal(actual_most_common_cluster, most_common_cluster{1})
            internal_counter = internal_counter + 1;
            actual_most_common_cluster = most_common_cluster{end-internal_counter};
        % else go to next count
        else
            value_counter = value_counter + 1;
            cur_value = maxvalue(value_counter);
            most_common_cluster = all_assemblies(all_assemblies_count==cur_value);
            actual_most_common_cluster = most_common_cluster{end};
            internal_counter = 0;
        end
        % if count is one and we went over all, just take the most common anyways
        if cur_value == 1 && isequal(actual_most_common_cluster, most_common_cluster{1})
            cur_value = maxvalue(1);
            most_common_cluster = all_assemblies(all_assemblies_count==cur_value);
            actual_most_common_cluster = most_common_cluster{end};
            break;
        end
    end

    % while ~isempty(find(cellfun(@(x) isequal(x, most_common_cluster), active_assemblies_before), 1))
    %     value_counter = value_counter+1;
    %     cur_value = maxvalue(value_counter);
    %     most_common_cluster = all_assemblies(all_assemblies_count==cur_value);
    %     if cur_value == maxvalue(value_counter-1)
    %         most_common_cluster = most_common_cluster{end-1};
    % 
    %     else
    %         most_common_cluster = most_common_cluster{end};
    %     end
    % end

    template_cluster{i} = actual_most_common_cluster;
    template_cluster_count{i} = cur_value;
    % in order to get vector, take average vector but fix sign ambiguity first
    maxindex = find(cellfun(@(x) isequal(x, actual_most_common_cluster), all_assemblies));
    template_vector{i} = all_vectors{maxindex}(:,1);

    %% go to each time this cluster is active, find peaks in activity and get candidate templates
    cur_template = zeros(numel(actual_most_common_cluster),bins_together);
    counter = 0;
    for j = 1:cur_last_interval
        for k = 1:numel(ica_assemblies{i,j})
            if isequal(ica_neurons_of_interest{i,j}(ica_assemblies{i,j}{k}), actual_most_common_cluster)
                if length(ica_activity{i,j}(:,k)) >= 3
                    [pks, locs] = findpeaks(abs(ica_activity{i,j}(:,k)),"NPeaks",intervals_together,"MinPeakHeight",0.4*max(abs(ica_activity{i,j}(:,k))));
                    for jj = 1:length(locs)
                        raw_data_index = (j-1)*intervals_together + ceil((locs(jj))/ceil(interval_size/bins_together));
                        position_in_data = mod(locs(jj)-1,ceil(interval_size/bins_together))*bins_together+1;
                        cur_raw_data = stimulus_data_m{i,raw_data_index};
                        cur_assembly_data = cur_raw_data(template_cluster{i}, position_in_data:min(position_in_data+bins_together-1, interval_size));
                        cur_template(:,1:size(cur_assembly_data,2)) = cur_template(:,1:size(cur_assembly_data,2)) + double(cur_assembly_data);
                        counter = counter + 1;
                    end
                end
            end
        end
    end
    cur_template = cur_template/counter;
    template{i} =  cur_template;


    % Gaussian kernel parameters
    kernel_size = 3; % Size of the kernel (odd number)
    sigma = 1; % Standard deviation of the Gaussian kernel

    % Create the Gaussian kernel
    kernel = gausswin(kernel_size, sigma);

    % Normalize the kernel
    kernel = kernel / sum(kernel);


    % Initialize a matrix to store smoothed data
    smoothed_template = zeros(size(cur_template));

    % Apply convolution to each row independently
    for j = 1:size(cur_template, 1)
        smoothed_template(j, :) = conv(cur_template(j, :), kernel, 'same');
    end
    template_smoothed{i} = flip(smoothed_template,2);
end

save('/scratch/mathiass-takeokalab/01/template_m_reverse_template.mat', 'template')
save('/scratch/mathiass-takeokalab/01/template_cluster_m_reverse_template.mat', 'template_cluster')
save('/scratch/mathiass-takeokalab/01/template_smoothed_3_m_reverse_template.mat', 'template_smoothed')
%% check for replay in before and after data
cur_correlation_before = zeros(size(stimulus_data_m,1),size(before_data_m,2)-1-size(cur_template,2)+1);
cur_correlation_after = cell(size(stimulus_data_m,1),1);
cur_correlation_between = cell(size(stimulus_data_m,1),1);% make a cell because different lengths
cur_activity_before = zeros(size(stimulus_data_m,1),size(before_data_m,2)-1-size(cur_template,2)+1);
cur_activity_after = zeros(size(stimulus_data_m,1),size(after_data_m,2)-1-size(cur_template,2)+1);
cur_activity_between = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_before = zeros(size(stimulus_data_m,1),size(before_data_m,2)-1-size(cur_template,2)+1);
adj_cur_correlation_after = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_between = cell(size(stimulus_data_m,1),1);% make a cell because different lengths
for i = setdiff(1:size(stimulus_data_m,1),mouse_to_exclude)
    cur_before_data = before_data_m(before_data_m(:,1) == i,:);
    cur_before_data = cur_before_data(:,2:end);
    %cur_after_data = after_data_m(after_data_m(:,1) == i,:);
    %cur_after_data = cur_after_data(:,2:end);
    cur_after_data = after_data_m{i};
    cur_after_data = cell2mat(cur_after_data);

    cur_template = template_smoothed{i};
    %cur_template = template{i};
    cur_cluster = template_cluster{i};

    for j = 1:size(cur_before_data,2)-size(cur_template,2)+1
        cur_correlation_before(i,j) = sum(cur_template.*cur_before_data(cur_cluster,j:j+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
    end
    cur_correlation_after{i} = zeros(1,size(after_data_m,2)-1-size(cur_template,2)+1);
    for j = 1:size(cur_after_data,2)-size(cur_template,2)+1
        cur_correlation_after{i}(j) = sum(cur_template.*cur_after_data(cur_cluster,j:j+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-size(cur_template,2)+1
            cur_correlation_between{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(cur_template.*double(stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1)),'all') / (size(cur_template,1) * size(cur_template,2));
        end
    end

    % check if at least 2 spikes in data, otherwise 0
    for j = 1:size(cur_before_data,2)-size(cur_template,2)+1
        if sum(sum(cur_before_data(cur_cluster,j:j+size(cur_template,2)-1), 2) > 0) >= 2 || numel(cur_cluster) == 2
            adj_cur_correlation_before(i,j) = sum(cur_template.*cur_before_data(cur_cluster,j:j+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
        end
    end
    adj_cur_correlation_after{i} = zeros(1,size(cur_after_data,2)-1-size(cur_template,2)+1);
    for j = 1:size(cur_after_data,2)-size(cur_template,2)+1
        if sum(sum(cur_after_data(cur_cluster,j:j+size(cur_template,2)-1), 2) > 0) >= 2 || numel(cur_cluster) == 2
            adj_cur_correlation_after{i}(j) = sum(cur_template.*cur_after_data(cur_cluster,j:j+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
        end
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-size(cur_template,2)+1
            if sum(sum(stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1), 2) > 0) >= 2 || numel(cur_cluster) == 2
                adj_cur_correlation_between{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(cur_template.*double(stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1)),'all') / (size(cur_template,1) * size(cur_template,2));
            end
        end
    end
    for j = 1:size(cur_before_data,2)-15
        cur_activity_before(i,j) = sum(cur_before_data(cur_cluster, j:j+14),2)'*template_vector{i};
    end
    for j = 1:size(cur_after_data,2)-15
        cur_activity_after(i,j) = sum(cur_after_data(cur_cluster,j:j+14),2)'*template_vector{i};
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-15
            cur_activity_between{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(double(stimulus_data_m{i,j}(cur_cluster,k:k+14)),2)'*template_vector{i};
        end
    end
end
% using activity is a worse measure, as it is only one column, so using
% this is actually a template of one column only
% but now i tried with 15 ms bins which makes much more sense so let's see
%%
threshold = zeros(numel(cur_correlation_after),1);
for i = setdiff(1:size(stimulus_data_m,1),mouse_to_exclude)
    threshold(i) = prctile(cur_correlation_after{i},99,2);
end

thresholded_before = cur_correlation_before;
thresholded_between = cur_correlation_between;
thresholded_after = cur_correlation_after;

thresholded_before(thresholded_before<threshold) = 0;
for i = setdiff(1:size(stimulus_data_m,1),mouse_to_exclude)
    thresholded_between{i}(thresholded_between{i}<threshold(i)) = 0;
    thresholded_after{i}(thresholded_after{i}<threshold(i)) = 0;
end

avg_cor_before = mean(thresholded_before,2);
avg_cor_between=zeros(numel(thresholded_between),1);
avg_cor_after=zeros(numel(thresholded_after),1);
avg_act_between=zeros(numel(cur_activity_between),1);
for i = setdiff(1:size(stimulus_data_m,1),mouse_to_exclude)
    avg_cor_between(i) = mean(thresholded_between{i});
    avg_cor_after(i) = mean(thresholded_after{i});
    avg_act_between(i) = mean(cur_activity_between{i});
end
avg_act_before = mean(abs(cur_activity_before),2);
avg_act_after = mean(abs(cur_activity_after),2);

%% create boxplot figures

figure;boxplot([avg_cor_before,avg_cor_between,avg_cor_after])
figure
boxplot([avg_cor_before, avg_cor_between, avg_cor_after], 'Labels', {'Baseline', 'Experiment', 'Rest'})
hold on
scatter(ones(size(avg_cor_before,1)),avg_cor_before, 'filled', 'blue')
scatter(ones(size(avg_cor_between,1))*2,avg_cor_between, 'filled', 'blue')
scatter(ones(size(avg_cor_after,1))*3,avg_cor_after, 'filled', 'blue')
line([ones(size(avg_cor_before)), ones(size(avg_cor_between))*2]',[avg_cor_before, avg_cor_between]','Color','green')
line([ones(size(avg_cor_between))*2, ones(size(avg_cor_after))*3]',[avg_cor_between, avg_cor_after]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_bba_reverse_template.png")

figure
boxplot([avg_cor_before, avg_cor_after], 'Labels', {'Baseline', 'Rest'})
hold on
scatter(ones(size(avg_cor_before,1)),avg_cor_before, 'filled', 'blue')
scatter(ones(size(avg_cor_after,1))*2,avg_cor_after, 'filled', 'blue')
line([ones(size(avg_cor_before)), ones(size(avg_cor_after))*2]',[avg_cor_before, avg_cor_after]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_ba_reverse_template.png")

avg_adj_cur_correlation_before = mean(adj_cur_correlation_before,2);
avg_adj_cur_correlation_between = zeros(numel(adj_cur_correlation_between),1);
avg_adj_cur_correlation_after = zeros(numel(adj_cur_correlation_after),1);
for i = setdiff(1:size(stimulus_data_m,1),mouse_to_exclude)
    avg_adj_cur_correlation_between(i) = mean(adj_cur_correlation_between{i});
    avg_adj_cur_correlation_after(i) = mean(adj_cur_correlation_after{i});
end

figure;boxplot([avg_adj_cur_correlation_before,avg_adj_cur_correlation_between,avg_adj_cur_correlation_after])
figure
boxplot([avg_adj_cur_correlation_before, avg_adj_cur_correlation_between, avg_adj_cur_correlation_after], 'Labels', {'Baseline', 'Experiment', 'Rest'})
hold on
scatter(ones(size(avg_adj_cur_correlation_before,1)),avg_adj_cur_correlation_before, 'filled', 'blue')
scatter(ones(size(avg_adj_cur_correlation_between,1))*2,avg_adj_cur_correlation_between, 'filled', 'blue')
scatter(ones(size(avg_adj_cur_correlation_after,1))*3,avg_adj_cur_correlation_after, 'filled', 'blue')
line([ones(size(avg_adj_cur_correlation_before)), ones(size(avg_adj_cur_correlation_between))*2]',[avg_adj_cur_correlation_before, avg_adj_cur_correlation_between]','Color','green')
line([ones(size(avg_adj_cur_correlation_between))*2, ones(size(avg_adj_cur_correlation_after))*3]',[avg_adj_cur_correlation_between, avg_adj_cur_correlation_after]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_adjusted_bba_reverse_template.png")

figure
boxplot([avg_adj_cur_correlation_before, avg_adj_cur_correlation_after], 'Labels', {'Baseline', 'Rest'})
hold on
scatter(ones(size(avg_adj_cur_correlation_before,1)),avg_adj_cur_correlation_before, 'filled', 'blue')
scatter(ones(size(avg_adj_cur_correlation_after,1))*2,avg_adj_cur_correlation_after, 'filled', 'blue')
line([ones(size(avg_adj_cur_correlation_before)), ones(size(avg_adj_cur_correlation_after))*2]',[avg_adj_cur_correlation_before, avg_adj_cur_correlation_after]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_adjusted_ba_reverse_template.png")

save("/scratch/mathiass-takeokalab/01/correlation_before_smoothed_width3_reverse_template.mat","avg_adj_cur_correlation_before")
save("/scratch/mathiass-takeokalab/01/correlation_between_smoothed_width3_reverse_template.mat","avg_adj_cur_correlation_between")
save("/scratch/mathiass-takeokalab/01/correlation_after_smoothed_width3_reverse_template.mat","avg_adj_cur_correlation_after")
%% determine 95% threshold based on before data ---->>> not working well
% threshold = prctile(cur_correlation_before,99,2);
% adjusted_cor_before = cur_correlation_before;
% adjusted_cor_between = cur_correlation_between;
% adjusted_cor_after= cur_correlation_after;
% adjusted_cor_before(adjusted_cor_before<threshold) = 0;
% for i = setdiff(1:size(stimulus_data_m,1),exclude)
%     adjusted_cor_between{i}(adjusted_cor_between{i}<threshold(i)) = 0;
% end
% adjusted_cor_after(adjusted_cor_after<threshold) = 0;
% 
% avg_adj_cor_before = mean(adjusted_cor_before,2);
% avg_adj_cor_between=zeros(numel(adjusted_cor_between),1);
% for i = setdiff(1:size(stimulus_data_m,1),exclude)
%     avg_adj_cor_between(i) = mean(adjusted_cor_between{i});
% end
% avg_adj_cor_after = mean(adjusted_cor_after,2);




%% plot correlations
% wilcoxin signed-rank test (no gaussian assumptions +  paired data)
% [p,h,stats] = signrank(avg_adj_cur_correlation_before,avg_adj_cur_correlation_after)

% wilcoxon rank-sum test (no gaussian assumptions + no paired data)
% [p,h,stats] = ranksum(avg_adj_cur_correlation_before,avg_adj_cur_correlation_after)

% test for gaussianity: shapiro wilk, kolmogorov smimov, qqplot


% for i = setdiff(1:size(stimulus_data_m,1),exclude)
% figure
% subplot(3,1,1);plot(cur_correlation_before(i,:))
% subplot(3,1,2);plot(cur_correlation_between{i})
% subplot(3,1,3);plot(cur_correlation_after{i})
% saveas(gcf,"/scratch/mathiass-takeokalab/01/correlation" + i + ".png")
% figure
% subplot(3,1,1);plot(adj_cur_correlation_before(i,:))
% subplot(3,1,2);plot(adj_cur_correlation_between{i})
% subplot(3,1,3);plot(adj_cur_correlation_after{i})
% saveas(gcf,"/scratch/mathiass-takeokalab/01/correlation_adj" + i + ".png")
% end

%% more neurons in cluster afterwards?
% 
% load("X:\Mathias\switch_data\clusters\cluster_matrices_before_m.mat")
% all_cluster_matrices_before = all_cluster_matrices;
% load("X:\Mathias\switch_data\clusters\cluster_matrices_between_m.mat")
% all_cluster_matrices_between = all_cluster_matrices;
% load("X:\Mathias\switch_data\clusters\cluster_matrices_after_m.mat")
% all_cluster_matrices_after = all_cluster_matrices;
% 
% avg_nb_cluster_before = zeros(numel(all_cluster_matrices_before),1);
% avg_nb_cluster_between = zeros(numel(all_cluster_matrices_between),1);
% avg_nb_cluster_after = zeros(numel(all_cluster_matrices_after),1);
% for i = setdiff(1:size(stimulus_data_m,1),exclude)
%     avg_nb_cluster_before(i) = sum(all_cluster_matrices_before{i},'all') / size(all_cluster_matrices_before{i},2);
%     avg_nb_cluster_between(i) = sum(all_cluster_matrices_between{i},'all') / size(all_cluster_matrices_between{i},2);
%     avg_nb_cluster_after(i) = sum(all_cluster_matrices_after{i},'all') / size(all_cluster_matrices_after{i},2);
% end
% figure
% boxplot([avg_nb_cluster_before, avg_nb_cluster_between, avg_nb_cluster_after], 'Labels', {'Baseline', 'Experiment', 'Rest'})
% hold on
% scatter(ones(size(avg_nb_cluster_before,1)),avg_nb_cluster_before, 'filled', 'blue')
% scatter(ones(size(avg_nb_cluster_between,1))*2,avg_nb_cluster_between, 'filled', 'blue')
% scatter(ones(size(avg_nb_cluster_after,1))*3,avg_nb_cluster_after, 'filled', 'blue')
% line([ones(size(avg_nb_cluster_before)), ones(size(avg_nb_cluster_between))*2]',[avg_nb_cluster_before, avg_nb_cluster_between]','Color','green')
% line([ones(size(avg_nb_cluster_between))*2, ones(size(avg_nb_cluster_after))*3]',[avg_nb_cluster_between, avg_nb_cluster_after]','Color','green')
% 
% 

