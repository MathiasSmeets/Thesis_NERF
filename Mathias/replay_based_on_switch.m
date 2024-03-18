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

switch_data_m = load(fullfile(volume_base2, path_to_data, "after_stimulus_data_m_switch.mat"));
switch_data_m = switch_data_m.after_stimulus_switch_m;

stimulus_data_m = load(fullfile(volume_base2, path_to_data, "after_stimulus_data_m_horridge.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;

output_m = load(fullfile(volume_base2, path_to_noi, "neurons_of_interest_horridge_m.mat"));
output_m = output_m.output_m;

inhibited_m = load(fullfile(volume_base2, path_to_noi, "inhibited_horridge_m.mat"));
inhibited_m = inhibited_m.inhibited_m;

after_data_m = load(fullfile(volume_base2, path_to_data, "waiting_data_m.mat"));
after_data_m = after_data_m.waiting_data;

before_data_m = load(fullfile(volume_base2, path_to_data, "before_data_m.mat"));
before_data_m = before_data_m.before_data;

ica_assemblies = load(fullfile(volume_base2,path_to_clusters, "assemblies_switch_m.mat")); ica_assemblies = ica_assemblies.total_assemblies;
ica_data = load(fullfile(volume_base2,path_to_clusters, "data_switch_m.mat")); ica_data = ica_data.total_data;
ica_neurons_of_interest = load(fullfile(volume_base2,path_to_clusters, "neurons_of_interest_switch_m.mat")); ica_neurons_of_interest = ica_neurons_of_interest.total_neurons_of_interest;
ica_activity = load(fullfile(volume_base2,path_to_clusters, "activity_switch_m.mat")); ica_activity = ica_activity.total_activity;
ica_vector = load(fullfile(volume_base2,path_to_clusters, "ica_vector_switch_m.mat")); ica_vector = ica_vector.total_vector;

ica_neurons_of_interest_before = load(fullfile(volume_base2, path_to_clusters, "neurons_of_interest_after_m.mat"));ica_neurons_of_interest_before = ica_neurons_of_interest_before.total_neurons_of_interest;
ica_assemblies_before = load(fullfile(volume_base2, path_to_clusters, "assemblies_after_m.mat")); ica_assemblies_before = ica_assemblies_before.total_assemblies;
%% get template

interval_size = 70;
intervals_together = 30;
bins_together = 15;
intervals_together_before = intervals_together*interval_size;
last_interval_data = zeros(1,size(switch_data_m,1));
template = cell(1,size(switch_data_m,1));
template_cluster = cell(1,size(switch_data_m,1));
template_vector = cell(1,size(switch_data_m,1));
template_cluster_count = cell(1,size(switch_data_m,1));
for i = 1:size(switch_data_m,1)
    % get last interval
    % get last interval
    counter = size(switch_data_m,2);
    last_interval_data(i) = counter;
    while isempty(switch_data_m{i,counter})
        counter = counter - 1;
        last_interval_data(i) = counter;
    end
    cur_avg = zeros(size(switch_data_m{i,last_interval_data(i)}));
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
    threshold_20percent = 0.2*last_interval_before;
    all_active_assemblies_before_before_threshold{i} = all_assemblies_before;
    all_act_as_count{i} = all_assemblies_count_before;
    active_assemblies_before = all_assemblies_before(all_assemblies_count_before>threshold_20percent);
    all_active_assemblies_before{i} = active_assemblies_before;
    
    % get most active assembly that is not active before experiment (20% of the time)
    maxvalue = sort(all_assemblies_count,'descend');
    value_counter = 1;
    cur_value = maxvalue(value_counter);
    most_common_cluster = all_assemblies(all_assemblies_count==cur_value);
    most_common_cluster = most_common_cluster{end};
    while ~isempty(find(cellfun(@(x) isequal(x, most_common_cluster), active_assemblies_before), 1))
        value_counter = value_counter+1;
        cur_value = maxvalue(value_counter);
        most_common_cluster = all_assemblies(all_assemblies_count==cur_value);
        if cur_value == maxvalue(value_counter-1)
            most_common_cluster = most_common_cluster{end-1};
        else
            most_common_cluster = most_common_cluster{end};
        end
    end

    template_cluster{i} = most_common_cluster;
    template_cluster_count{i} = cur_value;
    % in order to get vector, take average vector but fix sign ambiguity first
    %template_vector{i} = all_vectors{maxindex}(:,1);

    %% go to each time this cluster is active, find peaks in activity and get candidate templates
    cur_template = zeros(numel(most_common_cluster),bins_together);
    counter = 0;
    for j = 1:cur_last_interval
        for k = 1:numel(ica_assemblies{i,j})
            if isequal(ica_neurons_of_interest{i,j}(ica_assemblies{i,j}{k}), most_common_cluster)
                if length(ica_activity{i,j}(:,k)) >= 3
                    [pks, locs] = findpeaks(abs(ica_activity{i,j}(:,k)),"NPeaks",intervals_together,"MinPeakHeight",0.4*max(abs(ica_activity{i,j}(:,k))));
                    for jj = 1:length(locs)
                        raw_data_index = (j-1)*intervals_together + ceil((locs(jj))/ceil(interval_size/bins_together));
                        position_in_data = mod(locs(jj)-1,ceil(interval_size/bins_together))*bins_together+1;
                        cur_raw_data = stimulus_data_m{i,raw_data_index};
                        cur_assembly_data = cur_raw_data(template_cluster{i}, position_in_data:min(position_in_data+bins_together-1, interval_size));
                        cur_template(:,1:size(cur_assembly_data,2)) = cur_template(:,1:size(cur_assembly_data,2)) + cur_assembly_data;
                        counter = counter + 1;
                    end
                end
            end
        end
    end
    template{i} =  cur_template/counter;
end

save('/scratch/mathiass-takeokalab/01/template_m.mat', 'template')
save('/scratch/mathiass-takeokalab/01/template_cluster_m.mat', 'template_cluster')
%% check for replay in before and after data

adj_cur_correlation_before = zeros(size(stimulus_data_m,1),size(before_data_m,2)-1-size(cur_template,2)+1);
adj_cur_correlation_after = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_between = cell(size(stimulus_data_m,1),1);% make a cell because different lengths
adj_cur_correlation_switch = cell(size(switch_data_m,1),1);
for i = 1:size(switch_data_m,1)
    cur_before_data = before_data_m(before_data_m(:,1) == i,:);
    cur_before_data = cur_before_data(:,2:end);
    %cur_after_data = after_data_m(after_data_m(:,1) == i,:);
    %cur_after_data = cur_after_data(:,2:end);
    cur_after_data = after_data_m{i};
    cur_after_data = cell2mat(cur_after_data);

    cur_template = template{i};
    cur_cluster = template_cluster{i};

    % check if at least 2 spikes in data, otherwise 0
    for j = 1:size(cur_before_data,2)-size(cur_template,2)+1
        if sum(sum(cur_before_data(cur_cluster,j:j+size(cur_template,2)-1), 2) > 0) >= 2
            adj_cur_correlation_before(i,j) = sum(cur_template.*cur_before_data(cur_cluster,j:j+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
        end
    end
    adj_cur_correlation_after{i} = zeros(1,size(cur_after_data,2)-1-size(cur_template,2)+1);
    for j = 1:size(cur_after_data,2)-size(cur_template,2)+1
        if sum(sum(cur_after_data(cur_cluster,j:j+size(cur_template,2)-1), 2) > 0) >= 2
            adj_cur_correlation_after{i}(j) = sum(cur_template.*cur_after_data(cur_cluster,j:j+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
        end
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(switch_data_m{i,j},2)-size(cur_template,2)+1
            if sum(sum(switch_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1), 2) > 0) >= 2
                adj_cur_correlation_between{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(cur_template.*stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
            end
        end
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(switch_data_m{i,j},2)-size(cur_template,2)+1
            if sum(sum(switch_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1), 2) > 0) >= 2
                adj_cur_correlation_switch{i}((j-1)*size(switch_data_m{i,j},2)+k) = sum(cur_template.*switch_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
            end
        end
    end
end
% using activity is a worse measure, as it is only one column, so using
% this is actually a template of one column only


avg_adj_cur_correlation_before = mean(adj_cur_correlation_before,2);
avg_adj_cur_correlation_between = zeros(numel(adj_cur_correlation_between),1);
avg_adj_cur_correlation_after = zeros(numel(adj_cur_correlation_after),1);
avg_adj_cur_correlation_switch = zeros(numel(adj_cur_correlation_switch),1);
for i = 1:numel(adj_cur_correlation_between)
    avg_adj_cur_correlation_between(i) = mean(adj_cur_correlation_between{i});
    avg_adj_cur_correlation_after(i) = mean(adj_cur_correlation_after{i});
    avg_adj_cur_correlation_switch(i) = mean(adj_cur_correlation_switch{i});
end

figure
boxplot([avg_adj_cur_correlation_before, avg_adj_cur_correlation_between, avg_adj_cur_correlation_after,avg_adj_cur_correlation_switch], 'Labels', {'Baseline', 'Experiment', 'Rest','Switch'})
hold on
scatter(ones(size(avg_adj_cur_correlation_before,1)),avg_adj_cur_correlation_before, 'filled', 'blue')
scatter(ones(size(avg_adj_cur_correlation_between,1))*2,avg_adj_cur_correlation_between, 'filled', 'blue')
scatter(ones(size(avg_adj_cur_correlation_after,1))*3,avg_adj_cur_correlation_after, 'filled', 'blue')
scatter(ones(size(avg_adj_cur_correlation_switch,1))*4,avg_adj_cur_correlation_switch,'filled', 'blue')
line([ones(size(avg_adj_cur_correlation_before)), ones(size(avg_adj_cur_correlation_between))*2]',[avg_adj_cur_correlation_before, avg_adj_cur_correlation_between]','Color','green')
line([ones(size(avg_adj_cur_correlation_between))*2, ones(size(avg_adj_cur_correlation_after))*3]',[avg_adj_cur_correlation_between, avg_adj_cur_correlation_after]','Color','green')
line([ones(size(avg_adj_cur_correlation_after))*3, ones(size(avg_adj_cur_correlation_switch))*3]',[avg_adj_cur_correlation_after, avg_adj_cur_correlation_switch]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_adjusted_bbas_based_on_switch.png")

%% determine 95% threshold based on before data ---->>> not working well

% wilcoxin signed-rank test (no gaussian assumptions +  paired data)
% [p,h,stats] = signrank(x,y)

% wilcoxon rank-sum test (no gaussian assumptions + no paired data)
% [p,h,stats] = ranksum(x,y)

