stimulus_data_m = load("X:\Mathias\switch_data\data_after_stimulus\after_stimulus_data_m_horridge.mat");
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;

output_m = load("X:\Mathias\switch_data\neurons_of_interest\neurons_of_interest_horridge_m.mat");
output_m = output_m.output_m;

inhibited_m = load("X:\Mathias\switch_data\neurons_of_interest\inhibited_horridge_m.mat");
inhibited_m = inhibited_m.inhibited_m;

after_data_m = load("X:\Mathias\switch_data\data_after_stimulus\after_data_m.mat");
after_data_m = after_data_m.after_data;

before_data_m = load("X:\Mathias\switch_data\data_after_stimulus\before_data_m.mat");
before_data_m = before_data_m.before_data;

ica_assemblies = load("X:\Mathias\switch_data\clusters\assemblies_horridge_m.mat"); ica_assemblies = ica_assemblies.total_assemblies;
ica_data = load("X:\Mathias\switch_data\clusters\data_horridge_m.mat"); ica_data = ica_data.total_data;
ica_neurons_of_interest = load("X:\Mathias\switch_data\clusters\neurons_of_interest_horridge_m.mat"); ica_neurons_of_interest = ica_neurons_of_interest.total_neurons_of_interest;
ica_activity = load("X:\Mathias\switch_data\clusters\activity_horridge_m.mat"); ica_activity = ica_activity.total_activity;

%% get template

interval_size = 70;
intervals_together = 30;
bins_together = 15;
last_interval_data = zeros(1,size(stimulus_data_m,1));
template = cell(1,size(stimulus_data_m,1));
cluster_template = cell(1,size(stimulus_data_m,1));
for i = 1:size(stimulus_data_m,1)
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
    for j = 1:cur_last_interval
        for k = 1:numel(ica_assemblies{i,j})
            cur_assembly = ica_neurons_of_interest{i,j}(ica_assemblies{i,j}{k});
            idx = find(cellfun(@(x) isequal(x, cur_assembly), all_assemblies));
            if ~isempty(idx)
                all_assemblies_count(idx) = all_assemblies_count(idx) + 1;
            else
                all_assemblies{end+1} = cur_assembly;
                all_assemblies_count = [all_assemblies_count, 1];
            end
        end
    end
    maxvalue = max(all_assemblies_count);
    maxindex = find(all_assemblies_count == maxvalue,1, 'last');
    most_common_cluster = all_assemblies{maxindex};

    %% go to each time this cluster is active, find peaks in activity and get candidate templates
    cur_template = zeros(numel(most_common_cluster),bins_together);
    counter = 0;
    for j = 1:cur_last_interval
        for k = 1:numel(ica_assemblies{i,j})
            if isequal(ica_neurons_of_interest{i,j}(ica_assemblies{i,j}{k}), most_common_cluster)
                cluster_template{i} = most_common_cluster;
                if length(ica_activity{i,j}(:,k)) >= 3
                    [pks, locs] = findpeaks(abs(ica_activity{i,j}(:,k)),"NPeaks",intervals_together,"MinPeakHeight",0.4*max(abs(ica_activity{i,j}(:,k))));
                    for jj = 1:length(locs)
                        raw_data_index = (j-1)*intervals_together + ceil((locs(jj))/ceil(interval_size/bins_together));
                        position_in_data = mod(locs(jj)-1,ceil(interval_size/bins_together))*bins_together+1;
                        cur_raw_data = stimulus_data_m{i,raw_data_index};
                        cur_assembly_data = cur_raw_data(cluster_template{i}, position_in_data:min(position_in_data+bins_together-1, interval_size));
                        cur_template(:,1:size(cur_assembly_data,2)) = cur_template(:,1:size(cur_assembly_data,2)) + cur_assembly_data;
                        counter = counter + 1;
                    end
                end
            end
        end
    end
    template{i} =  cur_template/counter;
end


%% check for replay in before and after data
cur_correlation_before = zeros(size(stimulus_data_m,1),size(before_data_m,2)-1-size(cur_template,2)+1);
cur_correlation_after = zeros(size(stimulus_data_m,1),size(before_data_m,2)-1-size(cur_template,2)+1);
cur_correlation_between = cell(size(stimulus_data_m,1),1);% make a cell because different lengths
for i = 1:size(stimulus_data_m,1)
    cur_before_data = before_data_m(before_data_m(:,1) == i,:);
    cur_before_data = cur_before_data(:,2:end);
    cur_after_data = after_data_m(after_data_m(:,1) == i,:);
    cur_after_data = cur_after_data(:,2:end);

    cur_template = template{i};
    cur_cluster = cluster_template{i};

    for j = 1:size(cur_before_data,2)-size(cur_template,2)+1
        cur_correlation_before(i,j) = sum(cur_template.*cur_before_data(cur_cluster,j:j+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
    end
    for j = 1:size(cur_after_data,2)-size(cur_template,2)+1
        cur_correlation_after(i,j) = sum(cur_template.*cur_after_data(cur_cluster,j:j+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-size(cur_template,2)+1
            cur_correlation_between{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(cur_template.*stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
        end
    end
end
avg_cor_before = mean(cur_correlation_before,2);
avg_cor_after = mean(cur_correlation_after,2);
figure;boxplot([avg_cor_before,avg_cor_after])
    
