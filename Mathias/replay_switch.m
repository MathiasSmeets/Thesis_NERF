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

switch_data_m = load(fullfile(volume_base2, path_to_data, "after_stimulus_data_m_horridge.mat"));
switch_data_m = switch_data_m.after_stimulus_data_m;

output_m = load(fullfile(volume_base2, path_to_noi, "neurons_of_interest_horridge_m.mat"));
output_m = output_m.output_m;

inhibited_m = load(fullfile(volume_base2, path_to_noi, "inhibited_horridge_m.mat"));
inhibited_m = inhibited_m.inhibited_m;

template = load(fullfile(volume_base2, path_to_clusters, "template_m.mat"));template = template.template;
template_cluster = load(fullfile(volume_base2, path_to_clusters, "template_cluster_m.mat"));template_cluster = template_cluster.template_cluster;
%% 

interval_size = 70;
intervals_together = 30;
bins_together = 15;
intervals_together_before = intervals_together*interval_size;
last_interval_data = zeros(1,size(stimulus_data_m,1));

% cur_activity_before = zeros(size(stimulus_data_m,1),size(before_data_m,2)-1-size(cur_template,2)+1);
% cur_activity_after = zeros(size(stimulus_data_m,1),size(after_data_m,2)-1-size(cur_template,2)+1);
% cur_activity_between = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_horridge = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_switch = cell(size(stimulus_data_m,1),1);
for i = 1:size(stimulus_data_m,1)
    % get last interval
    % get last interval
    counter = size(stimulus_data_m,2);
    last_interval_data(i) = counter;
    while isempty(stimulus_data_m{i,counter})
        counter = counter - 1;
        last_interval_data(i) = counter;
    end

    % do a compression of 0.2; 0.5; 3; 5
    compressed_02 = create_compressed_template(template{i},0.1);
    compressed_05 = create_compressed_template(template{i},0.5);
    compressed_3 = create_compressed_template(template{i},3);
    compressed_5 = create_compressed_template(template{i},5);


    cur_template = template{i};
    cur_cluster = template_cluster{i};

    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-size(cur_template,2)+1
            if sum(sum(stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1), 2) > 0) >= 2
                adj_cur_correlation_horridge{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(cur_template.*stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
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

%% calculate average
% compression 1
avg_cur_cor_horrige = zeros(numel(adj_cur_correlation_horridge),1);
for i = 1:numel(avg_cur_cor_horrige)
    avg_cur_cor_horrige(i) = mean(adj_cur_correlation_horridge{i});
end
avg_cur_cor_switch = zeros(numel(adj_cur_correlation_switch),1);
for i = 1:numel(avg_cur_cor_switch)
    avg_cur_cor_switch(i) = mean(adj_cur_correlation_switch{i});
end

% create boxplots
figure
boxplot([avg_cur_cor_horrige, avg_cur_cor_switch], 'Labels', {'Horridge', 'Switch'})
hold on
scatter(ones(size(avg_cur_cor_horrige,1)),avg_cur_cor_horrige, 'filled', 'blue')
scatter(ones(size(avg_cur_cor_switch,1))*2,avg_cur_cor_switch, 'filled', 'blue')
line([ones(size(avg_cur_cor_horrige)), ones(size(avg_cur_cor_switch))*2]',[avg_cur_cor_horrige, avg_cur_cor_switch]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_horridge_switch.png")
