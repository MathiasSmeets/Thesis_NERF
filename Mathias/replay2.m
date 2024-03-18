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

output_m = load(fullfile(volume_base2, path_to_noi, "neurons_of_interest_horridge_m.mat"));
output_m = output_m.output_m;

inhibited_m = load(fullfile(volume_base2, path_to_noi, "inhibited_horridge_m.mat"));
inhibited_m = inhibited_m.inhibited_m;

after_data_m = load(fullfile(volume_base2, path_to_data, "waiting_data_m.mat"));
after_data_m = after_data_m.waiting_data;

before_data_m = load(fullfile(volume_base2, path_to_data, "before_data_m.mat"));
before_data_m = before_data_m.before_data;

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
adj_cur_correlation_before = zeros(size(stimulus_data_m,1),1);
adj_cur_correlation_after = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_between = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_before_02 = zeros(size(stimulus_data_m,1),1);
adj_cur_correlation_after_02 = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_between_02 = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_before_05 = zeros(size(stimulus_data_m,1),1);
adj_cur_correlation_after_05 = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_between_05 = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_before_3 = zeros(size(stimulus_data_m,1),1);
adj_cur_correlation_after_3 = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_between_3 = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_before_5 = zeros(size(stimulus_data_m,1),1);
adj_cur_correlation_after_5 = cell(size(stimulus_data_m,1),1);
adj_cur_correlation_between_5 = cell(size(stimulus_data_m,1),1);
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
    compressed_02 = create_compressed_template(template{i},0.2);
    compressed_05 = create_compressed_template(template{i},0.5);
    compressed_3 = create_compressed_template(template{i},3);
    compressed_5 = create_compressed_template(template{i},5);

    cur_before_data = before_data_m(before_data_m(:,1) == i,:);
    cur_before_data = cur_before_data(:,2:end);
    %cur_after_data = after_data_m(after_data_m(:,1) == i,:);
    %cur_after_data = cur_after_data(:,2:end);
    cur_after_data = after_data_m{i};
    cur_after_data = cell2mat(cur_after_data);

    cur_template = template{i};
    cur_cluster = template_cluster{i};

    % check if at least 2 spikes in data, otherwise 0
    % compression 1
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
        for k = 1:size(stimulus_data_m{i,j},2)-size(cur_template,2)+1
            if sum(sum(stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1), 2) > 0) >= 2
                adj_cur_correlation_between{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(cur_template.*stimulus_data_m{i,j}(cur_cluster,k:k+size(cur_template,2)-1),'all') / (size(cur_template,1) * size(cur_template,2));
            end
        end
    end

    % compression 0.2
    for j = 1:size(cur_before_data,2)-size(compressed_02,2)+1
        if sum(sum(cur_before_data(cur_cluster,j:j+size(compressed_02,2)-1), 2) > 0) >= 2
            adj_cur_correlation_before_02(i,j) = sum(compressed_02.*cur_before_data(cur_cluster,j:j+size(compressed_02,2)-1),'all') / (size(compressed_02,1) * size(compressed_02,2));
        end
    end
    adj_cur_correlation_after{i} = zeros(1,size(cur_after_data,2)-1-size(compressed_02,2)+1);
    for j = 1:size(cur_after_data,2)-size(compressed_02,2)+1
        if sum(sum(cur_after_data(cur_cluster,j:j+size(compressed_02,2)-1), 2) > 0) >= 2
            adj_cur_correlation_after_02{i}(j) = sum(compressed_02.*cur_after_data(cur_cluster,j:j+size(compressed_02,2)-1),'all') / (size(compressed_02,1) * size(compressed_02,2));
        end
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-size(compressed_02,2)+1
            if sum(sum(stimulus_data_m{i,j}(cur_cluster,k:k+size(compressed_02,2)-1), 2) > 0) >= 2
                adj_cur_correlation_between_02{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(compressed_02.*stimulus_data_m{i,j}(cur_cluster,k:k+size(compressed_02,2)-1),'all') / (size(compressed_02,1) * size(compressed_02,2));
            end
        end
    end

    % compression 0.5
    for j = 1:size(cur_before_data,2)-size(compressed_05,2)+1
        if sum(sum(cur_before_data(cur_cluster,j:j+size(compressed_05,2)-1), 2) > 0) >= 2
            adj_cur_correlation_before_05(i,j) = sum(compressed_05.*cur_before_data(cur_cluster,j:j+size(compressed_05,2)-1),'all') / (size(compressed_05,1) * size(compressed_05,2));
        end
    end
    adj_cur_correlation_after{i} = zeros(1,size(cur_after_data,2)-1-size(compressed_05,2)+1);
    for j = 1:size(cur_after_data,2)-size(compressed_05,2)+1
        if sum(sum(cur_after_data(cur_cluster,j:j+size(compressed_05,2)-1), 2) > 0) >= 2
            adj_cur_correlation_after_05{i}(j) = sum(compressed_05.*cur_after_data(cur_cluster,j:j+size(compressed_05,2)-1),'all') / (size(compressed_05,1) * size(compressed_05,2));
        end
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-size(compressed_05,2)+1
            if sum(sum(stimulus_data_m{i,j}(cur_cluster,k:k+size(compressed_05,2)-1), 2) > 0) >= 2
                adj_cur_correlation_between_05{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(compressed_05.*stimulus_data_m{i,j}(cur_cluster,k:k+size(compressed_05,2)-1),'all') / (size(compressed_05,1) * size(compressed_05,2));
            end
        end
    end

    % compression 3
    for j = 1:size(cur_before_data,2)-size(compressed_3,2)+1
        if sum(sum(cur_before_data(cur_cluster,j:j+size(compressed_3,2)-1), 2) > 0) >= 2
            adj_cur_correlation_before_3(i,j) = sum(compressed_3.*cur_before_data(cur_cluster,j:j+size(compressed_3,2)-1),'all') / (size(compressed_3,1) * size(compressed_3,2));
        end
    end
    adj_cur_correlation_after{i} = zeros(1,size(cur_after_data,2)-1-size(compressed_3,2)+1);
    for j = 1:size(cur_after_data,2)-size(compressed_3,2)+1
        if sum(sum(cur_after_data(cur_cluster,j:j+size(compressed_3,2)-1), 2) > 0) >= 2
            adj_cur_correlation_after_3{i}(j) = sum(compressed_3.*cur_after_data(cur_cluster,j:j+size(compressed_3,2)-1),'all') / (size(compressed_3,1) * size(compressed_3,2));
        end
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-size(compressed_3,2)+1
            if sum(sum(stimulus_data_m{i,j}(cur_cluster,k:k+size(compressed_3,2)-1), 2) > 0) >= 2
                adj_cur_correlation_between_3{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(compressed_3.*stimulus_data_m{i,j}(cur_cluster,k:k+size(compressed_3,2)-1),'all') / (size(compressed_3,1) * size(compressed_3,2));
            end
        end
    end

    % compression 5
    for j = 1:size(cur_before_data,2)-size(compressed_5,2)+1
        if sum(sum(cur_before_data(cur_cluster,j:j+size(compressed_5,2)-1), 2) > 0) >= 2
            adj_cur_correlation_before_5(i,j) = sum(compressed_5.*cur_before_data(cur_cluster,j:j+size(compressed_5,2)-1),'all') / (size(compressed_5,1) * size(compressed_5,2));
        end
    end
    adj_cur_correlation_after{i} = zeros(1,size(cur_after_data,2)-1-size(compressed_5,2)+1);
    for j = 1:size(cur_after_data,2)-size(compressed_5,2)+1
        if sum(sum(cur_after_data(cur_cluster,j:j+size(compressed_5,2)-1), 2) > 0) >= 2
            adj_cur_correlation_after_5{i}(j) = sum(compressed_5.*cur_after_data(cur_cluster,j:j+size(compressed_5,2)-1),'all') / (size(compressed_5,1) * size(compressed_5,2));
        end
    end
    for j = 1:last_interval_data(i)
        for k = 1:size(stimulus_data_m{i,j},2)-size(compressed_5,2)+1
            if sum(sum(stimulus_data_m{i,j}(cur_cluster,k:k+size(compressed_5,2)-1), 2) > 0) >= 2
                adj_cur_correlation_between_5{i}((j-1)*size(stimulus_data_m{i,j},2)+k) = sum(compressed_5.*stimulus_data_m{i,j}(cur_cluster,k:k+size(compressed_5,2)-1),'all') / (size(compressed_5,1) * size(compressed_5,2));
            end
        end
    end
end

%% calculate average
% compression 1
avg_adj_cor_before = mean(adj_cur_correlation_before,2);
avg_adj_cor_between = zeros(numel(adj_cur_correlation_between),1);
avg_adj_cor_after = zeros(numel(adj_cur_correlation_after),1);
for i = 1:numel(adj_cur_correlation_between)
    avg_adj_cor_between(i) = mean(adj_cur_correlation_between{i});
    avg_adj_cor_after(i) = mean(adj_cur_correlation_after{i});
end
% compression 0.2
avg_adj_cor_before_02 = mean(adj_cur_correlation_before_02,2);
avg_adj_cor_between_02 = zeros(numel(adj_cur_correlation_between_02),1);
avg_adj_cor_after_02 = zeros(numel(adj_cur_correlation_after_02),1);
for i = 1:numel(adj_cur_correlation_between_02)
    avg_adj_cor_between_02(i) = mean(adj_cur_correlation_between_02{i});
    avg_adj_cor_after_02(i) = mean(adj_cur_correlation_after_02{i});
end
% compression 0.5
avg_adj_cor_before_05 = mean(adj_cur_correlation_before_05,2);
avg_adj_cor_between_05 = zeros(numel(adj_cur_correlation_between_05),1);
avg_adj_cor_after_05 = zeros(numel(adj_cur_correlation_after_05),1);
for i = 1:numel(adj_cur_correlation_between_05)
    avg_adj_cor_between_05(i) = mean(adj_cur_correlation_between_05{i});
    avg_adj_cor_after_05(i) = mean(adj_cur_correlation_after_05{i});
end
% compression 3
avg_adj_cor_before_3 = mean(adj_cur_correlation_before_3,2);
avg_adj_cor_between_3 = zeros(numel(adj_cur_correlation_between_3),1);
avg_adj_cor_after_3 = zeros(numel(adj_cur_correlation_after_3),1);
for i = 1:numel(adj_cur_correlation_between_3)
    avg_adj_cor_between_3(i) = mean(adj_cur_correlation_between_3{i});
    avg_adj_cor_after_3(i) = mean(adj_cur_correlation_after_3{i});
end
% compression 5
avg_adj_cor_before_5 = mean(adj_cur_correlation_before_5,2);
avg_adj_cor_between_5 = zeros(numel(adj_cur_correlation_between_5),1);
avg_adj_cor_after_5 = zeros(numel(adj_cur_correlation_after_5),1);
for i = 1:numel(adj_cur_correlation_between_5)
    avg_adj_cor_between_5(i) = mean(adj_cur_correlation_between_5{i});
    avg_adj_cor_after_5(i) = mean(adj_cur_correlation_after_5{i});
end

% create boxplots
figure
boxplot([avg_adj_cor_after_5, avg_adj_cor_after_3, avg_adj_cor_after,avg_adj_cor_after_05,avg_adj_cor_after_02], 'Labels', {'5x Compression', '3x Compression', '1x Compression', '2x Expansion', '5x Expansion'})
hold on
scatter(ones(size(avg_adj_cor_after_5,1)),avg_adj_cor_after_5, 'filled', 'blue')
scatter(ones(size(avg_adj_cor_after_3,1))*2,avg_adj_cor_after_3, 'filled', 'blue')
scatter(ones(size(avg_adj_cor_after,1))*3,avg_adj_cor_after, 'filled', 'blue')
scatter(ones(size(avg_adj_cor_after_05,1))*4,avg_adj_cor_after_05, 'filled', 'blue')
scatter(ones(size(avg_adj_cor_after_02,1))*5,avg_adj_cor_after_02, 'filled', 'blue')
line([ones(size(avg_nb_cluster_before)), ones(size(avg_nb_cluster_between))*2]',[avg_adj_cor_after_5, avg_adj_cor_after_3]','Color','green')
line([ones(size(avg_nb_cluster_between))*2, ones(size(avg_nb_cluster_after))*3]',[avg_adj_cor_after_3, avg_adj_cor_after]','Color','green')
line([ones(size(avg_nb_cluster_between))*3, ones(size(avg_nb_cluster_after))*4]',[avg_adj_cor_after, avg_adj_cor_after_05]','Color','green')
line([ones(size(avg_nb_cluster_between))*4, ones(size(avg_nb_cluster_after))*5]',[avg_adj_cor_after_05, avg_adj_cor_after_02]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_scaled_template" + i + ".png")


