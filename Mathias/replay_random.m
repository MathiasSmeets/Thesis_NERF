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
path_to_correlations = "takeokalabwip2023/Mathias/switch_data/correlations";

stimulus_data_m = load(fullfile(volume_base2, path_to_data, "after_stimulus_data_m_horridge.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_m = stimulus_data_m(1:9,:);

raw_data_m = load(fullfile(volume_base2, "takeokalabwip2023/Mathias/switch_data/tabled_data", "horridge_data_m.mat"));
raw_data_m = raw_data_m.data;
raw_data_m(raw_data_m(:,1)>=10,:) = [];

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
template_smoothed = load(fullfile(volume_base2, path_to_correlations, "template_smoothed_3_m.mat"));template_smoothed = template_smoothed.template_smoothed;
template_cluster = load(fullfile(volume_base2, path_to_correlations, "template_cluster_m.mat"));template_cluster = template_cluster.teplate_cluster;
%mouse_to_exclude = 2;
%mouse_to_exclude = 4:9;

%% get template

interval_size = 60;
interval_size_before = 70;
intervals_together = 30;
bins_together = 15;
intervals_together_before = intervals_together*interval_size_before;
last_interval_data = zeros(1,size(stimulus_data_m,1));

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
    
    % get number of stimuli
    counter = size(stimulus_data_m,2);
    last_interval_data(i) = counter;
    while isempty(stimulus_data_m{i,counter})
        counter = counter - 1;
        last_interval_data(i) = counter;
    end

    cur_before_data = before_data_m(before_data_m(:,1) == i,:);
    cur_before_data = cur_before_data(:,2:end);
    %cur_after_data = after_data_m(after_data_m(:,1) == i,:);
    %cur_after_data = cur_after_data(:,2:end);
    cur_after_data = after_data_m{i};
    cur_after_data = cell2mat(cur_after_data);
    cur_raw_data = raw_data_m(raw_data_m(:,1) == i,:);
    cur_raw_data = cur_raw_data(:,2:end);

    cur_template = template_smoothed{i};
    %cur_template = template{i};
    cur_cluster = template_cluster{i};

    % do randomization 10 times
    for iteration = 1:10
        % create random template
        random_timepoints = randi([1,size(cur_raw_data,2)-size(cur_template)+1], [1,last_interval_data(i)]);
        % get template by averaging activity from each timepoint
        new_template = zeros(size(cur_template));
        for cur_timepoint = 1:numel(random_timepoints)
            new_template = new_template + cur_raw_data(cur_cluster,cur_timepoint:cur_timepoint+size(new_template,2)-1);
        end
        new_template = new_template / numel(random_timepoints);
        error("stop")
    end
    error("stopp")

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
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_bba_random.png")

figure
boxplot([avg_cor_before, avg_cor_after], 'Labels', {'Baseline', 'Rest'})
hold on
scatter(ones(size(avg_cor_before,1)),avg_cor_before, 'filled', 'blue')
scatter(ones(size(avg_cor_after,1))*2,avg_cor_after, 'filled', 'blue')
line([ones(size(avg_cor_before)), ones(size(avg_cor_after))*2]',[avg_cor_before, avg_cor_after]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_ba_random.png")

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
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_adjusted_bba_random.png")

figure
boxplot([avg_adj_cur_correlation_before, avg_adj_cur_correlation_after], 'Labels', {'Baseline', 'Rest'})
hold on
scatter(ones(size(avg_adj_cur_correlation_before,1)),avg_adj_cur_correlation_before, 'filled', 'blue')
scatter(ones(size(avg_adj_cur_correlation_after,1))*2,avg_adj_cur_correlation_after, 'filled', 'blue')
line([ones(size(avg_adj_cur_correlation_before)), ones(size(avg_adj_cur_correlation_after))*2]',[avg_adj_cur_correlation_before, avg_adj_cur_correlation_after]','Color','green')
saveas(gcf,"/scratch/mathiass-takeokalab/01/boxplot_adjusted_ba_random.png")

save("/scratch/mathiass-takeokalab/01/correlation_before_smoothed_width3_random.mat","avg_adj_cur_correlation_before")
save("/scratch/mathiass-takeokalab/01/correlation_between_smoothed_width3_random.mat","avg_adj_cur_correlation_between")
save("/scratch/mathiass-takeokalab/01/correlation_after_smoothed_width3_random.mat","avg_adj_cur_correlation_after")
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
% [p,h,stats] = signrank(avg_cor_before,avg_cor_after)

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

