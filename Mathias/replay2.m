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

ica_assemblies = load(fullfile(volume_base2,path_to_clusters, "assemblies_horridge_m.mat")); ica_assemblies = ica_assemblies.total_assemblies;
ica_data = load(fullfile(volume_base2,path_to_clusters, "data_horridge_m.mat")); ica_data = ica_data.total_data;
ica_neurons_of_interest = load(fullfile(volume_base2,path_to_clusters, "neurons_of_interest_horridge_m.mat")); ica_neurons_of_interest = ica_neurons_of_interest.total_neurons_of_interest;
ica_activity = load(fullfile(volume_base2,path_to_clusters, "activity_horridge_m.mat")); ica_activity = ica_activity.total_activity;
ica_vector = load(fullfile(volume_base2,path_to_clusters, "ica_vector_horridge_m.mat")); ica_vector = ica_vector.total_vector;

ica_neurons_of_interest_before = load(fullfile(volume_base2, path_to_clusters, "neurons_of_interest_after_m.mat"));ica_neurons_of_interest_before = ica_neurons_of_interest_before.total_neurons_of_interest;
ica_assemblies_before = load(fullfile(volume_base2, path_to_clusters, "assemblies_after_m.mat")); ica_assemblies_before = ica_assemblies_before.total_assemblies;

template = load(fullfile(volume_base2, path_to_clusters, "template_m.mat"));template = template.template;
template_cluster = load(fullfile(volume_base2, path_to_clusters, "template_m.mat"));template_cluster = template_cluster.template_cluster;
%% get template

interval_size = 70;
intervals_together = 30;
bins_together = 15;
intervals_together_before = intervals_together*interval_size;
last_interval_data = zeros(1,size(stimulus_data_m,1));

for i = 1:size(stimulus_data_m,1)
    % get last interval
    % get last interval
    counter = size(stimulus_data_m,2);
    last_interval_data(i) = counter;
    while isempty(stimulus_data_m{i,counter})
        counter = counter - 1;
        last_interval_data(i) = counter;
    end

end