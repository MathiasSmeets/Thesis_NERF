path1 = "X:\Mathias\switch_data\clusters";
path2 = "X:\Mathias\switch_data\data_after_stimulus";
stimulus_data_m = load(fullfile(path2, "after_stimulus_data_m_switch"));stimulus_data_m=stimulus_data_m.after_stimulus_switch_m;
total_assemblies = load(fullfile(path1,"assemblies_switch_m.mat"));total_assemblies = total_assemblies.total_assemblies;
total_neurons_of_interest = load(fullfile(path1, "neurons_of_interest_switch_m.mat"));total_neurons_of_interest = total_neurons_of_interest.total_neurons_of_interest;
total_data = load(fullfile(path1, "data_switch_m.mat"));total_data = total_data.total_data;

interval_step = 30;
population_similarities_start = zeros(2+1,size(total_assemblies,1));
population_similarities_end = zeros(2+1,size(total_assemblies,1));
all_cluster_matrices = cell(1,size(stimulus_data_m,1));
% loop over all mice
for i = 1:size(stimulus_data_m,1)
    % get last interval
    j = size(total_data,2);
    while isempty(total_data{i,j}) && j > 1
        j = j-1;
    end
    last_interval_index = j;

    % create clustermatrix
    cluster_matrix = zeros(size(stimulus_data_m{i,1},1), last_interval_index);
    for k = 1:last_interval_index
        for l = 1:length(total_assemblies{i,k})
            cluster_matrix(total_neurons_of_interest{i,k}(total_assemblies{i,k}{l}),k) = 1;
        end
    end
    figure
    imagesc(cluster_matrix)
    % xlabel("Time Bins (" + interval_step + " together)")
    % ylabel("Neurons")
    % title("Neurons in a cluster")
    % all_cluster_matrices{i} = cluster_matrix;
end