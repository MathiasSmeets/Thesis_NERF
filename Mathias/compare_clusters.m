cluster_matrices_m=load("X:\Mathias\switch_data\clusters\cluster_matrices_between_m.mat");cluster_matrices_m=cluster_matrices_m.all_cluster_matrices;
cluster_matrices_np2=load("X:\Mathias\switch_data\clusters\cluster_matrices_between_np2.mat");cluster_matrices_np2=cluster_matrices_np2.all_cluster_matrices;
cluster_matrices_y=load("X:\Mathias\switch_data\clusters\cluster_matrices_between_y.mat");cluster_matrices_y=cluster_matrices_y.all_cluster_matrices;

nb_neurons_m = size(1,9);
nb_neurons_y = size(1,9);

for i = 1:9
    nb_neurons_m(i) = sum(cluster_matrices_m{i},'all')/size(cluster_matrices_m{i},2);
    if ~isempty(cluster_matrices_y{i})
    nb_neurons_y(i) = sum(cluster_matrices_y{i},'all')/size(cluster_matrices_y{i},2);
    else
        nb_neurons_y(i) = 0;
    end
end

figure
boxplot([nb_neurons_m', nb_neurons_y'])