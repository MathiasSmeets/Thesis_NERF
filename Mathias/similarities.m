clear;clc;close all;
path1 = "X:\Mathias\switch_data\clusters";
path2 = "X:\Mathias\switch_data\correlations";
cluster_matrix_before = load(fullfile(path1,"cluster_matrices_before_m.mat")); cluster_matrix_before = cluster_matrix_before.all_cluster_matrices;
cluster_matrix_between = load(fullfile(path1,"cluster_matrices_between_m.mat")); cluster_matrix_between = cluster_matrix_between.all_cluster_matrices;
cluster_matrix_after = load(fullfile(path1,"cluster_matrices_after_m.mat")); cluster_matrix_after = cluster_matrix_after.all_cluster_matrices;
cluster_matrix_horridge = load(fullfile(path1,"cluster_matrices_switch_m.mat")); cluster_matrix_horridge = cluster_matrix_horridge.all_cluster_matrices;

clusters = load(fullfile(path2, "template_cluster_m.mat"));clusters = clusters.template_cluster;


similarities_before = cell(1,9);
similarities_between = cell(1,9);
similarities_after = cell(1,9);
similarities_horridge = cell(1,9);
for i = 1:numel(clusters)
    vector = zeros(size(cluster_matrix_before{i},1),1);
    vector(clusters{i}) = 1;

    % Compute cosine similarity between vector and each column of the matrix
    vector_magnitude = sqrt(sum(vector.^2)); % Magnitude of the vector
    matrix_magnitude_before = sqrt(sum(cluster_matrix_before{i}.^2)); % Magnitude of each column of the matrix
    dot_products_before = sum(vector .* cluster_matrix_before{i}); % Dot product between vector and each column of the matrix

    matrix_magnitude_between = sqrt(sum(cluster_matrix_between{i}.^2)); % Magnitude of each column of the matrix
    dot_products_between = sum(vector .* cluster_matrix_between{i}); % Dot product between vector and each column of the matrix

    matrix_magnitude_after = sqrt(sum(cluster_matrix_after{i}.^2)); % Magnitude of each column of the matrix
    dot_products_after = sum(vector .* cluster_matrix_after{i}); % Dot product between vector and each column of the matrix

    matrix_magnitude_horridge = sqrt(sum(cluster_matrix_horridge{i}.^2)); % Magnitude of each column of the matrix
    dot_products_horridge = sum(vector .* cluster_matrix_horridge{i}); % Dot product between vector and each column of the matrix

    similarities_before{i} = dot_products_before ./ (vector_magnitude * matrix_magnitude_before);
    similarities_between{i} = dot_products_between ./ (vector_magnitude * matrix_magnitude_between);
    similarities_after{i} = dot_products_after ./ (vector_magnitude * matrix_magnitude_after);
    similarities_horridge{i} = dot_products_horridge ./ (vector_magnitude * matrix_magnitude_horridge);

    zero_indices_before = (matrix_magnitude_before==0);
    zero_indices_between = (matrix_magnitude_between==0);
    zero_indices_after = (matrix_magnitude_after==0);
    zero_indices_horridge = (matrix_magnitude_horridge==0);

    similarities_before{i}(zero_indices_before) = 0;
    similarities_between{i}(zero_indices_between) = 0;
    similarities_after{i}(zero_indices_after) = 0;
    similarities_horridge{i}(zero_indices_horridge) = 0;
end


%% get averages

averages = zeros(9,4);
for i = 1:9
    averages(i,1) = mean(similarities_before{i});
    averages(i,2) = mean(similarities_between{i});
    averages(i,3) = mean(similarities_after{i});
    averages(i,4) = mean(similarities_horridge{i});
end

figure
boxplot([averages(:,1), averages(:,2), averages(:,3), averages(:,4)])
p = signrank(averages(:,1), averages(:,3));
