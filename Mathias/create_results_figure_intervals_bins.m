function create_results_figure_intervals_bins(results_table)

nb_iterations = max(results_table.NB_Iterations);

%% TP assembly plot with no missing neurons, based on amount of neurons
max_intervals_together = max(results_table.NB_Intervals);
max_bins_together = 26;
nb_rows =length(10:10:max_intervals_together);
nb_columns = length(1:5:max_bins_together);
input = zeros(nb_rows, nb_columns);
figure
for i = 2:5
    row_counter = 1;
    for cur_nb_intervals_together = 10:10:max_intervals_together
        column_counter = 1;
        for cur_nb_bins_together = 1:5:max_bins_together
            if cur_nb_bins_together ~= 1
                cur_cur_nb_bins_together = cur_nb_bins_together - 1;
            else
                cur_cur_nb_bins_together = 1;
            end
            correct_ica_nb_neurons_2 = sum(results_table.TP_assemblies_ICA(results_table.NB_Bins == cur_cur_nb_bins_together & results_table.NB_Neurons == i & results_table.NB_Intervals == cur_nb_intervals_together));

            correct_cpd_nb_neurons_2 = sum(results_table.TP_assemblies_CPD(results_table.NB_Bins == cur_cur_nb_bins_together & results_table.NB_Neurons == i & results_table.NB_Intervals == cur_nb_intervals_together));

            correct_nmf_nb_neurons_2 = sum(results_table.TP_assemblies_NMF(results_table.NB_Bins == cur_cur_nb_bins_together & results_table.NB_Neurons == i & results_table.NB_Intervals == cur_nb_intervals_together));

            % correct_pca_nb_neurons_2 = sum(results_table.TP_assemblies_PCA(results_table.NB_Bins == cur_nb_bins_together & results_table.NB_Neurons == i & results_table.NB_Intervals == cur_nb_intervals_together));
            % correct_pca_nb_neurons_3 = sum(results_table.TP_assemblies_PCA(results_table.NB_Bins == cur_nb_bins_together & results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
            % correct_pca_nb_neurons_4 = sum(results_table.TP_assemblies_PCA(results_table.NB_Bins == cur_nb_bins_together & results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
            % correct_pca_nb_neurons_5 = sum(results_table.TP_assemblies_PCA(results_table.NB_Bins == cur_nb_bins_together & results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

            nb_iterations = max(results_table.NB_Iterations);
            input(row_counter,column_counter) = correct_ica_nb_neurons_2/ nb_iterations;
            column_counter = column_counter + 1;
                
            % input = [correct_ica_nb_neurons_2 correct_cpd_nb_neurons_2 correct_nmf_nb_neurons_2; ...
            %     correct_ica_nb_neurons_3 correct_cpd_nb_neurons_3 correct_nmf_nb_neurons_3; ...
            %     correct_ica_nb_neurons_4 correct_cpd_nb_neurons_4 correct_nmf_nb_neurons_4; ...
            %     correct_ica_nb_neurons_5 correct_cpd_nb_neurons_5 correct_nmf_nb_neurons_5] / nb_iterations;
        end
        row_counter = row_counter + 1;
    end
    subplot(2,2,i-1)
    row_values = 10:10:max_intervals_together;
    column_values = [1 5 10 15 20 25];
    cur_h = heatmap(column_values, row_values, input);
    
    cur_title = "TP rate with " + i + " neurons";
    cur_h.Title = cur_title;
    
    cur_h.XLabel = 'Bin Size (ms)';
    cur_h.YLabel = 'Intervals';
end

%% FP assembly plot with no missing neurons, based on amount of neurons

figure
input = zeros(nb_rows, nb_columns);
for i = 2:5
    row_counter = 1;
    for cur_nb_intervals_together = 10:10:max_intervals_together
        column_counter = 1;
        for cur_nb_bins_together = 1:5:max_bins_together
            if cur_nb_bins_together ~= 1
                cur_cur_nb_bins_together = cur_nb_bins_together - 1;
            else
                cur_cur_nb_bins_together = 1;
            end
            correct_ica_nb_neurons_2 = sum(results_table.FP_assemblies_ICA(results_table.NB_Bins == cur_cur_nb_bins_together & results_table.NB_Neurons == i & results_table.NB_Intervals == cur_nb_intervals_together));

            correct_cpd_nb_neurons_2 = sum(results_table.FP_assemblies_CPD(results_table.NB_Bins == cur_cur_nb_bins_together & results_table.NB_Neurons == i & results_table.NB_Intervals == cur_nb_intervals_together));

            correct_nmf_nb_neurons_2 = sum(results_table.FP_assemblies_NMF(results_table.NB_Bins == cur_cur_nb_bins_together & results_table.NB_Neurons == i & results_table.NB_Intervals == cur_nb_intervals_together));

            % correct_pca_nb_neurons_2 = sum(results_table.TP_assemblies_PCA(results_table.NB_Bins == cur_nb_bins_together & results_table.NB_Neurons == i & results_table.NB_Intervals == cur_nb_intervals_together));
            % correct_pca_nb_neurons_3 = sum(results_table.TP_assemblies_PCA(results_table.NB_Bins == cur_nb_bins_together & results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
            % correct_pca_nb_neurons_4 = sum(results_table.TP_assemblies_PCA(results_table.NB_Bins == cur_nb_bins_together & results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
            % correct_pca_nb_neurons_5 = sum(results_table.TP_assemblies_PCA(results_table.NB_Bins == cur_nb_bins_together & results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

            nb_iterations = max(results_table.NB_Iterations);
            input(row_counter,column_counter) = correct_ica_nb_neurons_2/ nb_iterations;
            column_counter = column_counter + 1;
                
            % input = [correct_ica_nb_neurons_2 correct_cpd_nb_neurons_2 correct_nmf_nb_neurons_2; ...
            %     correct_ica_nb_neurons_3 correct_cpd_nb_neurons_3 correct_nmf_nb_neurons_3; ...
            %     correct_ica_nb_neurons_4 correct_cpd_nb_neurons_4 correct_nmf_nb_neurons_4; ...
            %     correct_ica_nb_neurons_5 correct_cpd_nb_neurons_5 correct_nmf_nb_neurons_5] / nb_iterations;
        end
        row_counter = row_counter + 1;
    end
    subplot(2,2,i-1)
    row_values = 10:10:max_intervals_together;
    column_values = [1 5 10 15 20 25];
    cur_h = heatmap(column_values, row_values, input);
    
    cur_title = "FP rate with " + i + " neurons";
    cur_h.Title = cur_title;
    
    cur_h.XLabel = 'Bin Size (ms)';
    cur_h.YLabel = 'Intervals';
end

end