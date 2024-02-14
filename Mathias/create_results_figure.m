function create_results_figure(results_table)

%% assembly plot with no missing neurons, based on amount of neurons
figure
for i = 1:4
    cur_missing_neurons = i -1;
    wrong_ica_nb_neurons_2 = sum(results_table.Wrong_Assembly_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_ica_nb_neurons_2 = sum(results_table.Correct_Assembly_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_ica_nb_neurons_3 = sum(results_table.Correct_Assembly_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_ica_nb_neurons_4 = sum(results_table.Correct_Assembly_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_ica_nb_neurons_5 = sum(results_table.Correct_Assembly_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_cpd_nb_neurons_2 = sum(results_table.Correct_Assembly_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_cpd_nb_neurons_3 = sum(results_table.Correct_Assembly_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_cpd_nb_neurons_4 = sum(results_table.Correct_Assembly_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_cpd_nb_neurons_5 = sum(results_table.Correct_Assembly_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_nmf_nb_neurons_2 = sum(results_table.Correct_Assembly_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_nmf_nb_neurons_3 = sum(results_table.Correct_Assembly_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_nmf_nb_neurons_4 = sum(results_table.Correct_Assembly_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_nmf_nb_neurons_5 = sum(results_table.Correct_Assembly_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_pca_nb_neurons_2 = sum(results_table.Correct_Assembly_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_pca_nb_neurons_3 = sum(results_table.Correct_Assembly_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_pca_nb_neurons_4 = sum(results_table.Correct_Assembly_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_pca_nb_neurons_5 = sum(results_table.Correct_Assembly_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    count_assemblies = correct_ica_nb_neurons_2 + wrong_ica_nb_neurons_2;

    input = [correct_ica_nb_neurons_2 correct_cpd_nb_neurons_2 correct_nmf_nb_neurons_2 correct_pca_nb_neurons_2; ...
        correct_ica_nb_neurons_3 correct_cpd_nb_neurons_3 correct_nmf_nb_neurons_3 correct_pca_nb_neurons_3; ...
        correct_ica_nb_neurons_4 correct_cpd_nb_neurons_4 correct_nmf_nb_neurons_4 correct_pca_nb_neurons_4; ...
        correct_ica_nb_neurons_5 correct_cpd_nb_neurons_5 correct_nmf_nb_neurons_5 correct_pca_nb_neurons_5] / count_assemblies;

    subplot(2,2,i)
    bar(2:5,input)
    hold on
    ylim([0 1])
    legend('ICA', 'CPD', 'NMF', 'PCA')
    if i == 1
        title('Performance Clustering No Missing Spikes')
    else
        cur_title = "Performance Clustering " + cur_missing_neurons + " Missing Spikes";
        title(cur_title)
    end
    xlabel('Number of Neurons in Cluster')
    ylabel('Succesful Detection Rate')

end

end