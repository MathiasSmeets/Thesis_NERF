function create_results_figure(results_table)

nb_iterations = max(results_table.NB_Iterations);

%% TP assembly plot with no missing neurons, based on amount of neurons
figure
for i = 1:4
    cur_missing_neurons = i -1;
    FN_assemblies_ica_2 = sum(results_table.FN_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_ica_nb_neurons_2 = sum(results_table.TP_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_ica_nb_neurons_3 = sum(results_table.TP_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_ica_nb_neurons_4 = sum(results_table.TP_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_ica_nb_neurons_5 = sum(results_table.TP_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_cpd_nb_neurons_2 = sum(results_table.TP_assemblies_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_cpd_nb_neurons_3 = sum(results_table.TP_assemblies_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_cpd_nb_neurons_4 = sum(results_table.TP_assemblies_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_cpd_nb_neurons_5 = sum(results_table.TP_assemblies_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_nmf_nb_neurons_2 = sum(results_table.TP_assemblies_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_nmf_nb_neurons_3 = sum(results_table.TP_assemblies_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_nmf_nb_neurons_4 = sum(results_table.TP_assemblies_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_nmf_nb_neurons_5 = sum(results_table.TP_assemblies_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_pca_nb_neurons_2 = sum(results_table.TP_assemblies_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_pca_nb_neurons_3 = sum(results_table.TP_assemblies_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_pca_nb_neurons_4 = sum(results_table.TP_assemblies_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_pca_nb_neurons_5 = sum(results_table.TP_assemblies_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    count_assemblies = correct_ica_nb_neurons_2 + FN_assemblies_ica_2;

    input = [correct_ica_nb_neurons_2 correct_cpd_nb_neurons_2 correct_nmf_nb_neurons_2 correct_pca_nb_neurons_2; ...
        correct_ica_nb_neurons_3 correct_cpd_nb_neurons_3 correct_nmf_nb_neurons_3 correct_pca_nb_neurons_3; ...
        correct_ica_nb_neurons_4 correct_cpd_nb_neurons_4 correct_nmf_nb_neurons_4 correct_pca_nb_neurons_4; ...
        correct_ica_nb_neurons_5 correct_cpd_nb_neurons_5 correct_nmf_nb_neurons_5 correct_pca_nb_neurons_5] / nb_iterations;
    
    % input = [correct_ica_nb_neurons_2 correct_cpd_nb_neurons_2 correct_nmf_nb_neurons_2; ...
    %     correct_ica_nb_neurons_3 correct_cpd_nb_neurons_3 correct_nmf_nb_neurons_3; ...
    %     correct_ica_nb_neurons_4 correct_cpd_nb_neurons_4 correct_nmf_nb_neurons_4; ...
    %     correct_ica_nb_neurons_5 correct_cpd_nb_neurons_5 correct_nmf_nb_neurons_5] / nb_iterations;


    subplot(2,2,i)
    bar(2:5,input)
    hold on
    ylim([0 1])
    legend('ICA', 'CPD', 'NMF', 'PCA')
    %legend('ICA', 'CPD', 'NMF')
    if i == 1
        title('True Positive Rate Clustering No Missing Spikes')
    else
        cur_title = "True Positive Rate Clustering " + cur_missing_neurons + " Missing Spikes";
        title(cur_title)
    end
    xlabel('Number of Neurons in Cluster')
    ylabel('Succesful Detection Rate')

end

%% FP assembly plot with no missing neurons, based on amount of neurons

figure
for i = 1:4
    cur_missing_neurons = i -1;
    FN_assemblies_ica_2 = sum(results_table.FN_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_ica_nb_neurons_2 = sum(results_table.FP_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_ica_nb_neurons_3 = sum(results_table.FP_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_ica_nb_neurons_4 = sum(results_table.FP_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_ica_nb_neurons_5 = sum(results_table.FP_assemblies_ICA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_cpd_nb_neurons_2 = sum(results_table.FP_assemblies_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_cpd_nb_neurons_3 = sum(results_table.FP_assemblies_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_cpd_nb_neurons_4 = sum(results_table.FP_assemblies_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_cpd_nb_neurons_5 = sum(results_table.FP_assemblies_CPD(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_nmf_nb_neurons_2 = sum(results_table.FP_assemblies_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_nmf_nb_neurons_3 = sum(results_table.FP_assemblies_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_nmf_nb_neurons_4 = sum(results_table.FP_assemblies_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_nmf_nb_neurons_5 = sum(results_table.FP_assemblies_NMF(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    correct_pca_nb_neurons_2 = sum(results_table.FP_assemblies_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 2));
    correct_pca_nb_neurons_3 = sum(results_table.FP_assemblies_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 3));
    correct_pca_nb_neurons_4 = sum(results_table.FP_assemblies_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 4));
    correct_pca_nb_neurons_5 = sum(results_table.FP_assemblies_PCA(results_table.NB_Missing_Neurons == cur_missing_neurons & results_table.NB_Neurons == 5));

    count_assemblies = correct_ica_nb_neurons_2 + FN_assemblies_ica_2;

    input = [correct_ica_nb_neurons_2 correct_cpd_nb_neurons_2 correct_nmf_nb_neurons_2 correct_pca_nb_neurons_2; ...
        correct_ica_nb_neurons_3 correct_cpd_nb_neurons_3 correct_nmf_nb_neurons_3 correct_pca_nb_neurons_3; ...
        correct_ica_nb_neurons_4 correct_cpd_nb_neurons_4 correct_nmf_nb_neurons_4 correct_pca_nb_neurons_4; ...
        correct_ica_nb_neurons_5 correct_cpd_nb_neurons_5 correct_nmf_nb_neurons_5 correct_pca_nb_neurons_5] / nb_iterations;

    % input = [correct_ica_nb_neurons_2 correct_cpd_nb_neurons_2 correct_nmf_nb_neurons_2; ...
    %     correct_ica_nb_neurons_3 correct_cpd_nb_neurons_3 correct_nmf_nb_neurons_3; ...
    %     correct_ica_nb_neurons_4 correct_cpd_nb_neurons_4 correct_nmf_nb_neurons_4; ...
    %     correct_ica_nb_neurons_5 correct_cpd_nb_neurons_5 correct_nmf_nb_neurons_5] / nb_iterations;

    subplot(2,2,i)
    bar(2:5,input)
    hold on
    legend('ICA', 'CPD', 'NMF', 'PCA')
    % legend('ICA', 'CPD', 'NMF')
    if i == 1
        title('False Positive Clustering No Missing Spikes')
    else
        cur_title = "False Positive Clustering " + cur_missing_neurons + " Missing Spikes";
        title(cur_title)
    end
    xlabel('Number of Neurons in Cluster')
    ylabel('Assemblies Wrongly Detected')

end

end