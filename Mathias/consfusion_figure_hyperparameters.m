clear;clc;close all;
load("X:\Mathias\results_table_hyperparameters.mat")

precision_pca = zeros(25,4);
precision_ica = zeros(25,4);
precision_nmf = zeros(25,4);
precision_cpd = zeros(25,4);
precision_cpca = zeros(25,4);
precision_cica = zeros(25,4);
recall_pca = zeros(25,4);
recall_ica = zeros(25,4);
recall_nmf = zeros(25,4);
recall_cpd = zeros(25,4);
recall_cpca = zeros(25,4);
recall_cica = zeros(25,4);
counter = 1;

for cur_nb_bins_together = [1,5,10,15,20]
    for cur_nb_intervals_together = [5,15,30,45,60]
        %% TP
        TP_ica_nb_2 = sum(results_table.TP_ICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_ica_nb_3 = sum(results_table.TP_ICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_ica_nb_4 = sum(results_table.TP_ICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_ica_nb_5 = sum(results_table.TP_ICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TP_cpd_nb_2 = sum(results_table.TP_CPD(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_cpd_nb_3 = sum(results_table.TP_CPD(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_cpd_nb_4 = sum(results_table.TP_CPD(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_cpd_nb_5 = sum(results_table.TP_CPD(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TP_nmf_nb_2 = sum(results_table.TP_NMF(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_nmf_nb_3 = sum(results_table.TP_NMF(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_nmf_nb_4 = sum(results_table.TP_NMF(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_nmf_nb_5 = sum(results_table.TP_NMF(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TP_pca_nb_2 = sum(results_table.TP_PCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_pca_nb_3 = sum(results_table.TP_PCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_pca_nb_4 = sum(results_table.TP_PCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_pca_nb_5 = sum(results_table.TP_PCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TP_cpca_nb_2 = sum(results_table.TP_CPCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cpca_nb_3 = sum(results_table.TP_CPCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cpca_nb_4 = sum(results_table.TP_CPCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cpca_nb_5 = sum(results_table.TP_CPCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        TP_cica_nb_2 = sum(results_table.TP_CICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cica_nb_3 = sum(results_table.TP_CICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cica_nb_4 = sum(results_table.TP_CICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cica_nb_5 = sum(results_table.TP_CICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        %% FP
        FP_ica_nb_2 = sum(results_table.FP_ICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_ica_nb_3 = sum(results_table.FP_ICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_ica_nb_4 = sum(results_table.FP_ICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_ica_nb_5 = sum(results_table.FP_ICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FP_cpd_nb_2 = sum(results_table.FP_CPD(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_cpd_nb_3 = sum(results_table.FP_CPD(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_cpd_nb_4 = sum(results_table.FP_CPD(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_cpd_nb_5 = sum(results_table.FP_CPD(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FP_nmf_nb_2 = sum(results_table.FP_NMF(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_nmf_nb_3 = sum(results_table.FP_NMF(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_nmf_nb_4 = sum(results_table.FP_NMF(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_nmf_nb_5 = sum(results_table.FP_NMF(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FP_pca_nb_2 = sum(results_table.FP_PCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_pca_nb_3 = sum(results_table.FP_PCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_pca_nb_4 = sum(results_table.FP_PCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_pca_nb_5 = sum(results_table.FP_PCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FP_cpca_nb_2 = sum(results_table.FP_CPCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cpca_nb_3 = sum(results_table.FP_CPCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cpca_nb_4 = sum(results_table.FP_CPCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cpca_nb_5 = sum(results_table.FP_CPCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        FP_cica_nb_2 = sum(results_table.FP_CICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cica_nb_3 = sum(results_table.FP_CICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cica_nb_4 = sum(results_table.FP_CICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cica_nb_5 = sum(results_table.FP_CICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        %% FN
        FN_ica_nb_2 = sum(results_table.FN_ICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_ica_nb_3 = sum(results_table.FN_ICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_ica_nb_4 = sum(results_table.FN_ICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_ica_nb_5 = sum(results_table.FN_ICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FN_cpd_nb_2 = sum(results_table.FN_CPD(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_cpd_nb_3 = sum(results_table.FN_CPD(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_cpd_nb_4 = sum(results_table.FN_CPD(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_cpd_nb_5 = sum(results_table.FN_CPD(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FN_nmf_nb_2 = sum(results_table.FN_NMF(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_nmf_nb_3 = sum(results_table.FN_NMF(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_nmf_nb_4 = sum(results_table.FN_NMF(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_nmf_nb_5 = sum(results_table.FN_NMF(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FN_pca_nb_2 = sum(results_table.FN_PCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_pca_nb_3 = sum(results_table.FN_PCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_pca_nb_4 = sum(results_table.FN_PCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FN_pca_nb_5 = sum(results_table.FN_PCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FN_cpca_nb_2 = sum(results_table.FN_CPCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        FN_cpca_nb_3 = sum(results_table.FN_CPCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        FN_cpca_nb_4 = sum(results_table.FN_CPCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        FN_cpca_nb_5 = sum(results_table.FN_CPCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        FN_cica_nb_2 = sum(results_table.FN_CICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        FN_cica_nb_3 = sum(results_table.FN_CICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        FN_cica_nb_4 = sum(results_table.FN_CICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        FN_cica_nb_5 = sum(results_table.FN_CICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        %% TN
        TN_ica_nb_2 = sum(results_table.TN_ICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_ica_nb_3 = sum(results_table.TN_ICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_ica_nb_4 = sum(results_table.TN_ICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_ica_nb_5 = sum(results_table.TN_ICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TN_cpd_nb_2 = sum(results_table.TN_CPD(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_cpd_nb_3 = sum(results_table.TN_CPD(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_cpd_nb_4 = sum(results_table.TN_CPD(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_cpd_nb_5 = sum(results_table.TN_CPD(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TN_nmf_nb_2 = sum(results_table.TN_NMF(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_nmf_nb_3 = sum(results_table.TN_NMF(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_nmf_nb_4 = sum(results_table.TN_NMF(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_nmf_nb_5 = sum(results_table.TN_NMF(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TN_pca_nb_2 = sum(results_table.TN_PCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_pca_nb_3 = sum(results_table.TN_PCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_pca_nb_4 = sum(results_table.TN_PCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TN_pca_nb_5 = sum(results_table.TN_PCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TN_cpca_nb_2 = sum(results_table.TN_CPCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        TN_cpca_nb_3 = sum(results_table.TN_CPCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        TN_cpca_nb_4 = sum(results_table.TN_CPCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        TN_cpca_nb_5 = sum(results_table.TN_CPCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        TN_cica_nb_2 = sum(results_table.TN_CICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        TN_cica_nb_3 = sum(results_table.TN_CICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        TN_cica_nb_4 = sum(results_table.TN_CICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        TN_cica_nb_5 = sum(results_table.TN_CICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        %% create figure

        precision_pca(counter,1) = TP_pca_nb_2 / (TP_pca_nb_2 + FP_pca_nb_2);
        precision_pca(counter,2) = TP_pca_nb_3 / (TP_pca_nb_3 + FP_pca_nb_3);
        precision_pca(counter,3) = TP_pca_nb_4 / (TP_pca_nb_4 + FP_pca_nb_4);
        precision_pca(counter,4) = TP_pca_nb_5 / (TP_pca_nb_5 + FP_pca_nb_5);
        precision_ica(counter,1) = TP_ica_nb_2 / (TP_ica_nb_2 + FP_ica_nb_2);
        precision_ica(counter,2) = TP_ica_nb_3 / (TP_ica_nb_3 + FP_ica_nb_3);
        precision_ica(counter,3) = TP_ica_nb_4 / (TP_ica_nb_4 + FP_ica_nb_4);
        precision_ica(counter,4) = TP_ica_nb_5 / (TP_ica_nb_5 + FP_ica_nb_5);
        precision_nmf(counter,1) = TP_nmf_nb_2 / (TP_nmf_nb_2 + FP_nmf_nb_2);
        precision_nmf(counter,2) = TP_nmf_nb_3 / (TP_nmf_nb_3 + FP_nmf_nb_3);
        precision_nmf(counter,3) = TP_nmf_nb_4 / (TP_nmf_nb_4 + FP_nmf_nb_4);
        precision_nmf(counter,4) = TP_nmf_nb_5 / (TP_nmf_nb_5 + FP_nmf_nb_5);
        precision_cpd(counter,1) = TP_cpd_nb_2 / (TP_cpd_nb_2 + FP_cpd_nb_2);
        precision_cpd(counter,2) = TP_cpd_nb_3 / (TP_cpd_nb_3 + FP_cpd_nb_3);
        precision_cpd(counter,3) = TP_cpd_nb_4 / (TP_cpd_nb_4 + FP_cpd_nb_4);
        precision_cpd(counter,4) = TP_cpd_nb_5 / (TP_cpd_nb_5 + FP_cpd_nb_5);
        precision_cpca(counter,1) = TP_cpca_nb_2 / (TP_cpca_nb_2 + FP_cpca_nb_2);
        precision_cpca(counter,2) = TP_cpca_nb_3 / (TP_cpca_nb_3 + FP_cpca_nb_3);
        precision_cpca(counter,3) = TP_cpca_nb_4 / (TP_cpca_nb_4 + FP_cpca_nb_4);
        precision_cpca(counter,4) = TP_cpca_nb_5 / (TP_cpca_nb_5 + FP_cpca_nb_5);
        precision_cica(counter,1) = TP_cica_nb_2 / (TP_cica_nb_2 + FP_cica_nb_2);
        precision_cica(counter,2) = TP_cica_nb_3 / (TP_cica_nb_3 + FP_cica_nb_3);
        precision_cica(counter,3) = TP_cica_nb_4 / (TP_cica_nb_4 + FP_cica_nb_4);
        precision_cica(counter,4) = TP_cica_nb_5 / (TP_cica_nb_5 + FP_cica_nb_5);

        recall_pca(counter,1) = TP_pca_nb_2 / (TP_pca_nb_2 + FN_pca_nb_2);
        recall_pca(counter,2) = TP_pca_nb_3 / (TP_pca_nb_3 + FN_pca_nb_3);
        recall_pca(counter,3) = TP_pca_nb_4 / (TP_pca_nb_4 + FN_pca_nb_4);
        recall_pca(counter,4) = TP_pca_nb_5 / (TP_pca_nb_5 + FN_pca_nb_5);
        recall_ica(counter,1) = TP_ica_nb_2 / (TP_ica_nb_2 + FN_ica_nb_2);
        recall_ica(counter,2) = TP_ica_nb_3 / (TP_ica_nb_3 + FN_ica_nb_3);
        recall_ica(counter,3) = TP_ica_nb_4 / (TP_ica_nb_4 + FN_ica_nb_4);
        recall_ica(counter,4) = TP_ica_nb_5 / (TP_ica_nb_5 + FN_ica_nb_5);
        recall_nmf(counter,1) = TP_nmf_nb_2 / (TP_nmf_nb_2 + FN_nmf_nb_2);
        recall_nmf(counter,2) = TP_nmf_nb_3 / (TP_nmf_nb_3 + FN_nmf_nb_3);
        recall_nmf(counter,3) = TP_nmf_nb_4 / (TP_nmf_nb_4 + FN_nmf_nb_4);
        recall_nmf(counter,4) = TP_nmf_nb_5 / (TP_nmf_nb_5 + FN_nmf_nb_5);
        recall_cpd(counter,1) = TP_cpd_nb_2 / (TP_cpd_nb_2 + FN_cpd_nb_2);
        recall_cpd(counter,2) = TP_cpd_nb_3 / (TP_cpd_nb_3 + FN_cpd_nb_3);
        recall_cpd(counter,3) = TP_cpd_nb_4 / (TP_cpd_nb_4 + FN_cpd_nb_4);
        recall_cpd(counter,4) = TP_cpd_nb_5 / (TP_cpd_nb_5 + FN_cpd_nb_5);
        recall_cpca(counter,1) = TP_cpca_nb_2 / (TP_cpca_nb_2 + FN_cpca_nb_2);
        recall_cpca(counter,2) = TP_cpca_nb_3 / (TP_cpca_nb_3 + FN_cpca_nb_3);
        recall_cpca(counter,3) = TP_cpca_nb_4 / (TP_cpca_nb_4 + FN_cpca_nb_4);
        recall_cpca(counter,4) = TP_cpca_nb_5 / (TP_cpca_nb_5 + FN_cpca_nb_5);
        recall_cica(counter,1) = TP_cica_nb_2 / (TP_cica_nb_2 + FN_cica_nb_2);
        recall_cica(counter,2) = TP_cica_nb_3 / (TP_cica_nb_3 + FN_cica_nb_3);
        recall_cica(counter,3) = TP_cica_nb_4 / (TP_cica_nb_4 + FN_cica_nb_4);
        recall_cica(counter,4) = TP_cica_nb_5 / (TP_cica_nb_5 + FN_cica_nb_5);
        
        counter = counter + 1;
    end
end

figure
hold on
counter = 0;
values_transparancy = [0.2 0.4 0.6 0.8 1];
blue_values = {"#2d3fa6ff", };
for i = 1:5:size(recall_ica,1)-1
    counter = counter +1;
    scatter(mean(recall_ica(i,:), 2), mean(precision_ica(i,:), 2), 200, "filled", 'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter))
    scatter(mean(recall_nmf(i,:), 2), mean(precision_nmf(i,:), 2), 200, "filled", 'MarkerFaceColor', [0.8500 0.3250 0.0980], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter))
    scatter(mean(recall_cpd(i,:), 2), mean(precision_cpd(i,:), 2), 200, "filled", 'MarkerFaceColor', [0.9290 0.6940 0.1250], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter))

    scatter(mean(recall_ica(i+1,:), 2), mean(precision_ica(i+4,:), 2), 200, 'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "Diamond")
    scatter(mean(recall_nmf(i+1,:), 2), mean(precision_nmf(i+4,:), 2), 200, 'MarkerFaceColor', [0.8500 0.3250 0.0980], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "Diamond")
    scatter(mean(recall_cpd(i+1,:), 2), mean(precision_cpd(i+4,:), 2), 200, 'MarkerFaceColor', [0.9290 0.6940 0.1250], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "Diamond")
    
    scatter(mean(recall_ica(i+2,:), 2), mean(precision_ica(i+1,:), 2), 200, 'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "Square")
    scatter(mean(recall_nmf(i+2,:), 2), mean(precision_nmf(i+1,:), 2), 200, 'MarkerFaceColor', [0.8500 0.3250 0.0980], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "Square")
    scatter(mean(recall_cpd(i+2,:), 2), mean(precision_cpd(i+1,:), 2), 200, 'MarkerFaceColor', [0.9290 0.6940 0.1250], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "Square")

    scatter(mean(recall_ica(i+3,:), 2), mean(precision_ica(i+2,:), 2), 200, 'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "v")
    scatter(mean(recall_nmf(i+3,:), 2), mean(precision_nmf(i+2,:), 2), 200, 'MarkerFaceColor', [0.8500 0.3250 0.0980], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "v")
    scatter(mean(recall_cpd(i+3,:), 2), mean(precision_cpd(i+2,:), 2), 200, 'MarkerFaceColor', [0.9290 0.6940 0.1250], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "v")

    scatter(mean(recall_ica(i+4,:), 2), mean(precision_ica(i+3,:), 2), 200, 'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "^")
    scatter(mean(recall_nmf(i+4,:), 2), mean(precision_nmf(i+3,:), 2), 200, 'MarkerFaceColor', [0.8500 0.3250 0.0980], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "^")
    scatter(mean(recall_cpd(i+4,:), 2), mean(precision_cpd(i+3,:), 2), 200, 'MarkerFaceColor', [0.9290 0.6940 0.1250], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', values_transparancy(counter), "Marker", "^")


end
scatter(mean(recall_cpca(1,:), 2), mean(precision_cpca(1,:), 2), 200, 'MarkerFaceColor', [0.4940 0.1840 0.5560], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2)
scatter(mean(recall_cpca(2,:), 2), mean(precision_cpca(2,:), 2), 200, 'MarkerFaceColor', [0.4940 0.1840 0.5560], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2, "Marker", "Diamond")
scatter(mean(recall_cpca(3,:), 2), mean(precision_cpca(3,:), 2), 200, 'MarkerFaceColor', [0.4940 0.1840 0.5560], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2, "Marker", "Square")
scatter(mean(recall_cpca(4,:), 2), mean(precision_cpca(4,:), 2), 200, 'MarkerFaceColor', [0.4940 0.1840 0.5560], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2, "Marker", "v")
scatter(mean(recall_cpca(5,:), 2), mean(precision_cpca(5,:), 2), 200, 'MarkerFaceColor', [0.4940 0.1840 0.5560], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2, "Marker", "^")

scatter(mean(recall_cica(1,:), 2), mean(precision_cica(1,:), 2), 200, 'MarkerFaceColor', [0.4660 0.6740 0.1880], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2)
scatter(mean(recall_cica(2,:), 2), mean(precision_cica(2,:), 2), 200, 'MarkerFaceColor', [0.4660 0.6740 0.1880], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2, "Marker", "Diamond")
scatter(mean(recall_cica(3,:), 2), mean(precision_cica(3,:), 2), 200, 'MarkerFaceColor', [0.4660 0.6740 0.1880], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2, "Marker", "Square")
scatter(mean(recall_cica(4,:), 2), mean(precision_cica(4,:), 2), 200, 'MarkerFaceColor', [0.4660 0.6740 0.1880], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2, "Marker", "v")
scatter(mean(recall_cica(5,:), 2), mean(precision_cica(5,:), 2), 200, 'MarkerFaceColor', [0.4660 0.6740 0.1880], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2, "Marker", "^")


s=gcf;
xlim([0 1])
ylim([0 1])
h = gca;
h.LineWidth = 2;

%% create dummy plot for legend
figure
hold on
scatter(1,2,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0])
scatter(1,3,200,"filled",'MarkerFaceColor', [0.8500 0.3250 0.0980], "MarkerEdgeColor", [0 0 0])
scatter(1,4,200,"filled",'MarkerFaceColor', [0.9290 0.6940 0.1250], "MarkerEdgeColor", [0 0 0])
scatter(1,5,200,"filled",'MarkerFaceColor', [0.4940 0.1840 0.5560], "MarkerEdgeColor", [0 0 0])
scatter(1,6,200,"filled",'MarkerFaceColor', [0.4660 0.6740 0.1880], "MarkerEdgeColor", [0 0 0])

scatter(1,7,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.2)
scatter(1,8,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.4)
scatter(1,9,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.6)
scatter(1,10,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 0.8)
scatter(1,11,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], 'MarkerFaceAlpha', 1)

scatter(1,12,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0])
scatter(1,13,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], "Marker", "Diamond")
scatter(1,14,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], "Marker", "Square")
scatter(1,15,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], "Marker", "v")
scatter(1,16,200,"filled",'MarkerFaceColor', [0 0.4470 0.7410], "MarkerEdgeColor", [0 0 0], "Marker", "^")


lgd = legend("PCA + ICA", "NMF", "CPD", "CPCA", "CICA", "1ms Bins", "5ms Bins", "10ms Bins", "15ms Bins", "20ms Bins", "5 Intervals", "15 Intervals", "30 Intervals", "45 Intervals", "60 Intervals", "Location", "Eastoutside", "NumColumns", 3);
fontsize(lgd, 14, 'points')
h = gca;
h.LineWidth = 2;




