clear;clc;close all;
load("X:\Mathias\results_table_assemblies.mat")

for cur_nb_bins_together = 15
    for cur_nb_intervals_together = 30
        %% TP
        TP_ica_nb_2 = sum(results_table.TP_assemblies_ICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_ica_nb_3 = sum(results_table.TP_assemblies_ICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_ica_nb_4 = sum(results_table.TP_assemblies_ICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_ica_nb_5 = sum(results_table.TP_assemblies_ICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TP_cpd_nb_2 = sum(results_table.TP_assemblies_CPD(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_cpd_nb_3 = sum(results_table.TP_assemblies_CPD(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_cpd_nb_4 = sum(results_table.TP_assemblies_CPD(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_cpd_nb_5 = sum(results_table.TP_assemblies_CPD(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TP_nmf_nb_2 = sum(results_table.TP_assemblies_NMF(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_nmf_nb_3 = sum(results_table.TP_assemblies_NMF(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_nmf_nb_4 = sum(results_table.TP_assemblies_NMF(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_nmf_nb_5 = sum(results_table.TP_assemblies_NMF(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TP_pca_nb_2 = sum(results_table.TP_assemblies_PCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_pca_nb_3 = sum(results_table.TP_assemblies_PCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_pca_nb_4 = sum(results_table.TP_assemblies_PCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        TP_pca_nb_5 = sum(results_table.TP_assemblies_PCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        TP_cpca_nb_2 = sum(results_table.TP_assemblies_CPCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cpca_nb_3 = sum(results_table.TP_assemblies_CPCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cpca_nb_4 = sum(results_table.TP_assemblies_CPCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cpca_nb_5 = sum(results_table.TP_assemblies_CPCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        TP_cica_nb_2 = sum(results_table.TP_assemblies_CICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cica_nb_3 = sum(results_table.TP_assemblies_CICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cica_nb_4 = sum(results_table.TP_assemblies_CICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        TP_cica_nb_5 = sum(results_table.TP_assemblies_CICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        %% FP
        FP_ica_nb_2 = sum(results_table.FP_assemblies_ICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_ica_nb_3 = sum(results_table.FP_assemblies_ICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_ica_nb_4 = sum(results_table.FP_assemblies_ICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_ica_nb_5 = sum(results_table.FP_assemblies_ICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FP_cpd_nb_2 = sum(results_table.FP_assemblies_CPD(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_cpd_nb_3 = sum(results_table.FP_assemblies_CPD(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_cpd_nb_4 = sum(results_table.FP_assemblies_CPD(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_cpd_nb_5 = sum(results_table.FP_assemblies_CPD(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FP_nmf_nb_2 = sum(results_table.FP_assemblies_NMF(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_nmf_nb_3 = sum(results_table.FP_assemblies_NMF(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_nmf_nb_4 = sum(results_table.FP_assemblies_NMF(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_nmf_nb_5 = sum(results_table.FP_assemblies_NMF(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FP_pca_nb_2 = sum(results_table.FP_assemblies_PCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_pca_nb_3 = sum(results_table.FP_assemblies_PCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_pca_nb_4 = sum(results_table.FP_assemblies_PCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));
        FP_pca_nb_5 = sum(results_table.FP_assemblies_PCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together & results_table.NB_Bins == cur_nb_bins_together));

        FP_cpca_nb_2 = sum(results_table.FP_assemblies_CPCA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cpca_nb_3 = sum(results_table.FP_assemblies_CPCA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cpca_nb_4 = sum(results_table.FP_assemblies_CPCA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cpca_nb_5 = sum(results_table.FP_assemblies_CPCA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        FP_cica_nb_2 = sum(results_table.FP_assemblies_CICA(results_table.NB_Neurons == 2 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cica_nb_3 = sum(results_table.FP_assemblies_CICA(results_table.NB_Neurons == 3 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cica_nb_4 = sum(results_table.FP_assemblies_CICA(results_table.NB_Neurons == 4 & results_table.NB_Intervals == cur_nb_intervals_together));
        FP_cica_nb_5 = sum(results_table.FP_assemblies_CICA(results_table.NB_Neurons == 5 & results_table.NB_Intervals == cur_nb_intervals_together));

        pca_detected = TP_pca_nb_2 + TP_pca_nb_3 + TP_pca_nb_4 + TP_pca_nb_5 + FP_pca_nb_2 + FP_pca_nb_3 + FP_pca_nb_4 + FP_pca_nb_5;
        ica_detected = TP_ica_nb_2 + TP_ica_nb_3 + TP_ica_nb_4 + TP_ica_nb_5 + FP_ica_nb_2 + FP_ica_nb_3 + FP_ica_nb_4 + FP_ica_nb_5;
        nmf_detected = TP_nmf_nb_2 + TP_nmf_nb_3 + TP_nmf_nb_4 + TP_nmf_nb_5 + FP_nmf_nb_2 + FP_nmf_nb_3 + FP_nmf_nb_4 + FP_nmf_nb_5;
        cpd_detected = TP_cpd_nb_2 + TP_cpd_nb_3 + TP_cpd_nb_4 + TP_cpd_nb_5 + FP_cpd_nb_2 + FP_cpd_nb_3 + FP_cpd_nb_4 + FP_cpd_nb_5;
        cpca_detected = TP_cpca_nb_2 + TP_cpca_nb_3 + TP_cpca_nb_4 + TP_cpca_nb_5 + FP_cpca_nb_2 + FP_cpca_nb_3 + FP_cpca_nb_4 + FP_cpca_nb_5;
        cica_detected = TP_cica_nb_2 + TP_cica_nb_3 + TP_cica_nb_4 + TP_cica_nb_5 + FP_cica_nb_2 + FP_cica_nb_3 + FP_cica_nb_4 + FP_cica_nb_5;

    end
end

bar([cpd_detected/400, cpca_detected/400, cica_detected/400])
hold on
errorbar(1:3,[cpd_detected/400, cpca_detected/400, cica_detected/400], [std(results_table.TP_assemblies_PCA), std(results_table.TP_assemblies_CPCA), std(results_table.TP_assemblies_CICA)])