load("X:\Mathias\results_table.mat")


%% TP
FN_ica_2 = sum(results_table.FN_ICA(results_table.NB_Neurons == 2));
TP_ica_nb_2 = sum(results_table.TP_ICA(results_table.NB_Neurons == 2));
TP_ica_nb_3 = sum(results_table.TP_ICA(results_table.NB_Neurons == 3));
TP_ica_nb_4 = sum(results_table.TP_ICA(results_table.NB_Neurons == 4));
TP_ica_nb_5 = sum(results_table.TP_ICA(results_table.NB_Neurons == 5));

TP_cpd_nb_2 = sum(results_table.TP_CPD(results_table.NB_Neurons == 2));
TP_cpd_nb_3 = sum(results_table.TP_CPD(results_table.NB_Neurons == 3));
TP_cpd_nb_4 = sum(results_table.TP_CPD(results_table.NB_Neurons == 4));
TP_cpd_nb_5 = sum(results_table.TP_CPD(results_table.NB_Neurons == 5));

TP_nmf_nb_2 = sum(results_table.TP_NMF(results_table.NB_Neurons == 2));
TP_nmf_nb_3 = sum(results_table.TP_NMF(results_table.NB_Neurons == 3));
TP_nmf_nb_4 = sum(results_table.TP_NMF(results_table.NB_Neurons == 4));
TP_nmf_nb_5 = sum(results_table.TP_NMF(results_table.NB_Neurons == 5));

TP_pca_nb_2 = sum(results_table.TP_PCA(results_table.NB_Neurons == 2));
TP_pca_nb_3 = sum(results_table.TP_PCA(results_table.NB_Neurons == 3));
TP_pca_nb_4 = sum(results_table.TP_PCA(results_table.NB_Neurons == 4));
TP_pca_nb_5 = sum(results_table.TP_PCA(results_table.NB_Neurons == 5));

TP_cpca_nb_2 = sum(results_table.TP_CPCA(results_table.NB_Neurons == 2));
TP_cpca_nb_3 = sum(results_table.TP_CPCA(results_table.NB_Neurons == 3));
TP_cpca_nb_4 = sum(results_table.TP_CPCA(results_table.NB_Neurons == 4));
TP_cpca_nb_5 = sum(results_table.TP_CPCA(results_table.NB_Neurons == 5));

TP_cica_nb_2 = sum(results_table.TP_CICA(results_table.NB_Neurons == 2));
TP_cica_nb_3 = sum(results_table.TP_CICA(results_table.NB_Neurons == 3));
TP_cica_nb_4 = sum(results_table.TP_CICA(results_table.NB_Neurons == 4));
TP_cica_nb_5 = sum(results_table.TP_CICA(results_table.NB_Neurons == 5));

%% FP
FP_ica_nb_2 = sum(results_table.FP_ICA(results_table.NB_Neurons == 2));
FP_ica_nb_3 = sum(results_table.FP_ICA(results_table.NB_Neurons == 3));
FP_ica_nb_4 = sum(results_table.FP_ICA(results_table.NB_Neurons == 4));
FP_ica_nb_5 = sum(results_table.FP_ICA(results_table.NB_Neurons == 5));

FP_cpd_nb_2 = sum(results_table.FP_CPD(results_table.NB_Neurons == 2));
FP_cpd_nb_3 = sum(results_table.FP_CPD(results_table.NB_Neurons == 3));
FP_cpd_nb_4 = sum(results_table.FP_CPD(results_table.NB_Neurons == 4));
FP_cpd_nb_5 = sum(results_table.FP_CPD(results_table.NB_Neurons == 5));

FP_nmf_nb_2 = sum(results_table.FP_NMF(results_table.NB_Neurons == 2));
FP_nmf_nb_3 = sum(results_table.FP_NMF(results_table.NB_Neurons == 3));
FP_nmf_nb_4 = sum(results_table.FP_NMF(results_table.NB_Neurons == 4));
FP_nmf_nb_5 = sum(results_table.FP_NMF(results_table.NB_Neurons == 5));

FP_pca_nb_2 = sum(results_table.FP_PCA(results_table.NB_Neurons == 2));
FP_pca_nb_3 = sum(results_table.FP_PCA(results_table.NB_Neurons == 3));
FP_pca_nb_4 = sum(results_table.FP_PCA(results_table.NB_Neurons == 4));
FP_pca_nb_5 = sum(results_table.FP_PCA(results_table.NB_Neurons == 5));

FP_cpca_nb_2 = sum(results_table.FP_CPCA(results_table.NB_Neurons == 2));
FP_cpca_nb_3 = sum(results_table.FP_CPCA(results_table.NB_Neurons == 3));
FP_cpca_nb_4 = sum(results_table.FP_CPCA(results_table.NB_Neurons == 4));
FP_cpca_nb_5 = sum(results_table.FP_CPCA(results_table.NB_Neurons == 5));

FP_cica_nb_2 = sum(results_table.FP_CICA(results_table.NB_Neurons == 2));
FP_cica_nb_3 = sum(results_table.FP_CICA(results_table.NB_Neurons == 3));
FP_cica_nb_4 = sum(results_table.FP_CICA(results_table.NB_Neurons == 4));
FP_cica_nb_5 = sum(results_table.FP_CICA(results_table.NB_Neurons == 5));

%% FN
FN_ica_nb_2 = sum(results_table.FN_ICA(results_table.NB_Neurons == 2));
FN_ica_nb_3 = sum(results_table.FN_ICA(results_table.NB_Neurons == 3));
FN_ica_nb_4 = sum(results_table.FN_ICA(results_table.NB_Neurons == 4));
FN_ica_nb_5 = sum(results_table.FN_ICA(results_table.NB_Neurons == 5));

FN_cpd_nb_2 = sum(results_table.FN_CPD(results_table.NB_Neurons == 2));
FN_cpd_nb_3 = sum(results_table.FN_CPD(results_table.NB_Neurons == 3));
FN_cpd_nb_4 = sum(results_table.FN_CPD(results_table.NB_Neurons == 4));
FN_cpd_nb_5 = sum(results_table.FN_CPD(results_table.NB_Neurons == 5));

FN_nmf_nb_2 = sum(results_table.FN_NMF(results_table.NB_Neurons == 2));
FN_nmf_nb_3 = sum(results_table.FN_NMF(results_table.NB_Neurons == 3));
FN_nmf_nb_4 = sum(results_table.FN_NMF(results_table.NB_Neurons == 4));
FN_nmf_nb_5 = sum(results_table.FN_NMF(results_table.NB_Neurons == 5));

FN_pca_nb_2 = sum(results_table.FN_PCA(results_table.NB_Neurons == 2));
FN_pca_nb_3 = sum(results_table.FN_PCA(results_table.NB_Neurons == 3));
FN_pca_nb_4 = sum(results_table.FN_PCA(results_table.NB_Neurons == 4));
FN_pca_nb_5 = sum(results_table.FN_PCA(results_table.NB_Neurons == 5));

FN_cpca_nb_2 = sum(results_table.FN_CPCA(results_table.NB_Neurons == 2));
FN_cpca_nb_3 = sum(results_table.FN_CPCA(results_table.NB_Neurons == 3));
FN_cpca_nb_4 = sum(results_table.FN_CPCA(results_table.NB_Neurons == 4));
FN_cpca_nb_5 = sum(results_table.FN_CPCA(results_table.NB_Neurons == 5));

FN_cica_nb_2 = sum(results_table.FN_CICA(results_table.NB_Neurons == 2));
FN_cica_nb_3 = sum(results_table.FN_CICA(results_table.NB_Neurons == 3));
FN_cica_nb_4 = sum(results_table.FN_CICA(results_table.NB_Neurons == 4));
FN_cica_nb_5 = sum(results_table.FN_CICA(results_table.NB_Neurons == 5));

%% TN
TN_ica_nb_2 = sum(results_table.TN_ICA(results_table.NB_Neurons == 2));
TN_ica_nb_3 = sum(results_table.TN_ICA(results_table.NB_Neurons == 3));
TN_ica_nb_4 = sum(results_table.TN_ICA(results_table.NB_Neurons == 4));
TN_ica_nb_5 = sum(results_table.TN_ICA(results_table.NB_Neurons == 5));

TN_cpd_nb_2 = sum(results_table.TN_CPD(results_table.NB_Neurons == 2));
TN_cpd_nb_3 = sum(results_table.TN_CPD(results_table.NB_Neurons == 3));
TN_cpd_nb_4 = sum(results_table.TN_CPD(results_table.NB_Neurons == 4));
TN_cpd_nb_5 = sum(results_table.TN_CPD(results_table.NB_Neurons == 5));

TN_nmf_nb_2 = sum(results_table.TN_NMF(results_table.NB_Neurons == 2));
TN_nmf_nb_3 = sum(results_table.TN_NMF(results_table.NB_Neurons == 3));
TN_nmf_nb_4 = sum(results_table.TN_NMF(results_table.NB_Neurons == 4));
TN_nmf_nb_5 = sum(results_table.TN_NMF(results_table.NB_Neurons == 5));

TN_pca_nb_2 = sum(results_table.TN_PCA(results_table.NB_Neurons == 2));
TN_pca_nb_3 = sum(results_table.TN_PCA(results_table.NB_Neurons == 3));
TN_pca_nb_4 = sum(results_table.TN_PCA(results_table.NB_Neurons == 4));
TN_pca_nb_5 = sum(results_table.TN_PCA(results_table.NB_Neurons == 5));

TN_cpca_nb_2 = sum(results_table.TN_CPCA(results_table.NB_Neurons == 2));
TN_cpca_nb_3 = sum(results_table.TN_CPCA(results_table.NB_Neurons == 3));
TN_cpca_nb_4 = sum(results_table.TN_CPCA(results_table.NB_Neurons == 4));
TN_cpca_nb_5 = sum(results_table.TN_CPCA(results_table.NB_Neurons == 5));

TN_cica_nb_2 = sum(results_table.TN_CICA(results_table.NB_Neurons == 2));
TN_cica_nb_3 = sum(results_table.TN_CICA(results_table.NB_Neurons == 3));
TN_cica_nb_4 = sum(results_table.TN_CICA(results_table.NB_Neurons == 4));
TN_cica_nb_5 = sum(results_table.TN_CICA(results_table.NB_Neurons == 5));

%% create figure

precision_pca = zeros(1,4);
precision_ica = zeros(1,4);
precision_nmf = zeros(1,4);
precision_cpd = zeros(1,4);
precision_cpca = zeros(1,4);
precision_cica = zeros(1,4);
recall_pca = zeros(1,4);
recall_ica = zeros(1,4);
recall_nmf = zeros(1,4);
recall_cpd = zeros(1,4);
recall_cpca = zeros(1,4);
recall_cica = zeros(1,4);

precision_pca(1) = TP_pca_nb_2 / (TP_pca_nb_2 + FP_pca_nb_2);
precision_pca(2) = TP_pca_nb_3 / (TP_pca_nb_3 + FP_pca_nb_3);
precision_pca(3) = TP_pca_nb_4 / (TP_pca_nb_4 + FP_pca_nb_4);
precision_pca(4) = TP_pca_nb_5 / (TP_pca_nb_5 + FP_pca_nb_5);
precision_ica(1) = TP_ica_nb_2 / (TP_ica_nb_2 + FP_ica_nb_2);
precision_ica(2) = TP_ica_nb_3 / (TP_ica_nb_3 + FP_ica_nb_3);
precision_ica(3) = TP_ica_nb_4 / (TP_ica_nb_4 + FP_ica_nb_4);
precision_ica(4) = TP_ica_nb_5 / (TP_ica_nb_5 + FP_ica_nb_5);
precision_nmf(1) = TP_nmf_nb_2 / (TP_nmf_nb_2 + FP_nmf_nb_2);
precision_nmf(2) = TP_nmf_nb_3 / (TP_nmf_nb_3 + FP_nmf_nb_3);
precision_nmf(3) = TP_nmf_nb_4 / (TP_nmf_nb_4 + FP_nmf_nb_4);
precision_nmf(4) = TP_nmf_nb_5 / (TP_nmf_nb_5 + FP_nmf_nb_5);
precision_cpd(1) = TP_cpd_nb_2 / (TP_cpd_nb_2 + FP_cpd_nb_2);
precision_cpd(2) = TP_cpd_nb_3 / (TP_cpd_nb_3 + FP_cpd_nb_3);
precision_cpd(3) = TP_cpd_nb_4 / (TP_cpd_nb_4 + FP_cpd_nb_4);
precision_cpd(4) = TP_cpd_nb_5 / (TP_cpd_nb_5 + FP_cpd_nb_5);
precision_cpca(1) = TP_cpca_nb_2 / (TP_cpca_nb_2 + FP_cpca_nb_2);
precision_cpca(2) = TP_cpca_nb_3 / (TP_cpca_nb_3 + FP_cpca_nb_3);
precision_cpca(3) = TP_cpca_nb_4 / (TP_cpca_nb_4 + FP_cpca_nb_4);
precision_cpca(4) = TP_cpca_nb_5 / (TP_cpca_nb_5 + FP_cpca_nb_5);
precision_cica(1) = TP_cica_nb_2 / (TP_cica_nb_2 + FP_cica_nb_2);
precision_cica(2) = TP_cica_nb_3 / (TP_cica_nb_3 + FP_cica_nb_3);
precision_cica(3) = TP_cica_nb_4 / (TP_cica_nb_4 + FP_cica_nb_4);
precision_cica(4) = TP_cica_nb_5 / (TP_cica_nb_5 + FP_cica_nb_5);

recall_pca(1) = TP_pca_nb_2 / (TP_pca_nb_2 + FN_pca_nb_2);
recall_pca(2) = TP_pca_nb_3 / (TP_pca_nb_3 + FN_pca_nb_3);
recall_pca(3) = TP_pca_nb_4 / (TP_pca_nb_4 + FN_pca_nb_4);
recall_pca(4) = TP_pca_nb_5 / (TP_pca_nb_5 + FN_pca_nb_5);
recall_ica(1) = TP_ica_nb_2 / (TP_ica_nb_2 + FN_ica_nb_2);
recall_ica(2) = TP_ica_nb_3 / (TP_ica_nb_3 + FN_ica_nb_3);
recall_ica(3) = TP_ica_nb_4 / (TP_ica_nb_4 + FN_ica_nb_4);
recall_ica(4) = TP_ica_nb_5 / (TP_ica_nb_5 + FN_ica_nb_5);
recall_nmf(1) = TP_nmf_nb_2 / (TP_nmf_nb_2 + FN_nmf_nb_2);
recall_nmf(2) = TP_nmf_nb_3 / (TP_nmf_nb_3 + FN_nmf_nb_3);
recall_nmf(3) = TP_nmf_nb_4 / (TP_nmf_nb_4 + FN_nmf_nb_4);
recall_nmf(4) = TP_nmf_nb_5 / (TP_nmf_nb_5 + FN_nmf_nb_5);
recall_cpd(1) = TP_cpd_nb_2 / (TP_cpd_nb_2 + FN_cpd_nb_2);
recall_cpd(2) = TP_cpd_nb_3 / (TP_cpd_nb_3 + FN_cpd_nb_3);
recall_cpd(3) = TP_cpd_nb_4 / (TP_cpd_nb_4 + FN_cpd_nb_4);
recall_cpd(4) = TP_cpd_nb_5 / (TP_cpd_nb_5 + FN_cpd_nb_5);
recall_cpca(1) = TP_cpca_nb_2 / (TP_cpca_nb_2 + FN_cpca_nb_2);
recall_cpca(2) = TP_cpca_nb_3 / (TP_cpca_nb_3 + FN_cpca_nb_3);
recall_cpca(3) = TP_cpca_nb_4 / (TP_cpca_nb_4 + FN_cpca_nb_4);
recall_cpca(4) = TP_cpca_nb_5 / (TP_cpca_nb_5 + FN_cpca_nb_5);
recall_cica(1) = TP_cica_nb_2 / (TP_cica_nb_2 + FN_cica_nb_2);
recall_cica(2) = TP_cica_nb_3 / (TP_cica_nb_3 + FN_cica_nb_3);
recall_cica(3) = TP_cica_nb_4 / (TP_cica_nb_4 + FN_cica_nb_4);
recall_cica(4) = TP_cica_nb_5 / (TP_cica_nb_5 + FN_cica_nb_5);

figure
hold on
plot(recall_pca)
plot(recall_ica)
plot(recall_nmf)
plot(recall_cpd)
plot(recall_cpca)
plot(recall_cica)
title("Recall")

figure
hold on
plot(precision_pca)
plot(precision_ica)
plot(precision_nmf)
plot(precision_cpd)
plot(precision_cpca)
plot(precision_cica)
title("Precision")
legend("PCA","ICA", "NMF", "CPD", "CPCA", "CICA","Location", "eastoutside")

figure
hold on
scatter(recall_ica, precision_ica, 100, "filled")
scatter(recall_nmf, precision_nmf, 100, "filled")
scatter(recall_cpd, precision_cpd, 100, "filled")
scatter(recall_cpca, precision_cpca, 100, "filled")
scatter(recall_cica, precision_cica, 100, "filled")
s=gcf;
xlim([0 1])
ylim([0 1])

legend("ICA", "NMF", "CPD", "CPCA", "CICA", "Location", "Eastoutside")








