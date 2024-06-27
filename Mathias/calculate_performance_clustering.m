clear; clc; close all;

%% get data

if ispc
    volume_base = '\\nerffs13';
    volume_base2 = '\\nerffs17';
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/';
    volume_base2 = '/mnt/takeokalab/';
end

path_to_code = "takeokalabwip2023/Mathias/data";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

stimulus_data_m = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_m.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_y = load(fullfile(volume_base2, path_to_code,"data_after_stimulus_y.mat"));
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

neurons_of_interest_m = load(fullfile(volume_base2, path_to_code, "output_m.mat"));
neurons_of_interest_m = neurons_of_interest_m.output_m;
inhibited_neurons_m = load(fullfile(volume_base2, path_to_code, "inhibited_m.mat"));
inhibited_neurons_m = inhibited_neurons_m.inhibited_m;

sos_results_m = load(fullfile(volume_base2, path_to_code, "sos_results_m.mat"));
sos_results_m = sos_results_m.sos_results_m;

folder = fileparts(which("calculate_performance_clustering.m"));
addpath(genpath(folder))

%% initialize variables
nb_iterations = 100;
nb_intervals = 10;
bins_together = 10;
interval_length = 7;
nb_neurons_of_interest = 10;
max_bins_together = 26;
max_intervals_together = 50;
missing_neurons = 0;

largest_delay = 6;
iterations_nb_neurons = 5;
iterations_missing_neurons = 0;

nb_assembly_array = [];
nb_neurons_array = [];
nb_missing_neurons_array = [];
nb_iterations_array = [];
background_strength = [];
nb_bins_array = [];
nb_intervals_array = [];
TP_assemblies_ica = [];
FP_assemblies_ica = [];
FN_assemblies_ica = [];
TP_ica = [];
FP_ica = [];
FN_ica = [];
TN_ica = [];
TP_assemblies_cpd = [];
FP_assemblies_cpd = [];
FN_assemblies_cpd = [];
TP_cpd = [];
FP_cpd = [];
FN_cpd = [];
TN_cpd = [];
TP_assemblies_nmf = [];
FP_assemblies_nmf = [];
FN_assemblies_nmf = [];
TP_nmf = [];
FP_nmf = [];
FN_nmf = [];
TN_nmf = [];
TP_assemblies_pca = [];
FP_assemblies_pca = [];
FN_assemblies_pca = [];
TP_pca = [];
FP_pca = [];
FN_pca = [];
TN_pca = [];
TP_assemblies_cpca = [];
FP_assemblies_cpca = [];
FN_assemblies_cpca = [];
TP_cpca = [];
FP_cpca = [];
FN_cpca = [];
TN_cpca = [];
TP_assemblies_cica = [];
FP_assemblies_cica = [];
FN_assemblies_cica = [];
TP_cica = [];
FP_cica = [];
FN_cica = [];
TN_cica = [];
%% loop over sets of hyperparameters and calculate performance
for nb_assemblies = 1:1
    for nb_neurons = 4:iterations_nb_neurons
        %for missing_neurons = 0:iterations_missing_neurons
        for cur_nb_bins_together = 15%[1,5,10,15,20]
            for cur_nb_intervals_together = 30%[5,15,30,45,60]
                disp("neurons: " + nb_neurons + ", bins: " + cur_nb_bins_together + ", intervals: " + cur_nb_intervals_together)
                for cur_background_strength = 1:1
                    for i = 1:nb_iterations
                        [synthetic_data, neurons_in_assembly, activations, synthetic_data_non_zscore] = generate_synthetic_data(stimulus_data_m, sos_results_m, nb_assemblies, nb_neurons, missing_neurons, cur_nb_intervals_together, cur_nb_bins_together, nb_neurons_of_interest, 1, cur_background_strength);
                        
                        % pca
                        [pca_predicted_nbr_assemblies, pca_predicted_nbr_neurons, pca_assemblies, pca_activity, ~] = ica_assembly_detection(synthetic_data', 0);

                        % ica
                        [ica_predicted_nbr_assemblies, ica_predicted_nbr_neurons, ica_assemblies, ica_activity, ~] = ica_assembly_detection(synthetic_data', 0);
                        if ica_predicted_nbr_assemblies == 1
                            disp("test1")
                        elseif ica_predicted_nbr_assemblies == 2
                            disp("test2")
                        end
 
                        % cpd
                        interval_length = ceil(70/cur_nb_bins_together);
                        synthetic_tensor = create_tensor(synthetic_data, interval_length);
                        [cpd_predicted_nbr_assemblies, cpd_predicted_nbr_neurons, cpd_assemblies, cpd_activity] = cpd_assembly_detection(synthetic_tensor, ica_assemblies);

                        % nmf
                        [nmf_predicted_nbr_assemblies, nmf_predicted_nbr_neurons, nmf_assemblies, nmf_activity] = nmf_assembly_detection(synthetic_data_non_zscore, ica_assemblies, 50);

                        % convolutive pca
                        [synthetic_data_cpca, neurons_in_assembly_cpca, activations_cpca, ~] = generate_synthetic_data(stimulus_data_m, sos_results_m, nb_assemblies, nb_neurons, missing_neurons, cur_nb_intervals_together, 1, nb_neurons_of_interest, 1, cur_background_strength);
                        synthetic_data_delayed_cpca = create_convolutive_data(synthetic_data_cpca, largest_delay);
                        [predicted_nbr_assemblies_cpca, predicted_nbr_neurons_cpca, original_assemblies_cpca,activity_cpca] = pca_assembly_detection(synthetic_data_delayed_cpca', 0);
                        cpca_assemblies = assemblies_delays(original_assemblies_cpca, largest_delay);

                        % convolutive ica
                        [predicted_nbr_assemblies_cica, predicted_nbr_neurons_cica, original_assemblies_cica,activity_cica] = ica_assembly_detection(synthetic_data_delayed_cpca', 0);
                        cica_assemblies = assemblies_delays(original_assemblies_cica, largest_delay);

                        % check performance
                        [ica_nb_assembly_TP, ica_nb_assembly_FP, ica_nb_assembly_FN, ica_nb_neurons_TP, ica_nb_neurons_FP, ica_nb_neurons_FN, ica_nb_neurons_TN] = get_performance_clustering(neurons_in_assembly, ica_assemblies);
                        [cpd_nb_assembly_TP, cpd_nb_assembly_FP, cpd_nb_assembly_FN, cpd_nb_neurons_TP, cpd_nb_neurons_FP, cpd_nb_neurons_FN, cpd_nb_neurons_TN] = get_performance_clustering(neurons_in_assembly, cpd_assemblies);
                        [nmf_nb_assembly_TP, nmf_nb_assembly_FP, nmf_nb_assembly_FN, nmf_nb_neurons_TP, nmf_nb_neurons_FP, nmf_nb_neurons_FN, nmf_nb_neurons_TN] = get_performance_clustering(neurons_in_assembly, nmf_assemblies);
                        [pca_nb_assembly_TP, pca_nb_assembly_FP, pca_nb_assembly_FN, pca_nb_neurons_TP, pca_nb_neurons_FP, pca_nb_neurons_FN, pca_nb_neurons_TN] = get_performance_clustering(neurons_in_assembly, pca_assemblies);
                        [cpca_nb_assembly_TP, cpca_nb_assembly_FP, cpca_nb_assembly_FN, cpca_nb_neurons_TP, cpca_nb_neurons_FP, cpca_nb_neurons_FN, cpca_nb_neurons_TN] = get_performance_clustering(neurons_in_assembly_cpca, cpca_assemblies);
                        [cica_nb_assembly_TP, cica_nb_assembly_FP, cica_nb_assembly_FN, cica_nb_neurons_TP, cica_nb_neurons_FP, cica_nb_neurons_FN, cica_nb_neurons_TN] = get_performance_clustering(neurons_in_assembly_cpca, cica_assemblies);
    
                        % if cpca_nb_assembly_TP ~= cica_nb_assembly_TP
                        %     disp("test")
                        % end

                        % assign everything to proper arrays
                        nb_assembly_array = [nb_assembly_array, nb_assemblies];
                        nb_neurons_array = [nb_neurons_array, nb_neurons];
                        %nb_missing_neurons_array = [nb_missing_neurons_array, missing_neurons];
                        nb_bins_array = [nb_bins_array, cur_nb_bins_together];
                        nb_intervals_array = [nb_intervals_array, cur_nb_intervals_together];
                        nb_iterations_array = [nb_iterations_array, i];
                        background_strength = [background_strength, cur_background_strength];
                        TP_assemblies_ica = [TP_assemblies_ica, ica_nb_assembly_TP];
                        FP_assemblies_ica = [FP_assemblies_ica, ica_nb_assembly_FP];
                        FN_assemblies_ica = [FN_assemblies_ica, ica_nb_assembly_FN];
                        TP_ica = [TP_ica, ica_nb_neurons_TP];
                        FP_ica = [FP_ica, ica_nb_neurons_FP];
                        FN_ica = [FN_ica, ica_nb_neurons_FN];
                        TN_ica = [TN_ica, ica_nb_neurons_TN];
                        TP_assemblies_cpd = [TP_assemblies_cpd, cpd_nb_assembly_TP];
                        FP_assemblies_cpd = [FP_assemblies_cpd, cpd_nb_assembly_FP];
                        FN_assemblies_cpd = [FN_assemblies_cpd, cpd_nb_assembly_FN];
                        TP_cpd = [TP_cpd, cpd_nb_neurons_TP];
                        FP_cpd = [FP_cpd, cpd_nb_neurons_FP];
                        FN_cpd = [FN_cpd, cpd_nb_neurons_FN];
                        TN_cpd = [TN_cpd, cpd_nb_neurons_TN];
                        TP_assemblies_nmf = [TP_assemblies_nmf, nmf_nb_assembly_TP];
                        FP_assemblies_nmf = [FP_assemblies_nmf, nmf_nb_assembly_FP];
                        FN_assemblies_nmf = [FN_assemblies_nmf, nmf_nb_assembly_FN];
                        TP_nmf = [TP_nmf, nmf_nb_neurons_TP];
                        FP_nmf = [FP_nmf, nmf_nb_neurons_FP];
                        FN_nmf = [FN_nmf, nmf_nb_neurons_FN];
                        TN_nmf = [TN_nmf, nmf_nb_neurons_TN];
                        TP_assemblies_pca = [TP_assemblies_pca, pca_nb_assembly_TP];
                        FP_assemblies_pca = [FP_assemblies_pca, pca_nb_assembly_FP];
                        FN_assemblies_pca = [FN_assemblies_pca, pca_nb_assembly_FN];
                        TP_pca = [TP_pca, pca_nb_neurons_TP];
                        FP_pca = [FP_pca, pca_nb_neurons_FP];
                        FN_pca = [FN_pca, pca_nb_neurons_FN];
                        TN_pca = [TN_pca, pca_nb_neurons_TN];
                        TP_assemblies_cpca = [TP_assemblies_cpca, cpca_nb_assembly_TP];
                        FP_assemblies_cpca = [FP_assemblies_cpca, cpca_nb_assembly_FP];
                        FN_assemblies_cpca = [FN_assemblies_cpca, cpca_nb_assembly_FN];
                        TP_cpca = [TP_cpca, cpca_nb_neurons_TP];
                        FP_cpca = [FP_cpca, cpca_nb_neurons_FP];
                        FN_cpca = [FN_cpca, cpca_nb_neurons_FN];
                        TN_cpca = [TN_cpca, cpca_nb_neurons_TN];
                        TP_assemblies_cica = [TP_assemblies_cica, cica_nb_assembly_TP];
                        FP_assemblies_cica = [FP_assemblies_cica, cica_nb_assembly_FP];
                        FN_assemblies_cica = [FN_assemblies_cica, cica_nb_assembly_FN];
                        TP_cica = [TP_cica, cica_nb_neurons_TP];
                        FP_cica = [FP_cica, cica_nb_neurons_FP];
                        FN_cica = [FN_cica, cica_nb_neurons_FN];
                        TN_cica = [TN_cica, cica_nb_neurons_TN];
                    end
                end
            end
        end
    end
end
%TP_assemblies_pca', FP_assemblies_pca', FN_assemblies_pca', TP_pca', FP_pca', FN_pca'
results_table = table(nb_assembly_array', nb_neurons_array', nb_iterations_array', background_strength', nb_intervals_array', nb_bins_array',...
    TP_assemblies_pca', FP_assemblies_pca', FN_assemblies_pca', TP_pca', FP_pca', FN_pca', TN_pca', ...
    TP_assemblies_ica', FP_assemblies_ica', FN_assemblies_ica', TP_ica', FP_ica', FN_ica', TN_ica', ...
    TP_assemblies_cpd', FP_assemblies_cpd', FN_assemblies_cpd', TP_cpd', FP_cpd', FN_cpd', TN_cpd', ...
    TP_assemblies_nmf', FP_assemblies_nmf', FN_assemblies_nmf', TP_nmf', FP_nmf', FN_nmf', TN_nmf', ...
    TP_assemblies_cpca', FP_assemblies_cpca', FN_assemblies_cpca', TP_cpca', FP_cpca', FN_cpca', TN_cpca', ...
    TP_assemblies_cica', FP_assemblies_cica', FN_assemblies_cica', TP_cica', FP_cica', FN_cica', TN_cica', ...
    'VariableNames',{'NB_Assemblies','NB_Neurons','NB_Iterations', 'Background_Strength', 'NB_Intervals', 'NB_Bins',...
    'TP_assemblies_PCA','FP_assemblies_PCA', 'FN_assemblies_PCA', 'TP_PCA', 'FP_PCA', 'FN_PCA', 'TN_PCA',...
    'TP_assemblies_ICA','FP_assemblies_ICA', 'FN_assemblies_ICA', 'TP_ICA', 'FP_ICA', 'FN_ICA', 'TN_ICA',...
    'TP_assemblies_CPD','FP_assemblies_CPD', 'FN_assemblies_CPD', 'TP_CPD', 'FP_CPD', 'FN_CPD', 'TN_CPD',...
    'TP_assemblies_NMF','FP_assemblies_NMF', 'FN_assemblies_NMF', 'TP_NMF', 'FP_NMF', 'FN_NMF', 'TN_NMF',...
    'TP_assemblies_CPCA', 'FP_assemblies_CPCA', 'FN_assemblies_CPCA', 'TP_CPCA', 'FP_CPCA', 'FN_CPCA', 'TN_CPCA',...
    'TP_assemblies_CICA', 'FP_assemblies_CICA', 'FN_assemblies_CICA', 'TP_CICA', 'FP_CICA', 'FN_CICA', 'TN_CICA'});
    %'TP_assemblies_PCA','FP_assemblies_PCA', 'FN_assemblies_PCA', 'TP_PCA', 'FP_PCA', 'FN_PCA'
%save("X:\Mathias\cluster_output\results_table", "results_table", "-v7.3");
% save("/scratch/mathiass-takeokalab/01/results_table", "results_table", "-v7.3");
% create_results_figure_intervals_bins(results_table)
% create_results_figure(results_table)
