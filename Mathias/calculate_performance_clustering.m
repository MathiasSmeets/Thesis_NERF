%clear; clc; close all;

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

%%
nb_iterations = 30;
nb_intervals = 10;
bins_together = 10;
interval_length = 7;
nb_neurons_of_interest = 10;

largest_delay = 6;
iterations_nb_neurons = 5;
iterations_missing_neurons = 5;

nb_assembly_array = [];
nb_neurons_array = [];
nb_missing_neurons_array = [];
nb_iterations_array = [];
correct_assemblies_ica = [];
wrong_assemblies_ica = [];
TP_ica = [];
FP_ica = [];
FN_ica = [];
correct_assemblies_cpd = [];
wrong_assemblies_cpd = [];
TP_cpd = [];
FP_cpd = [];
FN_cpd = [];
correct_assemblies_nmf = [];
wrong_assemblies_nmf = [];
TP_nmf = [];
FP_nmf = [];
FN_nmf = [];
correct_assemblies_pca = [];
wrong_assemblies_pca = [];
TP_pca = [];
FP_pca = [];
FN_pca = [];

for nb_assemblies = 1:1
    for nb_neurons = 2:iterations_nb_neurons
        disp(nb_neurons)
        for missing_neurons = 0:iterations_missing_neurons
            for i = 1:nb_iterations
                [synthetic_data, neurons_in_assembly, activations, synthetic_data_non_zscore] = generate_synthetic_data(stimulus_data_m, sos_results_m, nb_assemblies, nb_neurons, missing_neurons, nb_intervals, bins_together, nb_neurons_of_interest, 0);

                % ica
                [ica_predicted_nbr_assemblies, ica_predicted_nbr_neurons, ica_assemblies, ica_activity] = ica_assembly_detection(synthetic_data', 0);

                % cpd
                synthetic_tensor = create_tensor(synthetic_data, interval_length);
                [cpd_predicted_nbr_assemblies, cpd_predicted_nbr_neurons, cpd_assemblies, cpd_activity] = cpd_assembly_detection(synthetic_tensor, ica_assemblies);

                % nmf
                [nmf_predicted_nbr_assemblies, nmf_predicted_nbr_neurons, nmf_assemblies, nmf_activity] = nmf_assembly_detection(synthetic_data_non_zscore, ica_assemblies, 50);

                % convolutive pca
                [synthetic_data_pca, neurons_in_assembly_pca, activations_pca, ~] = generate_synthetic_data(stimulus_data_m, sos_results_m, nb_assemblies, nb_neurons, missing_neurons, nb_intervals, 1, nb_neurons_of_interest, 1);
                synthetic_data_delayed_pca = create_convolutive_data(synthetic_data_pca, largest_delay);
                [predicted_nbr_assemblies_pca, predicted_nbr_neurons_pca, original_assemblies_pca,activity_pca] = pca_assembly_detection(synthetic_data_delayed_pca', 0);
                pca_assemblies = assemblies_delays(original_assemblies_pca, largest_delay);
                
                % check performance
                [ica_nb_assembly_correct, ica_nb_assembly_wrong, ica_nb_neurons_TP, ica_nb_neurons_FP, ica_nb_neurons_FN] = get_performance_clustering({neurons_in_assembly}, ica_assemblies);
                [cpd_nb_assembly_correct, cpd_nb_assembly_wrong, cpd_nb_neurons_TP, cpd_nb_neurons_FP, cpd_nb_neurons_FN] = get_performance_clustering({neurons_in_assembly}, cpd_assemblies);
                [nmf_nb_assembly_correct, nmf_nb_assembly_wrong, nmf_nb_neurons_TP, nmf_nb_neurons_FP, nmf_nb_neurons_FN] = get_performance_clustering({neurons_in_assembly}, nmf_assemblies);
                [pca_nb_assembly_correct, pca_nb_assembly_wrong, pca_nb_neurons_TP, pca_nb_neurons_FP, pca_nb_neurons_FN] = get_performance_clustering({neurons_in_assembly_pca}, pca_assemblies);

                % assign everything to proper arrays
                nb_assembly_array = [nb_assembly_array, nb_assemblies];
                nb_neurons_array = [nb_neurons_array, nb_neurons];
                nb_missing_neurons_array = [nb_missing_neurons_array, missing_neurons];
                nb_iterations_array = [nb_iterations_array, i];
                correct_assemblies_ica = [correct_assemblies_ica, ica_nb_assembly_correct];
                wrong_assemblies_ica = [wrong_assemblies_ica, ica_nb_assembly_wrong];
                TP_ica = [TP_ica, ica_nb_neurons_TP];
                FP_ica = [FP_ica, ica_nb_neurons_FP];
                FN_ica = [FN_ica, ica_nb_neurons_FN];
                correct_assemblies_cpd = [correct_assemblies_cpd, cpd_nb_assembly_correct];
                wrong_assemblies_cpd = [wrong_assemblies_cpd, cpd_nb_assembly_wrong];
                TP_cpd = [TP_cpd, cpd_nb_neurons_TP];
                FP_cpd = [FP_cpd, cpd_nb_neurons_FP];
                FN_cpd = [FN_cpd, cpd_nb_neurons_FN];
                correct_assemblies_nmf = [correct_assemblies_nmf, nmf_nb_assembly_correct];
                wrong_assemblies_nmf = [wrong_assemblies_nmf, nmf_nb_assembly_wrong];
                TP_nmf = [TP_nmf, nmf_nb_neurons_TP];
                FP_nmf = [FP_nmf, nmf_nb_neurons_FP];
                FN_nmf = [FN_nmf, nmf_nb_neurons_FN];
                correct_assemblies_pca = [correct_assemblies_pca, pca_nb_assembly_correct];
                wrong_assemblies_pca = [wrong_assemblies_pca, pca_nb_assembly_wrong];
                TP_pca = [TP_pca, pca_nb_neurons_TP];
                FP_pca = [FP_pca, pca_nb_neurons_FP];
                FN_pca = [FN_pca, pca_nb_neurons_FN];
            end
        end
    end
end
results_table = table(nb_assembly_array', nb_neurons_array', nb_missing_neurons_array', nb_iterations_array', ...
    correct_assemblies_ica', wrong_assemblies_ica', TP_ica', FP_ica', FN_ica', ...
    correct_assemblies_cpd', wrong_assemblies_cpd', TP_cpd', FP_cpd', FN_cpd', ...
    correct_assemblies_nmf', wrong_assemblies_nmf', TP_nmf', FP_nmf', FN_nmf', ...
    correct_assemblies_pca', wrong_assemblies_pca', TP_pca', FP_pca', FN_pca', ...
    'VariableNames',{'NB_Assemblies','NB_Neurons','NB_Missing_Neurons','NB_Iterations', ...
    'Correct_Assembly_ICA','Wrong_Assembly_ICA', 'TP_ICA', 'FP_ICA', 'FN_ICA', ...
    'Correct_Assembly_CPD','Wrong_Assembly_CPD', 'TP_CPD', 'FP_CPD', 'FN_CPD', ...
    'Correct_Assembly_NMF','Wrong_Assembly_NMF', 'TP_NMF', 'FP_NMF', 'FN_NMF', ...
    'Correct_Assembly_PCA','Wrong_Assembly_PCA', 'TP_PCA', 'FP_PCA', 'FN_PCA'});

create_results_figure(results_table)
