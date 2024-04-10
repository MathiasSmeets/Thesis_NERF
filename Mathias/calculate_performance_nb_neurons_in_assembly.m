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

%%
nb_iterations = 50;
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
TP_assemblies_cpd = [];
FP_assemblies_cpd = [];
FN_assemblies_cpd = [];
TP_cpd = [];
FP_cpd = [];
FN_cpd = [];
TP_assemblies_nmf = [];
FP_assemblies_nmf = [];
FN_assemblies_nmf = [];
TP_nmf = [];
FP_nmf = [];
FN_nmf = [];
TP_assemblies_pca = [];
FP_assemblies_pca = [];
FN_assemblies_pca = [];
TP_pca = [];
FP_pca = [];
FN_pca = [];

correct_nb = [];
detected_nb = [];
detected_nb_original = [];
for nb_assemblies = 1:1
    for nb_neurons = 2:iterations_nb_neurons
        for i = 1:nb_iterations
            cur_cur_nb_bins_together = 15;
            cur_nb_intervals_together = 30;
            [synthetic_data, neurons_in_assembly, activations, synthetic_data_non_zscore] = generate_synthetic_data(stimulus_data_m, sos_results_m, nb_assemblies, nb_neurons, missing_neurons, nb_intervals, cur_cur_nb_bins_together, nb_neurons_of_interest, 1, 1);

            % ica
            [ica_predicted_nbr_assemblies, ica_predicted_nbr_neurons, ica_assemblies, ica_activity, ~] = ica_assembly_detection(synthetic_data', 0);
            [ica_predicted_nbr_assemblies_original, ica_predicted_nbr_neurons_original, ica_assemblies_original, ica_activity_original, ~] = ica_assembly_detection_original(synthetic_data', 0);

            correct_nb = [correct_nb, nb_neurons];
            detected_nb = [detected_nb, ica_predicted_nbr_neurons];
            detected_nb_original = [detected_nb_original, ica_predicted_nbr_neurons_original];
        end
    end
end

