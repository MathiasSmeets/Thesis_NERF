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
            [synthetic_data, neurons_in_assembly, activations, synthetic_data_non_zscore] = generate_synthetic_data(stimulus_data_m, sos_results_m, nb_assemblies, nb_neurons, missing_neurons, cur_nb_intervals_together, cur_cur_nb_bins_together, nb_neurons_of_interest, 1, 1);

            % ica
            [ica_predicted_nbr_assemblies, ica_predicted_nbr_neurons, ica_assemblies, ica_activity, ~] = ica_assembly_detection(synthetic_data', 0);
            [ica_predicted_nbr_assemblies_original, ica_predicted_nbr_neurons_original, ica_assemblies_original, ica_activity_original, ~] = ica_assembly_detection_original(synthetic_data', 0);

            correct_nb = [correct_nb, nb_neurons];
            detected_nb = [detected_nb, ica_predicted_nbr_neurons];
            detected_nb_original = [detected_nb_original, ica_predicted_nbr_neurons_original];


        end
    end
end

figure;hold on
actual = zeros(10,iterations_nb_neurons-1);
original = zeros(10,iterations_nb_neurons-1);
for i = 1:iterations_nb_neurons-1
    for k = 0:10
        actual(k+1,i) = numel(find(detected_nb((i-1)*nb_iterations+1:i*nb_iterations)==k));
        original(k+1,i) = numel(find(detected_nb_original((i-1)*nb_iterations+1:i*nb_iterations)==k));
    end
end

figure
boxplot([detected_nb_original(1:100)', detected_nb(1:100)',detected_nb_original(101:200)', detected_nb(101:200)',detected_nb_original(201:300)', detected_nb(201:300)',detected_nb_original(301:400)', detected_nb(301:400)'])

figure
boxplot([detected_nb_original(1:100)', detected_nb(1:100)',detected_nb_original(101:200)', detected_nb(101:200)',detected_nb_original(201:300)', detected_nb(201:300)',detected_nb_original(301:400)', detected_nb(301:400)'])
%
total_matrix = zeros(size(actual,1), size(actual,2)*2);
for i = 1:size(actual,2)
    total_matrix(:,1+(i-1)*2) = actual(:,i);
    total_matrix(:,2+(i-1)*2) = original(:,i);
end

% Define the positions for the bars
numBars = size(actual, 2); % Assuming both matrices have the same number of columns
groupWidth = 0.8; % Width of each group
barWidth = groupWidth / 2; % Width of each bar
spaceWidth = groupWidth / 100; % Width of space between groups
positions = 1:numBars; % Positions for bars
positions1 = positions - barWidth/2 - spaceWidth;
positions2 = positions + barWidth/2 + spaceWidth;
actual_positions = sort([positions1, positions2]);

% Create the grouped stacked bar plot
figure;
h=bar(actual_positions, total_matrix, 'stacked');
hold on;

% Customize the plot
xlabel('Group');
ylabel('Value');
title('Grouped Stacked Bar Plot');
legend('Matrix1', 'Matrix2');

% get numbers inside bar
% Compute bar segment centers
xbarCnt = vertcat(h.XEndPoints); 
ybarTop = vertcat(h.YEndPoints);
ybarCnt = ybarTop - total_matrix/2; 

yp = total_matrix / nb_iterations;
% Create text strings
txt = compose('%.1f%%',yp);

% Add text 
th = text(xbarCnt(:), ybarCnt(:), txt(:), ...
    'HorizontalAlignment', 'center', ....
    'VerticalAlignment', 'middle', ...
    'Color', 'w',....
    'FontSize', 8);

legend({"0 neurons", "1 neurons", "2 neurons", "3 neurons", "4 neurons", "5 neurons", "6 neurons"},"Location","Eastoutside")

h(4).FaceColor = [0.4660 0.6740 0.1880];
h(5).FaceColor = [0.6350 0.0780 0.1840];
h(6).FaceColor = [0.4940 0.1840 0.5560];
h(7).FaceColor = [0 0 0];
