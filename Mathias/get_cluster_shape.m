clear; clc; close all;

if ispc
    volume_base = '\\nerffs13';
    volume_base2 = '\\nerffs17';
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/';
    volume_base2 = '/mnt/takeokalab/';
end

path_to_code = "takeokalabwip2023/Mathias/switch_data/data_after_stimulus";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
%path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

stimulus_data_m = load(fullfile(volume_base2, path_to_code,"after_stimulus_data_m_horridge.mat"));
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_y = load(fullfile(volume_base2, path_to_code,"after_stimulus_data_y_horridge.mat"));
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

neurons_of_interest_path = "\takeokalabwip2023\Mathias\switch_data\neurons_of_interest";
neurons_of_interest_m = load(fullfile(volume_base2, neurons_of_interest_path, "neurons_of_interest_horridge_m.mat"));
neurons_of_interest_m = neurons_of_interest_m.output_m;
inhibited_neurons_m = load(fullfile(volume_base2, neurons_of_interest_path, "inhibited_horridge_m.mat"));
inhibited_neurons_m = inhibited_neurons_m.inhibited_m;

folder = fileparts(which("clusters_cpd.m"));
addpath(genpath(folder))

%% load in ica data

path_to_data = "takeokalabwip2023/Mathias/switch_data/clusters";
ica_activity = load(fullfile(volume_base2, path_to_data,"activity_horridge_m.mat")); ica_activity = ica_activity.total_activity;
ica_nb_neurons = load(fullfile(volume_base2, path_to_data,"nb_neurons_horridge_m.mat")); ica_nb_neurons = ica_nb_neurons.total_nb_neurons;
ica_assemblies = load(fullfile(volume_base2, path_to_data,"assemblies_horridge_m.mat")); ica_assemblies = ica_assemblies.total_assemblies;
ica_nb_assemblies = load(fullfile(volume_base2, path_to_data,"nb_assemblies_horridge_m.mat")); ica_nb_assemblies = ica_nb_assemblies.total_nb_assemblies;
ica_data = load(fullfile(volume_base2, path_to_data,"data_horridge_m.mat")); ica_data = ica_data.total_data;
ica_neurons_of_interest = load(fullfile(volume_base2, path_to_data,"neurons_of_interest_horridge_m.mat")); ica_neurons_of_interest = ica_neurons_of_interest.total_neurons_of_interest;

%% initialize data
interval_step = 30; % depends on how data was created
bins_together = 15;
interval_size = 70;
actual_connections = cell(size(ica_neurons_of_interest));
correlation_strengths = cell(size(ica_neurons_of_interest));

%% get cluster shape for each assembly

orders = cell(size(ica_activity));
% loop over each mouse
for k = 1:size(orders,1)
    % loop over each interval
    disp(k)
    for i = 1:size(orders,2)
        cur_assembly = ica_assemblies{k,i};
        if isempty(cur_assembly)
            actual_connections{k,i} = [];
        else
            cur_raw_connections = cell(size(cur_assembly));
            cur_actual_connections = cell(size(cur_assembly));
            cur_correlation_strength = cell(size(cur_assembly));
            for j = 1:length(cur_assembly)
                curcur_assembly = ica_neurons_of_interest{k,i}(cur_assembly{j});
                cur_data = ica_data{k,i};
                if length(ica_activity{k,i}(:,j)) >= 3
                    [pks, locs] = findpeaks(abs(ica_activity{k,i}(:,j)),"NPeaks",interval_step,"MinPeakHeight",0.4*max(abs(ica_activity{k,i}(:,j))));
                    candidate_templates = cell(1,length(locs));
                    TOTAL = [];
                    for jj = 1:length(locs)
                        raw_data_index = (i-1)*interval_step + ceil((locs(jj))/ceil(interval_size/bins_together));
                        position_in_data = mod(locs(jj)-1,ceil(interval_size/bins_together))*bins_together+1;
                        cur_raw_data = stimulus_data_m{k,raw_data_index};
                        cur_assembly_data = cur_raw_data(curcur_assembly, position_in_data:min(position_in_data+bins_together-1, interval_size));

                        % find the indices of the first and last non-zero elements in each row
                        first_nonzero = find(any(cur_assembly_data, 1), 1, 'first');
                        last_nonzero = find(any(cur_assembly_data, 1), 1, 'last');
                        candidate_templates{1,jj} = cur_assembly_data(:, first_nonzero:last_nonzero);

                    end
                    % find the order of the neurons firing
                    %orders{k,i} = find_order_neurons(candidate_templates);

                    % detect connections
                    if ~isempty(candidate_templates)
                        [cur_raw_connections, curcur_correlation_strength] = detect_connections(candidate_templates, bins_together);
                        curcur_actual_connections = cell(size(cur_raw_connections));
                        for ii = 1:length(cur_raw_connections)
                            curcur_actual_connections{ii} = curcur_assembly(cur_raw_connections{ii});
                        end
                    else
                        curcur_actual_connections = [];
                        curcur_correlation_strength = [];
                    end
                else
                    curcur_actual_connections = [];
                    curcur_correlation_strength = [];
                end
                cur_actual_connections{j} = curcur_actual_connections;
                cur_correlation_strength{j} = curcur_correlation_strength;
            end
            actual_connections{k,i} = cur_actual_connections;
            correlation_strengths{k,i} = cur_correlation_strength;
        end
    end
end

save("X:\Mathias\switch_data\connections\connections_horridge_m.mat","actual_connections", "-v7.3");
save("X:\Mathias\switch_data\connections\correlations_horridge_m.mat","correlation_strengths", "-v7.3");


%% create_figure

% for mouse_nb = 1:size(stimulus_data_m)
%     total_nb_neurons = size(stimulus_data_m{mouse_nb,1},1);
%     counter = size(stimulus_data_m,2);
    % while isempty(stimulus_data_m{mouse_nb,counter})
    %     counter = counter - 1;
    %     last_interval_data = counter;
    % end
%     while isempty(correlation_strengths{mouse_nb,counter}) || isempty(correlation_strengths{mouse_nb,counter}{1,1})
%         counter = counter - 1;
%         last_interval_correlations = counter;
%     end
%     [reduced_figure_connections, indices_non_removed_rows] = create_figure_connections(actual_connections, correlation_strengths, mouse_nb, total_nb_neurons, last_interval_correlations, last_interval_data, interval_step);
% end

%% create figure

for i = 1:size(stimulus_data_m,1)
    % get last interval
    counter = size(stimulus_data_m,2);
    last_interval_data = counter;
    while isempty(stimulus_data_m{i,counter})
        counter = counter - 1;
        last_interval_data = counter;
    end
    % create figure_data
    figure_data = zeros(size(stimulus_data_m{i,1},1), ceil(last_interval_data/interval_step));
    for j = 1:ceil(last_interval_data/interval_step)
        figure_data(:,j) = get_labels_from_connections(actual_connections{i,j}, size(stimulus_data_m{i,1},1));
    end
    figure;imagesc(figure_data)
end
