clear;clc;close all;

path_to_data = "X:\Mathias\switch_data\data_after_stimulus";
before_data_m = load("X:\Mathias\switch_data\data_after_stimulus\before_data_m.mat"); before_data_m = before_data_m.before_data;
after_data_m = load("X:\Mathias\switch_data\data_after_stimulus\after_data_m.mat"); after_data_m = after_data_m.after_data;
connections_m = load("X:\Mathias\switch_data\connections\connections_horridge_m.mat"); connections_m = connections_m.actual_connections;
stimulus_data_m = load("X:\Mathias\switch_data\data_after_stimulus\after_stimulus_data_m_horridge.mat"); stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
intervals_together = 30;
delay_size = 10;
total_neurons = max(before_data_m(:,1));

%% get last intervals and remove stimulus data
last_connection = zeros(total_neurons,1);
for i = 1:total_neurons
    last_interval = size(stimulus_data_m,2);
    while isempty(stimulus_data_m{i,last_interval})
        last_interval = last_interval-1;
    end
    last_connection(i) = ceil(last_interval/intervals_together);
end
clearvars stimulus_data_m

avg_correlations = zeros(total_neurons,2);
for i = 1:total_neurons
    %% initialize variables
    cur_before_data = before_data_m(before_data_m(:,1) == i,:);
    cur_before_data = cur_before_data(:,2:end);
    cur_after_data = after_data_m(after_data_m(:,1) == i,:);
    cur_after_data = cur_after_data(:,2:end);
    cur_nb_neurons = size(cur_before_data,1);
    nb_connections = 0;
    for j = 1:cur_nb_neurons-1
        nb_connections = nb_connections + j;
    end


    %% get connections
    connection_count = zeros(nb_connections,1);
    unique_connections = {};
    for j = 1:last_connection(i)
        for k = 1:numel(connections_m{i,j})
            for l = 1:numel(connections_m{i,j}{k})
                cur_connection = connections_m{i,j}{k}{l};
                cur_index = get_index(cur_connection, cur_nb_neurons);
                connection_count(cur_index) = connection_count(cur_index) + 1;
                if ~any(cellfun(@(x) isequal(x, cur_connection), unique_connections))
                    unique_connections{end+1} = cur_connection;
                end
            end
        end
    end
    connection_count = connection_count / last_connection(i);
    unique_connections = cell2mat(unique_connections.');
    unique_neurons = unique(reshape(unique_connections,[size(unique_connections,1)*2,1]));
    correlations = zeros(size(unique_connections,1),1);
    %% for each connection, get correlation between these neurons in pre and rest
    % get correlation for each connection
    candidate_templates_before = {};
    candidate_templates_after = {};
    for k = 1:size(unique_connections,1)
        % for each one in one of those neurons, add this to my candidate templates
        cur_connection = unique_connections(k,:);
        [~, spikes_col_before] = find(cur_before_data(cur_connection,:));
        [~, spikes_col_after] = find(cur_after_data(cur_connection,:));
        % for each column, create a new candidate template
        previous_col_end_before = 0;
        for j = 1:numel(spikes_col_before)
            col_start_before = max(spikes_col_before(j)-delay_size,1);
            col_end_before = min(spikes_col_before(j)+delay_size, size(cur_before_data,2));
            % make sure this spike was not in previous template
            if col_start_before > previous_col_end_before
                candidate_templates_before{end+1} = cur_before_data(cur_connection,[col_start_before:col_end_before]);
            end
            previous_col_end_before = col_end_before;
        end
        % repeat for after data
        previous_col_end_after = 0;
        for j = 1:numel(spikes_col_after)
            col_start_after = max(spikes_col_after(j)-delay_size,1);
            col_end_after = min(spikes_col_after(j)+delay_size,size(cur_after_data,2));
            if col_start_after > previous_col_end_after
                candidate_templates_after{end+1} = cur_after_data(cur_connection,[col_start_after:col_end_after]);
            end
            previous_col_end_after = col_end_after;
        end
        % create cc
        if ~isempty(candidate_templates_before)
            cc = create_cross_correlogram(candidate_templates_before, 1, 2, 1+2*delay_size);
            cc_reverse = create_cross_correlogram(candidate_templates_before, 2, 1, 1+2*delay_size);
            correlations(k,1) = min(get_correlation(cc, 0), get_correlation(cc_reverse, 0));
        end
        if ~isempty(candidate_templates_after)
            cc = create_cross_correlogram(candidate_templates_after, 1, 2, 1+2*delay_size);
            cc_reverse = create_cross_correlogram(candidate_templates_after, 2, 1, 1+2*delay_size);
            correlations(k,2) = min(get_correlation(cc, 0), get_correlation(cc_reverse, 0));
        end

    end

    %% compare
    avg_correlations(i,:) = mean(correlations);
end