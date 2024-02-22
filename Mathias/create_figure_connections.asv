function [reduced_figure_connections, indices_non_removed_rows] = create_figure_connections(connections, correlations, mouse_nb, total_nb_neurons, last_interval_connections, last_interval_data, interval_size)

nb_connections = 0;
for i = 1:total_nb_neurons
    nb_connections = nb_connections + i;
end

%last_interval = ceil(last_interval_data / interval_size);
figure_connections = zeros(nb_connections, last_interval_data);

% loop over all intervals
for i = 1:last_interval_connections
    if ~isempty(connections{mouse_nb, i})
        % loop over all clusters in this interval
        for j = 1:length(connections{mouse_nb,i})
            if ~isempty(connections{mouse_nb,i}{j})
                % loop over all connections in this cluster
                for k = 1:length(connections{mouse_nb,i}{j})
                    cur_connection = connections{mouse_nb,i}{j}{k};
                    cur_index = get_index(cur_connection, total_nb_neurons);
                    cur_correlation = correlations{mouse_nb,i}{j}{k};
                    figure_connections(cur_index,1+interval_size*(i-1):interval_size*i) = cur_correlation;
                end
            end
        end
    end
end

zero_rows = all(figure_connections == 0, 2);
reduced_figure_connections = figure_connections(~zero_rows, :);
indices_non_removed_rows = find(~zero_rows);

figure
plot(smooth(reduced_figure_connections'))

end