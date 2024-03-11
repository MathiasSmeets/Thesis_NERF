function neuron_labels = get_labels_from_connections(connections, nb_neurons)

% example connections:
%connections = { [1,2]; [2,3]; [4,5]; [6,7]; [1,6]};
% will return neurons 1,2,3,6,7 & 4,5 in clusters

original_connections = connections;
connections = cell(1,1);
for i = 1:numel(original_connections)
    connections{1} = [connections{1}, original_connections{i}];
end
connections = connections{1};

% Initialize clusters
clusters = {};

for i = 1:numel(connections)
    connection = connections{i};
    found = false;
    
    % Check if any cluster contains one of the neurons in the current connection
    for j = 1:numel(clusters)
        cluster = clusters{j};
        if any(ismember(cluster, connection))
            % Add both neurons of the connection to the cluster
            clusters{j} = unique([cluster, connection]);
            found = true;
            break;
        end
    end
    
    % If not found in any existing cluster, create a new cluster with both neurons
    if ~found
        clusters{end+1} = connection;
    end
end

% Combine clusters if necessary
combined_clusters = {};
for i = 1:numel(clusters)
    current_cluster = clusters{i};
    combined = false;
    for j = 1:numel(clusters)
        if i ~= j
            other_cluster = clusters{j};
            if any(ismember(current_cluster, other_cluster))
                new_cluster = unique([current_cluster, other_cluster]);
                % check if already new cluster is already in combined_clusters
                already_present = false;
                for k = 1:numel(combined_clusters)
                    if any(ismember(new_cluster, combined_clusters{k}))
                        already_present = true;
                    end
                end
                if ~already_present
                    combined_clusters{end+1} = unique([current_cluster, other_cluster]);
                end
                combined = true;
                break;
            end
        end
    end
    if ~combined
        combined_clusters{end+1} = current_cluster;
    end
end

% Assign labels to neurons based on clusters
neuron_labels = zeros(1, nb_neurons);
for i = 1:numel(combined_clusters)
    cluster = combined_clusters{i};
    label = i;
    for j = 1:numel(cluster)
        neuron_labels(cluster(j)) = label;
    end
end


end