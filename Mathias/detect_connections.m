function connections = detect_connections(candidate_templates)

connections = 0;
nb_neurons = size(candidate_templates{1},1);
p_values = zeros(1,factorial(nb_neurons));
counter = 1;
most_spikes = 0;
for i = 1:nb_neurons
    for j = i:nb_neurons 
        if i ~= j
            cc = create_cross_correlogram(candidate_templates, i, j);
            p_values(1,counter) = get_p_value(cc);
            counter = counter + 1;
            if max(sum(cc)) > most_spikes
                most_spikes = max(sum(cc));
            end
        end

    end
end

[~,nb_columns] = cellfun(@size,candidate_templates);
max_nb_columns = max(nb_columns);
array_p_values = get_array_p_values(max_nb_columns, most_spikes);

end
