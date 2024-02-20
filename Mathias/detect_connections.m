function connections = detect_connections(candidate_templates)

nb_neurons = size(candidate_templates{1},1);
options = 0;
for i = 1:nb_neurons
    options = options + i;
end
p_values = zeros(1,options);
counter = 1;
connections = cell(1, options);
for i = 1:nb_neurons
    for j = i:nb_neurons 
        if i ~= j
            cc = create_cross_correlogram(candidate_templates, i, j);
            p_values(1,counter) = get_p_value(cc, 'fisher');
            if p_values(1,counter) < 0.05
                connections{counter} = [i,j];
            end
            counter = counter + 1;
        end

    end
end
connections = connections(~cellfun('isempty',connections));
end
