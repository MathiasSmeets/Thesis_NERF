function [connections, correlations] = detect_connections(candidate_templates, bins_together)

nb_neurons = size(candidate_templates{1},1);
options = 0;
for i = 1:nb_neurons-1
    options = options + i;
end
counter = 1;
total_candidates = length(candidate_templates);
connections = cell(1, options);
correlations = cell(1, options);
for i = 1:nb_neurons
    for j = i:nb_neurons 
        if i ~= j
            cc = create_cross_correlogram(candidate_templates, i, j, bins_together);
            cc_reverse = create_cross_correlogram(candidate_templates, j, i, bins_together);
            if isempty(find(cc(1,:), 1)) || isempty(find(cc_reverse(1,:), 1))
                cur_p_value = 1;
                cur_reverse_p_value = 1;
            else
                cur_p_value = get_p_value(cc, 'fisher');
                cur_reverse_p_value = get_p_value(cc_reverse, 'fisher');
            end
            if cur_p_value < 0.05 && cur_reverse_p_value < 0.05
                connections{counter} = [i,j];
                correlations{counter} = min(get_correlation(cc, total_candidates), get_correlation(cc_reverse, total_candidates));
            end
            counter = counter + 1;
        end

    end
end
connections = connections(~cellfun('isempty',connections));
correlations = correlations(~cellfun('isempty', correlations));
end
