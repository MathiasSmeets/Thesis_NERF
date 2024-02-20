function cc = create_cross_correlogram(candidate_templates, reference_neuron, event_neuron)

[~,nb_columns] = cellfun(@size,candidate_templates);
max_nb_columns = max(nb_columns);
cc = zeros(2,max_nb_columns*2-1);

for i = 1:length(candidate_templates)
    cur_data = candidate_templates{i};
    
    if ~isempty(find(cur_data(reference_neuron,:) == 1, 1))

        cur_xcorr = xcorr(cur_data(event_neuron,:), cur_data(reference_neuron,:));
        size_difference = size(cc,2)-size(cur_xcorr,2);
        cc(1,size_difference/2+1:end-size_difference/2) = cc(1,size_difference/2+1:end-size_difference/2) + cur_xcorr;
        inverse_cur_xcorr = zeros(1,size(cc,2));
        inverse_cur_xcorr(1,size_difference/2+1:end-size_difference/2) = cur_xcorr;
        cc(2,:) = cc(2,:) + ~round(inverse_cur_xcorr);
    end
end
cc = round(cc);

end