function neurons_of_interest = get_neurons_of_interest(data, output_m, inhibited_neurons_m, neuron_counter)

    % get neurons of interest
    cur_neurons_of_interest = [];
    for ii = 1:size(output_m,2)-1 % do not use this for inhibited neurons
        for jj = 1:size(output_m{1,ii}{1,1},1)
            if output_m{1,ii}{1,1}{jj,1} >= neuron_counter && output_m{1,ii}{1,1}{jj,1} < neuron_counter + size(data,1)
                cur_neurons_of_interest = [cur_neurons_of_interest, output_m{1,ii}{1,1}{jj,1}];
            end
        end
    end
    inhibited_neurons_of_interest = inhibited_neurons_m(inhibited_neurons_m >= neuron_counter & inhibited_neurons_m < neuron_counter + size(data,1));
    cur_neurons_of_interest = [cur_neurons_of_interest, inhibited_neurons_of_interest];
    cur_neurons_of_interest = unique(cur_neurons_of_interest);
    neurons_of_interest = cur_neurons_of_interest - neuron_counter + 1;


end