function [synthetic_data, neurons_in_assembly, activations, synthetic_data_non_zscore] = generate_synthetic_data(reference_data, sos_results_m, nb_assemblies, nb_neurons_assembly, missing_neurons, nb_intervals, bins_together, total_neurons, random_wanted)
% create synthetic data based on statistics of reference data
% add assemblies with specified amount of neurons
% add additional noise by removing neurons from assemblies
interval_size = size(reference_data{1,1},2);
%avg total neurons of interest: 9.3514

%% calculate average over reference data
cur_counter = zeros(total_neurons,interval_size);
interval_counter = 0;
for i = 1%:size(reference_data,1)
    for j = 1:size(reference_data,2)
        if ~isempty(reference_data{i,j})
            if size(reference_data{i,j},1) >= total_neurons
                cur_counter = cur_counter + reference_data{i,j}(1:total_neurons,:);
                interval_counter = interval_counter + 1;
            end
        end
    end
end
lambda = cur_counter / interval_counter;

%% create background noise
synthetic_data = zeros(total_neurons,nb_intervals*interval_size);
for i = 1:nb_intervals
    synthetic_data(:,(i-1)*interval_size+1:i*interval_size) = poissrnd(lambda);
end

%% add assemblies

activations = zeros(1,nb_assemblies);
for i = 1:nb_assemblies
    neurons_in_assembly = sort(randperm(total_neurons, nb_neurons_assembly));
    index = randperm(55-12,1) + 11;
    activations(i) = index;
    random_factor = randperm(11,length(neurons_in_assembly))-6;
    for j = 1:nb_intervals
        if random_wanted
            for k = 1:length(neurons_in_assembly)
                synthetic_data(neurons_in_assembly(k), index + (j-1)*interval_size + random_factor(k)) = 1;
            end
        else
            random_factor = 0;
            synthetic_data(neurons_in_assembly, index + (j-1)*interval_size + random_factor) = 1;
        end
    end
end

%% remove neurons dependent on missing_neurons variable

if missing_neurons ~= 0
    assemblies_to_remove = randperm(nb_intervals, missing_neurons);
    for i = 1:missing_neurons
        index_neuron_to_remove = randperm(length(neurons_in_assembly),1);
        synthetic_data(neurons_in_assembly(index_neuron_to_remove),assemblies_to_remove) = 0;
    end
end

%% put bins together

indices = ceil((1:interval_size*nb_intervals)/bins_together);
new_synthetic_data = zeros(total_neurons, nb_intervals*interval_size / bins_together);
for i = 1:size(synthetic_data,1)
    new_synthetic_data(i,:) = accumarray(indices',synthetic_data(i,:)',[],@sum)';
end
synthetic_data = new_synthetic_data;

%% remove zero rows
for jj = 1:total_neurons
    if all(synthetic_data(jj,:)==synthetic_data(jj,1))
        synthetic_data(jj,randperm((nb_intervals*interval_size)/bins_together,1)) = 1;
    end
end


%% calculate zscores

synthetic_data_non_zscore = synthetic_data;
cur_std = sos_results_m{1,2};
cur_std(cur_std==0) = 0.01;
synthetic_data = (synthetic_data - sos_results_m{1,1}(1:total_neurons,:)) ./ cur_std(1:total_neurons,:);

end

