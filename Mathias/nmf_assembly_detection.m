function [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = nmf_assembly_detection(data, assembly_ica, iterations)

% remove neurons with no activity
neurons_to_remove = all(diff(data,[],2) == 0,2);
neurons_to_keep = find(~neurons_to_remove);
data = data(neurons_to_keep,:);

% start function
predicted_nbr_assemblies = size(assembly_ica,1);
if ~iscell(assembly_ica)
    predicted_nbr_assemblies = 0;
end
predicted_nbr_neurons = 0;

if predicted_nbr_assemblies ~= 0
    [W, H, ~] = best_nmf(data, iterations);
    assemblies = cell(predicted_nbr_assemblies, 1);

    if predicted_nbr_assemblies == 1
        predicted_nbr_neurons = length(assembly_ica{1,1});
        [~, cur_assembly] = maxk(abs(W), predicted_nbr_neurons);
        [~, column_index] = max(sum(abs(H),2),[],2);
        assemblies{1,1} = neurons_to_keep(unique(cur_assembly(:,column_index))')';
    else
        for i = 1:predicted_nbr_assemblies
            cur_neurons = length(assembly_ica{i,1});
            [~, cur_assembly] = maxk(abs(W), cur_neurons);
            [~, column_index] = max(sum(abs(H),2),[],2);
            assemblies{i,1} = neurons_to_keep(unique(cur_assembly(:,column_index))');
            predicted_nbr_neurons = predicted_nbr_neurons + cur_neurons;
        end
    end

    if length(assemblies) > 1
        for i = 1:length(assemblies)-1
            if isequal(assemblies{i},assemblies{i+1})
                assemblies{i} = {};
            end
        end
    end
    assemblies = assemblies(~cellfun('isempty',assemblies));

else
    predicted_nbr_neurons = 0;
    assemblies = 0;
end
activity = [];
end
