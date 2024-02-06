function [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = cpd_assembly_detection(tensor, assembly_ica)

predicted_nbr_assemblies = size(assembly_ica,1);
predicted_nbr_neurons = 0;

if predicted_nbr_assemblies ~= 0
    U = cpd(tensor, predicted_nbr_assemblies+2);
    assemblies = cell(predicted_nbr_assemblies, 1);
    if predicted_nbr_assemblies == 1
        predicted_nbr_neurons = length(assembly_ica{1,1});
        [~, cur_assembly] = maxk(abs(U{1}), predicted_nbr_neurons);
        [~, column_index] = max(sum(abs(U{3})));
        assemblies{1,1} = unique(cur_assembly(:,column_index)');
    else
        for i = 1:predicted_nbr_assemblies
            cur_neurons = length(assembly_ica{i,1});
            [~, cur_assembly] = maxk(abs(U{1}), cur_neurons);
            [~, column_index] = max(sum(abs(U{3})));
            assemblies{i,1} = unique(cur_assembly(:,column_index)');
            predicted_nbr_neurons = predicted_nbr_neurons + cur_neurons;
        end
    end
else
    predicted_nbr_neurons = 0;
    assemblies = 0;
end
activity = [];
end