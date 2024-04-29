function [detected_assemblies,nb_assembly_FP, nb_assembly_FN, nb_neurons_TP, nb_neurons_FP, nb_neurons_FN,nb_neurons_TN] = get_performance_clustering(true_assemblies, found_assemblies)
% 
detected_assemblies = numel(found_assemblies);
nb_assembly_FP = 0;
nb_assembly_FN = 0;

nb_neurons_TP = 0;
nb_neurons_FP = 0;
nb_neurons_FN = 0;
nb_neurons_TN = 0;

if isa(found_assemblies, 'double')
    if found_assemblies == 0
        cur_TP = 0;
        cur_FN = numel(true_assemblies);
        cur_FP = 0;
        cur_TN = 10 - cur_FP - cur_FN - cur_TP;

        nb_neurons_TP = cur_TP + nb_neurons_TP;
        nb_neurons_FN = cur_FN + nb_neurons_FN;
        nb_neurons_FP = cur_FP + nb_neurons_FP;
        nb_neurons_TN = cur_TN + nb_neurons_TN;
    end
else
    for i = 1:numel(found_assemblies)
        cur_assembly = found_assemblies{i};

        cur_TP = numel(intersect(true_assemblies, cur_assembly));
        cur_FN = numel(setdiff(true_assemblies, intersect(true_assemblies, cur_assembly)));
        cur_FP = numel(setdiff(cur_assembly, intersect(true_assemblies, cur_assembly)));
        cur_TN = 10 - cur_FP - cur_FN - cur_TP;

        nb_neurons_TP = cur_TP + nb_neurons_TP;
        nb_neurons_FN = cur_FN + nb_neurons_FN;
        nb_neurons_FP = cur_FP + nb_neurons_FP;
        nb_neurons_TN = cur_TN + nb_neurons_TN;
    end
end


% nb_assembly_TP = 0;
% nb_neurons_TP = 0;
% nb_found_assemblies = length(found_assemblies);
% 
% %% check if entire assemblies are correct
% if iscell(found_assemblies)
%     for i = 1:nb_found_assemblies
%         cur_assembly_correct = false;
%         cur_neurons_same = 0;
%         for j = 1:length(true_assemblies)
%             if isequal(sort(found_assemblies{i}), sort(true_assemblies{j})) % entire assembly is correct
%                 cur_assembly_correct = true;
%             else                                                % check for individual neurons
%                 cur_cur_neurons_same = length(intersect(found_assemblies{i}, true_assemblies{j}));
%                 if cur_cur_neurons_same > cur_neurons_same
%                     cur_neurons_same = cur_cur_neurons_same;
%                 end
%             end
%         end
%         if cur_assembly_correct
%             nb_assembly_TP = nb_assembly_TP + 1;
%             nb_neurons_TP = nb_neurons_TP + length(found_assemblies{i});
%         else
%             nb_neurons_TP = nb_neurons_TP + cur_neurons_same;
%         end
%     end    
% end
% 
% %% calculate neurons FN, FP, TN
% total_true_neurons = 0;
% for i = 1:length(true_assemblies)
%     total_true_neurons = total_true_neurons + length(true_assemblies{i});
% end
% total_detected_neurons = 0;
% for i = 1:length(found_assemblies)
%     if iscell(found_assemblies)
%         total_detected_neurons = total_detected_neurons + length(found_assemblies{i});
%     end
% end
% 
% nb_neurons_FN = total_true_neurons - nb_neurons_TP;
% nb_neurons_FP = total_detected_neurons - nb_neurons_TP;
% nb_assembly_FN = length(true_assemblies) - nb_assembly_TP;
% nb_assembly_FP = nb_found_assemblies - nb_assembly_TP;
% nb_neurons_TN = 10 - nb_neurons_FP - nb_neurons_FN - nb_neurons_TP;


