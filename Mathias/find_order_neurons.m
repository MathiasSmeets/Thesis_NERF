function order = find_order_neurons(candidate_templates)

amount_of_neurons = size(candidate_templates{1,1},1);
delay = zeros(size(candidate_templates{1,1},1));
order = zeros(size(candidate_templates{1,1},1));
counter = 1;

for i = 1:size(candidate_templates,2)
    cur_template = candidate_templates{1,i};
    % check if all rows have a spike
    
    if all(any(cur_template,2))
        % for these, calculate the average order and delay
        % for order, you first remove the columns that are all zeros
        nonzero_columns = any(cur_template, 1);
        cur_template_nonzero = cur_template(:,nonzero_columns);

        for j = 1:size(cur_template,1)
            delay(j,counter) = find(cur_template(j,:) == 1,1);
            order(j,counter) = find(cur_template_nonzero(j,:)==1,1);
        end
        counter = counter + 1;
    end
end
orders = mean(order,2);
[order_values, indices] = sort(orders);
order = zeros(size(order_values));
order(indices) = 1:numel(orders);

delay_values = mean(delay,2);
end