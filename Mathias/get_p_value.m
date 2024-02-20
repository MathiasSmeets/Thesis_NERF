function p_value = get_p_value(table)

test = 'chi_square';

if strcmp(test,'chi_square')
    [~, p_value] = chi2gof(table(1,:));
else
    small_n = table(1,1) + table(2,1); % total values of each column
    large_N = small_n * size(table,2); % large_N = sum small_n for each column
    % y_1j = spikes in column j
    % r1 = sum y_1j over all columns

    cur_prob = 1;
    for k = 1:size(table,2)
        cur_prob = cur_prob * nchoosek(small_n,table(1,k));
    end
    p_value = cur_prob / nchoosek(large_N,sum(table(1,:)));
end

end