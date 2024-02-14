function cur_total_mouse_zscore_delayed = create_convolutive_data(cur_total_mouse_zscore, largest_delay)


if largest_delay == 0
    cur_total_mouse_zscore_delayed = cur_total_mouse_zscore;
else
    cur_total_mouse_zscore_delayed = zeros(size(cur_total_mouse_zscore,1)*(largest_delay+1),size(cur_total_mouse_zscore,2));
    for i = 1:largest_delay+1
        cur_rows_delayed = i:largest_delay+1:size(cur_total_mouse_zscore,1)*(largest_delay+1);
        cur_total_mouse_zscore_delayed(cur_rows_delayed,i:end) = cur_total_mouse_zscore(:,1:end-(i-1));
    end
end

end