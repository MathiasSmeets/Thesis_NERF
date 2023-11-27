function [m_data, y_data] = assemble_data(total_frM, total_frY)
% this function will concatenate the early, mid and late periods into a
% matrix that is easier to work with

m_data = zeros(size(total_frM.zscore_early_full,1),size(total_frM.zscore_early_full{1,1},2) + size(total_frM.zscore_mid_full{1,1},2) + size(total_frM.zscore_late_full{1,1},2));
y_data = zeros(size(total_frY.zscore_early_full,1),size(total_frY.zscore_early_full{1,1},2) + size(total_frY.zscore_mid_full{1,1},2) + size(total_frY.zscore_late_full{1,1},2));

for i = 1:size(total_frY.zscore_early_full,1)

    temp_y_data = [total_frY.zscore_early_full{i,1}, total_frY.zscore_mid_full{i,1}, total_frY.zscore_late_full{i,1}];
    if size(temp_y_data,2) > size(y_data,2)
        temp_y_data = temp_y_data(1:size(y_data,2));
    elseif size(temp_y_data,2) < size(y_data,2)
        extra_zeros_y = zeros(1,size(y_data,2) - size(temp_y_data,2));
        temp_y_data = [temp_y_data , extra_zeros_y];
    end

    y_data(i,:) = temp_y_data;
end

for i = 1:size(total_frM.zscore_early_full,1)

    temp_m_data = [total_frM.zscore_early_full{i,1}, total_frM.zscore_mid_full{i,1}, total_frM.zscore_late_full{i,1}];
    if size(temp_m_data,2) > size(m_data,2)
        temp_m_data = temp_m_data(1:size(m_data,2));
    elseif size(temp_m_data,2) < size(m_data,2)
        extra_zeros_m = zeros(1,size(m_data,2) - size(temp_m_data,2));
        temp_m_data = [temp_m_data , extra_zeros_m];
    end

    m_data(i,:) = temp_m_data;
end

end