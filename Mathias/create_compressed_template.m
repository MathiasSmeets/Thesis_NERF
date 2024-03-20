function compressed_template = create_compressed_template(template, factor)

    % Get the size of the input matrix
    [rows, cols] = size(template);
    
    % Compute the new number of columns based on the compression factor
    new_cols = cols / factor;
        
    if factor > 1
        % Initialize the transformed matrix
        compressed_template = zeros(rows, new_cols);
        % Compress the matrix
        for i = 1:new_cols
            % Define the range of columns from the original matrix to compress
            start_col = (i - 1) * factor + 1;
            end_col = i * factor;
            
            % Compute the mean of the columns in the range
            compressed_template(:, i) = mean(template(:, start_col:end_col), 2);
        end
    else
        % Initialize the transformed matrix
        compressed_template = zeros(rows, numel(1:cols/new_cols:cols));
        % Expand the matrix
        for i = 1:rows
            % Interpolate values for the new column using linear interpolation
            compressed_template(i, :) = interp1(1:cols, template(i,:), 1:cols/new_cols:cols, 'linear');
        end
    end

compressed_template = (compressed_template./sum(compressed_template,2)).*sum(template,2);
end