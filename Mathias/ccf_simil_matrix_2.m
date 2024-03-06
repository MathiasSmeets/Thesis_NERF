function ccf_simil_matrix_2(candidate_templates)

%% Create aussian kernel
% Gaussian kernel parameters
kernel_size = 5; % Size of the kernel (odd number)
sigma = 1; % Standard deviation of the Gaussian kernel

% Create the Gaussian kernel
kernel = gausswin(kernel_size, sigma);

% Normalize the kernel
kernel = kernel / sum(kernel);

%% Cross-correlation between neurons

for i = 1:length(candidate_templates)
    cur_candidate_template = candidate_templates{i};

    % Initialize a matrix to store smoothed data
    smoothed_data = zeros(size(cur_candidate_template));

    % Apply convolution to each row independently
    for j = 1:size(cur_candidate_template, 1)
        smoothed_data(j, :) = conv(cur_candidate_template(j, :), kernel, 'same');
    end
    
    % do cross-correlation
    for j = 1:size(cur_candidate_template,1)-1
        for k = i+1:size(cur_candidate_template,1)
            if ~all(smoothed_data(j, :)==0,2) && ~all(smoothed_data(k, :)==0,2)
                cur_cor = xcorr(smoothed_data(j, :), smoothed_data(k, :));
            end
        end
    end

    % detect spike in cross_correlation

end

%% Put in similarity matrix

%% Perform clustering in similarity matrix

end