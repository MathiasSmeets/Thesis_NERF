clear;clc;close all
stimulus_data_m = load("X:\Mathias\switch_data\data_after_stimulus\after_stimulus_data_m_horridge.mat");stimulus_data_m=stimulus_data_m.after_stimulus_data_m;

[avg_param, std_param] = calculate_distribution_params(stimulus_data_m{2,1});
total_params = calculate_params(stimulus_data_m{2,1});

function [avg_param, std_param] = calculate_distribution_params(data)
total_ones = sum(data,'all');
% Get the size of the matrix
[rows, cols] = size(data);

all_params = zeros(rows*1000,1);
for i = 1:1000
% Create a new matrix with the same size filled with zeros
random_matrix = zeros(rows, cols);
% Randomly scatter ones in the new matrix
indices = randperm(rows*cols, total_ones);
[r, c] = ind2sub([rows, cols], indices);
random_matrix(sub2ind([rows, cols], r, c)) = 1;
all_params((i-1)*rows + 1:i*rows) = calculate_params(random_matrix);
end
avg_param = mean(all_params);
std_param = std(all_params);
end

function total_params = calculate_params(data)
% input: vector with rows = neurons; columns = time

%% add gaussian kernel
% Gaussian kernel parameters
kernel_size = 5; % Size of the kernel (odd number)
sigma = 1; % Standard deviation of the Gaussian kernel

% Create the Gaussian kernel
kernel = gausswin(kernel_size, sigma);

% Normalize the kernel
kernel = kernel / sum(kernel);

% Initialize a matrix to store smoothed data
smoothed_data = zeros(size(data));

% Apply convolution to each row independently
for i = 1:size(data, 1)
    smoothed_data(i, :) = conv(data(i, :), kernel, 'same');
end

%% go over each row and try to predict it

neurons_in_cluster = [];
total_params = zeros(size(smoothed_data,1),1);
for i = 1:size(smoothed_data,1)
    cur_target_row = smoothed_data(i,:)';
    cur_input = smoothed_data;
    cur_input(i, :) = [];
    cur_input = [ones(1,size(smoothed_data,2)); cur_input]';

    [parameters,~,~,~,stats] = regress(cur_target_row, cur_input);
    parameters = parameters(2:end);
    
    if stats(3) < 0.05
        neurons_in_cluster = [neurons_in_cluster i];
        if i == 1
            total_params(i+1:end) = total_params(i+1:end) + parameters;
        elseif i == length(total_params)
            total_params(1:end-1) = total_params(1:end-1) + parameters;
        else
            total_params(1:i-1) = total_params(1:i-1) + parameters(1:i-1);
            total_params(i+1:end) = total_params(i+1:end) + parameters(i:end);
        end
    end
end

end



%% problem with second method:
% if i only take one interval after stimulus, there is no such thing as
% doing this guassian kernel and measuring max, i will only have one peak
% in general along with some noise
% therefore, i have to take multiple intervals together and then this
% approach is diffeicult because it does not make sense to take multiple
% intervals and then only take one maximum