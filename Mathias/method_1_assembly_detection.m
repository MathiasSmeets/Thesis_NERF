function [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = method_2_assembly_detection(data)
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