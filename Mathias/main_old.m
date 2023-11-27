clear;clc;close all;

total_frM = load("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\frM_total.mat");
total_frM = struct2array(total_frM);
total_frY = load("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\frY_total.mat");
total_frY = struct2array(total_frY);

[m_data, y_data] = assemble_data(total_frM, total_frY);



%% determine features

% total average
average_m = mean(m_data, 2);
average_y = mean(y_data, 2);

% total variance
variance_m = var(m_data, 0, 2);
variance_y = var(y_data, 0, 2);

% moving average
moving_averaged_m = movmean(m_data, 1000, 2); % fs = 0.01s --> over 10s
moving_averaged_y = movmean(y_data, 1000, 2);

% variance
%moving_variance_m = movvar(m_data, 1000, 0, 2); not enough memory to do all
%moving_variance_y = movvar(y_data, 1000, 0, 2);

% freq analysis on moving average
%fft_ma_m = fft(moving_averaged_m, [], 2);
%fft_ma_y = fft(moving_averaged_y, [], 2);

% dimensionality reduction
moving_averaged_m = downsample(moving_averaged_m',10);
moving_averaged_y = downsample(moving_averaged_y',10);
moving_averaged_m = moving_averaged_m';
moving_averaged_y = moving_averaged_y';








