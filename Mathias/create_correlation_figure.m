clear;clc;close all;

correlation_path = "X:\Mathias\switch_data\correlations";
correlation_before_m = load(fullfile(correlation_path, "correlation_before_smoothed_width3.mat"));correlation_before_m = correlation_before_m.avg_cor_before;
correlation_between_m = load(fullfile(correlation_path, "correlation_between_smoothed_width3.mat"));correlation_between_m = correlation_between_m.avg_cor_between;
correlation_after_m = load(fullfile(correlation_path, "correlation_after_smoothed_width3.mat"));correlation_after_m = correlation_after_m.avg_cor_after;
correlation_before_y = load(fullfile(correlation_path, "correlation_before_smoothed_width3_y.mat"));correlation_before_y = correlation_before_y.avg_adj_cur_correlation_before;
correlation_between_y = load(fullfile(correlation_path, "correlation_between_smoothed_width3_y.mat"));correlation_between_y = correlation_between_y.avg_adj_cur_correlation_between;
correlation_after_y = load(fullfile(correlation_path, "correlation_after_smoothed_width3_y.mat"));correlation_after_y = correlation_after_y.avg_adj_cur_correlation_after;

