clear;clc;close all;

correlation_path = "X:\Mathias\switch_data\correlations";
correlation_before_m = load(fullfile(correlation_path, "thresholded_before_m.mat"));correlation_before_m = correlation_before_m.thresholded_before;
correlation_between_m = load(fullfile(correlation_path, "thresholded_between_m.mat"));correlation_between_m = correlation_between_m.thresholded_between;
correlation_after_m = load(fullfile(correlation_path, "thresholded_after_m.mat"));correlation_after_m = correlation_after_m.thresholded_after;
correlation_before_y = load(fullfile(correlation_path, "thresholded_before_y.mat"));correlation_before_y = correlation_before_y.thresholded_before;
correlation_between_y = load(fullfile(correlation_path, "thresholded_between_y.mat"));correlation_between_y = correlation_between_y.thresholded_between;
correlation_after_y = load(fullfile(correlation_path, "thresholded_after_y.mat"));correlation_after_y = correlation_after_y.thresholded_after;

%% learner mice

for i=1:9
    
end