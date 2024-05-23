clear;clc;close all;

m = load("X:\Mathias\switch_data\correlations\results_m.mat"); m = m.normalized_averages;
np2 = load("X:\Mathias\switch_data\correlations\results_np2.mat"); np2 = np2.normalized_averages;

nm = load("X:\Mathias\switch_data\data_after_stimulus\last_interval_data.mat"); nm = nm.last_interval_data;
nnp2 = load("X:\Mathias\switch_data\data_after_stimulus\last_interval_data_np2.mat"); nnp2 = nnp2.last_interval_data;


replay = [m(:,3)-m(:,1);np2(:,3)-np2(:,1)];

last_interval = [nm, nnp2];

scatter(replay, last_interval, "filled")
ylim([0 1100])
