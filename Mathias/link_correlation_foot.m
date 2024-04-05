clear;close all;clc;


foot_positions_m = zeros(1,9);
foot_positions_y = zeros(1,9);

threshold = 3;

% mouse 2
cur_foot = load("X:\Mathias\switch_data\foot_position\switch\pair6\ZmasterD1P6.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_m(2) = sum(cur_foot>threshold,"all")/numel(cur_foot);

cur_foot = load("X:\Mathias\switch_data\foot_position\switch\pair6\ZyokedD1P6.mat");
cur_foot = cur_foot.ZYoked(:,1);
foot_positions_y(2) = sum(cur_foot>threshold,"all")/numel(cur_foot);

% mouse 3
cur_foot = load("X:\Mathias\switch_data\foot_position\switch\pair7\ZmasterD1P7.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_m(3) = sum(cur_foot>threshold,"all")/numel(cur_foot);

cur_foot = load("X:\Mathias\switch_data\foot_position\switch\pair7\ZyokedD1P7.mat");
cur_foot = cur_foot.ZYoked(:,1);
foot_positions_y(3) = sum(cur_foot>threshold,"all")/numel(cur_foot);

% mouse 4
cur_foot = load("X:\Mathias\switch_data\foot_position\switch\pair8\ZmasterD1P8.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_m(4) = sum(cur_foot>threshold,"all")/numel(cur_foot);

cur_foot = load("X:\Mathias\switch_data\foot_position\switch\pair8\ZyokedD1P8.mat");
cur_foot = cur_foot.ZYoked(:,1);
foot_positions_y(4) = sum(cur_foot>threshold,"all")/numel(cur_foot);

% mouse 5
cur_foot = load("X:\Mathias\switch_data\foot_position\20220923_switch.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_m(5) = sum(cur_foot>threshold,"all")/numel(cur_foot);

cur_foot = load("X:\Mathias\switch_data\foot_position\20220927_switch_c.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_y(5) = sum(cur_foot>threshold,"all")/numel(cur_foot);

% mouse 6
cur_foot = load("X:\Mathias\switch_data\foot_position\20221011_switch.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_m(6) = sum(cur_foot>threshold,"all")/numel(cur_foot);

cur_foot = load("X:\Mathias\switch_data\foot_position\20221006_switch_c.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_y(6) = sum(cur_foot>threshold,"all")/numel(cur_foot);

% mouse 7
cur_foot = load("X:\Mathias\switch_data\foot_position\20221116_switch.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_m(7) = sum(cur_foot>threshold,"all")/numel(cur_foot);

cur_foot = load("X:\Mathias\switch_data\foot_position\20221018_switch_c.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_y(7) = sum(cur_foot>threshold,"all")/numel(cur_foot);

% mouse 8
cur_foot = load("X:\Mathias\switch_data\foot_position\20221124_switch.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_m(8) = sum(cur_foot>threshold,"all")/numel(cur_foot);

cur_foot = load("X:\Mathias\switch_data\foot_position\20220906_switch_c.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_y(8) = sum(cur_foot>threshold,"all")/numel(cur_foot);

% mouse 9
cur_foot = load("X:\Mathias\switch_data\foot_position\20221128_switch.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_m(9) = sum(cur_foot>threshold,"all")/numel(cur_foot);

cur_foot = load("X:\Mathias\switch_data\foot_position\20220413_switch_c.mat");
cur_foot = cur_foot.ZMaster(:,1);
foot_positions_y(8) = sum(cur_foot>threshold,"all")/numel(cur_foot);


%% calculations
results_m = load("X:\Mathias\switch_data\correlations\results_m.mat"); results_m = results_m.normalized_averages;
differences_m = (results_m(:,3)-results_m(:,1))./results_m(:,1);

figure
scatter(differences_m(2:9), foot_positions_m(2:9), 'filled')
hold on

% fit line to the data
p = polyfit(differences_m(2:9), foot_positions_m(2:9),1);
fit_line = polyval(p,differences_m(2:9));

% r-value
ss_tot = sum((foot_positions_m(2:9) - mean(foot_positions_m(2:9))).^2);
ss_res = sum((foot_positions_m(2:9) - fit_line).^2);
r_squared = 1 - ss_res/ss_tot;

plot(differences_m(2:9), fit_line)

