clear;close all;clc;

% missing:
% learner: switch 1
% control: switch 1 & 9
% np2 no switch so not included here

learner_indices = 2:9;
control_indices = 3:8;

trainer_m = zeros(1,9);
trainer_y = zeros(1,9);
switch_m = zeros(1,9);
switch_c = zeros(1,9);

threshold = 3;

%% learner
%% 20220725
trainer_1m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair9\ZmasterD1P9.mat");
trainer_1m = trainer_1m.ZMaster(:,1);
trainer_m(1) = sum(trainer_1m>threshold,"all")/numel(trainer_1m);

%switch_1m = load("")

%% 20220819
trainer_2m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair1\ZmasterD1P1.mat");
trainer_2m = trainer_2m.ZMaster(:,1);
trainer_m(2) = sum(trainer_2m>threshold,"all")/numel(trainer_2m);

switch_2m = load("X:\Mathias\switch_data\foot_positions\switch\pair6\ZmasterD1P6.mat");
switch_2m = switch_2m.ZMaster(:,1);
switch_m(2) = sum(switch_2m>threshold,"all")/numel(switch_2m);


%% 20220831
trainer_3m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair2\ZmasterD1P2.mat");
trainer_3m = trainer_3m.ZMaster(:,1);
trainer_m(3) = sum(trainer_3m>threshold,"all")/numel(trainer_3m);

switch_3m = load("X:\Mathias\switch_data\foot_positions\switch\pair7\ZmasterD1P7.mat");
switch_3m = switch_3m.ZMaster(:,1);
switch_m(3) = sum(switch_3m>threshold,"all")/numel(switch_3m);

%% 20220921

trainer_4m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair3\ZmasterD1P3.mat");
trainer_4m = trainer_4m.ZMaster(:,1);
trainer_m(4) = sum(trainer_4m>threshold,"all")/numel(trainer_4m);

switch_4m = load("X:\Mathias\switch_data\foot_positions\switch\pair8\ZmasterD1P8.mat");
switch_4m = switch_4m.ZMaster(:,1);
switch_m(4) = sum(switch_4m>threshold,"all")/numel(switch_4m);

%% 20220923

trainer_5m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair4\ZmasterD1P4.mat");
trainer_5m = trainer_5m.ZMaster(:,1);
trainer_m(5) = sum(trainer_5m>threshold,"all")/numel(trainer_5m);

switch_5m = load("X:\Mathias\switch_data\foot_positions\switch\20220923_switch.mat");
switch_5m = switch_5m.ZMaster(:,1);
switch_m(5) = sum(switch_5m>threshold,"all")/numel(switch_5m);


%% 20221011

trainer_6m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair5\ZmasterD1P5.mat");
trainer_6m = trainer_6m.ZMaster(:,1);
trainer_m(6) = sum(trainer_6m>threshold,"all")/numel(trainer_6m);

switch_6m = load("X:\Mathias\switch_data\foot_positions\switch\20221011_switch.mat");
switch_6m = switch_6m.ZMaster(:,1);
switch_m(6) = sum(switch_6m>threshold,"all")/numel(switch_6m);

%% 20221116

trainer_7m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair6\ZmasterD1P6.mat");
trainer_7m = trainer_7m.ZMaster(:,1);
trainer_m(7) = sum(trainer_7m>threshold,"all")/numel(trainer_7m);

switch_7m = load("X:\Mathias\switch_data\foot_positions\switch\20221116_switch.mat");
switch_7m = switch_7m.ZMaster(:,1);
switch_m(7) = sum(switch_7m>threshold,"all")/numel(switch_7m);

%% 20221124

trainer_8m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair7\ZmasterD1P7.mat");
trainer_8m = trainer_8m.ZMaster(:,1);
trainer_m(8) = sum(trainer_8m>threshold,"all")/numel(trainer_8m);

switch_8m = load("X:\Mathias\switch_data\foot_positions\switch\20221124_switch.mat");
switch_8m = switch_8m.ZMaster(:,1);
switch_m(8) = sum(switch_8m>threshold,"all")/numel(switch_8m);

%% 20221128

trainer_9m = load("X:\Mathias\switch_data\foot_positions\training\training_master\pair8\ZmasterD1P8.mat");
trainer_9m = trainer_9m.ZMaster(:,1);
trainer_m(9) = sum(trainer_9m>threshold,"all")/numel(trainer_9m);

switch_9m = load("X:\Mathias\switch_data\foot_positions\switch\20221128_switch.mat");
switch_9m = switch_9m.ZMaster(:,1);
switch_m(9) = sum(switch_9m>threshold,"all")/numel(switch_9m);

%% control
%% 20220705

%trainer_1c =

%% 20220726 (pair 8)

trainer_2c = load("X:\Mathias\switch_data\foot_positions\training\training_control\pair8\ZmasterD1P8.mat");
trainer_2c = trainer_2c.ZMaster(:,1);
trainer_c(2) = sum(trainer_2c>threshold,"all")/numel(trainer_2c);


%% 20220905 (pair 1)

trainer_3c = load("X:\Mathias\switch_data\foot_positions\training\training_control\pair1\ZmasterD1P1.mat");
trainer_3c = trainer_3c.ZMaster(:,1);
trainer_c(3) = sum(trainer_3c>threshold,"all")/numel(trainer_3c);

switch_3c = load("X:\Mathias\switch_data\foot_positions\switch\20220905_switch_c.mat");
switch_3c = switch_3c.ZMaster(:,1);
switch_c(3) = sum(switch_3c>threshold,"all")/numel(switch_3c);

%% 20220922 (pair 2)

trainer_4c = load("X:\Mathias\switch_data\foot_positions\training\training_control\pair2\ZmasterD1P2.mat");
trainer_4c = trainer_4c.ZMaster(:,1);
trainer_c(4) = sum(trainer_4c>threshold,"all")/numel(trainer_4c);

switch_4c = load("X:\Mathias\switch_data\foot_positions\switch\20220922_switch_c.mat");
switch_4c = switch_4c.ZMaster(:,1);
switch_c(4) = sum(switch_4c>threshold,"all")/numel(switch_4c);

%% 20220927 (pair 3)

trainer_5c = load("X:\Mathias\switch_data\foot_positions\training\training_control\pair3\ZmasterD1P3.mat");
trainer_5c = trainer_5c.ZMaster(:,1);
trainer_c(5) = sum(trainer_5c>threshold,"all")/numel(trainer_5c);

switch_5c = load("X:\Mathias\switch_data\foot_positions\switch\20220927_switch_c.mat");
switch_5c = switch_5c.ZMaster(:,1);
switch_c(5) = sum(switch_5c>threshold,"all")/numel(switch_5c);

%% 20221006 (pair 4)

trainer_6c = load("X:\Mathias\switch_data\foot_positions\training\training_control\pair4\ZmasterD1P4.mat");
trainer_6c = trainer_6c.ZMaster(:,1);
trainer_c(6) = sum(trainer_6c>threshold,"all")/numel(trainer_6c);

switch_6c = load("X:\Mathias\switch_data\foot_positions\switch\20221006_switch_c.mat");
switch_6c = switch_6c.ZMaster(:,1);
switch_c(6) = sum(switch_6c>threshold,"all")/numel(switch_6c);

%% 20221018 (pair 5)

trainer_7c = load("X:\Mathias\switch_data\foot_positions\training\training_control\pair5\ZmasterD1P5.mat");
trainer_7c = trainer_7c.ZMaster(:,1);
trainer_c(7) = sum(trainer_7c>threshold,"all")/numel(trainer_7c);

switch_7c = load("X:\Mathias\switch_data\foot_positions\switch\20221018_switch_c.mat");
switch_7c = switch_7c.ZMaster(:,1);
switch_c(7) = sum(switch_7c>threshold,"all")/numel(switch_7c);

%% 20220906 (pair 6)

trainer_8c = load("X:\Mathias\switch_data\foot_positions\training\training_control\pair6\ZmasterD1P6.mat");
trainer_8c = trainer_8c.ZMaster(:,1);
trainer_c(8) = sum(trainer_8c>threshold,"all")/numel(trainer_8c);

switch_8c = load("X:\Mathias\switch_data\foot_positions\switch\20220906_switch_c.mat");
switch_8c = switch_8c.ZMaster(:,1);
switch_c(8) = sum(switch_8c>threshold,"all")/numel(switch_8c);

%% 20220413 (pair 7)

trainer_9c = load("X:\Mathias\switch_data\foot_positions\training\training_control\pair7\ZmasterD1P7.mat");
trainer_9c = trainer_9c.ZMaster(:,1);
trainer_c(9) = sum(trainer_9c>threshold,"all")/numel(trainer_9c);

%% calculations
results_m = load("X:\Mathias\switch_data\correlations\results_m.mat"); results_m = results_m.normalized_averages;
differences_m = (results_m(:,3)-results_m(:,1));%./results_m(:,1);

figure
scatter(differences_m(learner_indices), switch_m(learner_indices)-trainer_m(learner_indices), 'filled')
hold on

% fit line to the data
p = polyfit(differences_m(learner_indices), switch_m(learner_indices)-trainer_m(learner_indices),1);
fit_line = polyval(p,differences_m(learner_indices));

% r-value
ss_tot = sum((switch_m(learner_indices)-trainer_m(learner_indices) - mean(switch_m(learner_indices)-trainer_m(learner_indices))).^2);
ss_res = sum((switch_m(learner_indices)-trainer_m(learner_indices) - fit_line').^2);
r_squared = 1 - ss_res/ss_tot;

plot(differences_m(learner_indices), fit_line)

correlation_coefficient_m = corrcoef(transpose(switch_m(learner_indices)-trainer_m(learner_indices)),differences_m(learner_indices));
disp("Correlation Learner: "+correlation_coefficient_m(1,2))
% test for significance
n = numel(learner_indices);
t_statistic = correlation_coefficient_m(1,2) * sqrt((n - 2) / (1 - correlation_coefficient_m(1,2)^2));
alpha = 0.05;
critical_t = tinv(1 - alpha/2, n - 2);
if abs(t_statistic) > critical_t
    disp('The correlation coefficient is statistically significant.');
else
    disp('The correlation coefficient is not statistically significant.');
end

% y
results_y = load("X:\Mathias\switch_data\correlations\results_y.mat"); results_y = results_y.normalized_averages;
differences_y = (results_y(:,3)-results_y(:,1));%./results_y(:,1);

figure
scatter(differences_y(control_indices), switch_c(control_indices)-trainer_c(control_indices), 'filled')
hold on

% fit line to the data
p = polyfit(differences_y(control_indices), switch_c(control_indices)-trainer_c(control_indices),1);
fit_line = polyval(p,differences_y(control_indices));

% r-value
ss_tot = sum((switch_c(control_indices)-trainer_c(control_indices) - mean(switch_c(control_indices)-trainer_c(control_indices))).^2);
ss_res = sum((switch_c(control_indices)-trainer_c(control_indices) - fit_line').^2);
r_squared = 1 - ss_res/ss_tot;

plot(differences_y(control_indices), fit_line)

correlation_coefficient_c = corrcoef(transpose(switch_m(control_indices)-trainer_m(control_indices)),differences_y(control_indices));
disp("Correlation Control: "+correlation_coefficient_c(1,2))

% test for significance
n = numel(control_indices);
t_statistic = correlation_coefficient_m(1,2) * sqrt((n - 2) / (1 - correlation_coefficient_m(1,2)^2));
p_value = 2 * (1 - tcdf(abs(t_statistic), n - 2));

disp("P-value for correlation: " + p_value)
