clear; clc; close all; % make sure to reference paper lab if i end up using this code

stimulus_data_m = load("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\data_after_stimulus_m.mat");
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_y = load("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\data_after_stimulus_y.mat");
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

interval_size = 70; 
%% get total amount of spikes + amount of spikes of each neuron
total_m = zeros(1,interval_size);
total_y = zeros(1,interval_size);

total_nb_neurons_m = 0;
total_nb_neurons_y = 0;
for i = 1:size(stimulus_data_m,1)
    total_nb_neurons_m = total_nb_neurons_m + size(stimulus_data_m{i,1},1);
end
for i = 1:size(stimulus_data_y,1)
    total_nb_neurons_y = total_nb_neurons_y + size(stimulus_data_y{i,1},1);
end
neuron_spikes_m = zeros(total_nb_neurons_m, interval_size);
neuron_spikes_y = zeros(total_nb_neurons_y, interval_size);

cur_neuron_nb_m = 1;
cur_neuron_nb_y = 1;
for i = 1:size(stimulus_data_m,1)
    for j = 1:size(stimulus_data_m,2)

        if ~isempty(stimulus_data_m{i,j})
            if size(stimulus_data_m{i,j},2) < size(total_m,2)
                total_m(1:size(stimulus_data_m{i,j},2)) = total_m(1:size(stimulus_data_m{i,j},2)) + sum(stimulus_data_m{i,j});
                neuron_spikes_m(cur_neuron_nb_m:cur_neuron_nb_m + size(stimulus_data_m{i,j},1) - 1, 1:size(stimulus_data_m{i,j},2)) = neuron_spikes_m(cur_neuron_nb_m:cur_neuron_nb_m + size(stimulus_data_m{i,j}, 1) - 1 ,1:size(stimulus_data_m{i,j},2)) + stimulus_data_m{i,j};
            else
                total_m = total_m + sum(stimulus_data_m{i,j});
                neuron_spikes_m(cur_neuron_nb_m:cur_neuron_nb_m + size(stimulus_data_m{i,j},1) - 1, :) = neuron_spikes_m(cur_neuron_nb_m:cur_neuron_nb_m + size(stimulus_data_m{i,j}, 1) - 1, :) + stimulus_data_m{i,j};
            end
        end
        if ~isempty(stimulus_data_y{i,j})
            if size(stimulus_data_y{i,j},2) < size(total_y,2)
                total_y(1:size(stimulus_data_y{i,j},2)) = total_y(1:size(stimulus_data_y{i,j},2)) + sum(stimulus_data_y{i,j});
                neuron_spikes_y(cur_neuron_nb_y:cur_neuron_nb_y + size(stimulus_data_y{i,j},1) - 1, 1:size(stimulus_data_y{i,j},2)) = neuron_spikes_y(cur_neuron_nb_y:cur_neuron_nb_y + size(stimulus_data_y{i,j}, 1) - 1 ,1:size(stimulus_data_y{i,j},2)) + stimulus_data_y{i,j};
            else
                total_y = total_y + sum(stimulus_data_y{i,j});
                neuron_spikes_y(cur_neuron_nb_y:cur_neuron_nb_y + size(stimulus_data_y{i,j}, 1) - 1, :) = neuron_spikes_y(cur_neuron_nb_y:cur_neuron_nb_y + size(stimulus_data_y{i,j}, 1) - 1, :) + stimulus_data_y{i,j};
            end
        end
    end
    cur_neuron_nb_m = cur_neuron_nb_m + size(stimulus_data_m{i,1},1);
    cur_neuron_nb_y = cur_neuron_nb_y + size(stimulus_data_y{i,1},1);
end

%% make figure for total spikes
figure
plot(total_m)
hold on
plot(total_y)
title("Total Amount of Spikes")
xlabel("Time After Stimulation (ms)")
ylabel("Amount of Spikes (over all tests)")
legend("Learner", "Control")

% figure
% hold on
% for i = 1:total_nb_neurons_m
%     plot(neuron_spikes_m(i,:))
% end
% figure
% hold on
% for i = 1:total_nb_neurons_y
%     plot(neuron_spikes_y(i,:))
% end

%% perform PCA on neuron_spikes
% normalized_spikes_neurons_m = zscore(neuron_spikes_m);
% 
% coeff_m = pca(normalized_spikes_neurons_m);
% 
% covariance_m = cov(normalized_spikes_neurons_m);
% eigenvalues_m = eig(covariance_m);
% figure
% stem(eigenvalues_m);
% 
% P = coeff_m(:,1:12);
% reduced_spikes_m = P' * normalized_spikes_neurons_m';
% 
% [indices_m, c_m] = kmeans(reduced_spikes_m',5,'Replicates',10);

%% plot spikes on the different figures corresponding to different classes

% t = tiledlayout(5,1);
% ax1 = nexttile;
% ax2 = nexttile;
% ax3 = nexttile;
% ax4 = nexttile;
% ax5 = nexttile;
% 
% hold([ax1, ax2, ax3, ax4, ax5] ,'on')
% for i = 1:size(indices_m,1)
%     if indices_m(i,1) == 1
%         plot(ax1, neuron_spikes_m(i,:))
%     elseif indices_m(i,1) == 2
%         plot(ax2, neuron_spikes_m(i,:))
%     elseif indices_m(i,1) == 3
%         plot(ax3, neuron_spikes_m(i,:))
%     elseif indices_m(i,1) == 4
%         plot(ax4, neuron_spikes_m(i,:))
%     elseif indices_m(i,1) == 5
%         plot(ax5, neuron_spikes_m(i,:))
%     end
% end

%% paper lab: page 8 start, page 8 annex paper lab and page 17 extra figure 

% compute z-scores of the PeriStimulus Time Histograms (PSTH)

% identify presence of increase in firing rate (>5 SD compared to 10 ms window before stimulation)

% criteria for second-order units: consistency of response to conditioning
% cues is over 5% of the entire trial and SD of response shorter than 1ms




% PSTH already computed in previous section (neuron_spikes_m/y)
normalised_m = zscore(neuron_spikes_m, [], 2);
normalised_y = zscore(neuron_spikes_y, [], 2);

mean_before_m = mean(normalised_m(:,1:10), 2);
mean_before_y = mean(normalised_y(:,1:10), 2);
mean_after_m = mean(normalised_m(:,11:70), 2);
mean_after_y = mean(normalised_y(:,11:70), 2);

std_before_m = std(normalised_m(:,1:10), [], 2);
std_before_y = std(normalised_y(:,1:10), [], 2);
std_after_m = std(normalised_m(:,11:70), [] , 2);
std_after_y = std(normalised_y(:,11:70), [] , 2);

increased_FR_idx_m = find(mean_after_m > mean_before_m + 5*std_before_m);
increased_FR_idx_y = find(mean_after_y > mean_before_y + 5*std_before_y);

% problem at the moment: if no spikes in before: will detect after as
% significant even when we only have like 4 spikes

%% find peaks
[peaks_m_max, loc_peaks_m_max] = max(neuron_spikes_m, [], 2);
[peaks_y_max, loc_peaks_y_max] = max(neuron_spikes_y, [], 2);
peak_interval_m = loc_peaks_m_max - 1 : loc_peaks_m_max + 1;
peak_interval_y = loc_peaks_y_max - 1 : loc_peaks_y_max + 1;


% peaks_locs_m = cell(size(neuron_spikes_m,2),1);
% peaks_locs_y = cell(size(neuron_spikes_y,2),1);
% for i = 1:size(neuron_spikes_m,1)
%     [peaks_locs_m{i,1}, peaks_locs_m{i,2}] = findpeaks(neuron_spikes_m(i,:),'MinPeakProminence',20, 'MinPeakHeight', 20);                                       % TODO
% end
% for i = 1:size(neuron_spikes_y,1)
%     
%     [peaks_locs_y{i,1}, peaks_locs_y{i,2}] = findpeaks(neuron_spikes_y(i,:),'MinPeakProminence',20, 'MinPeakHeight', 20);                                       % TODO
% end

%% secondary neurons:
%  get peaks + check consistency of response to conditioning cues (>5%) + calculate response latency (<1ms)

inside_interval_indices_m = [];
response_latency_sd_low_m = [];
neuron_counter_m = 0;
% loop over all mice
for i = 1:size(stimulus_data_m,1)
    % if this mouse is not empty
    if ~isempty(stimulus_data_m{i,1})
        % loop over all its neurons
        for j = 1:size(stimulus_data_m{i,1},1)
            count_inside_interval_m = 0;
            total_counter_m = 0;
            peak_locations = [];
            % loop over all stimuli
            for k = 1:size(stimulus_data_m,2)
                if ~isempty(stimulus_data_m{i,k})
                    total_counter_m = total_counter_m + 1;
                    [peak_m, loc_peak_m] = findpeaks(stimulus_data_m{i,k}(j,:));

                    % check if a peak is within wanted interval
                    inside_interval = false;
                    for ii = length(loc_peak_m)
                        if ii ~= 0
                            if ismember(loc_peak_m(ii), peak_interval_m) 
                                inside_interval = true;
                            end
                        end
                    end
                    if inside_interval
                        count_inside_interval_m = count_inside_interval_m + 1;
                    end

                    % find location closest to overall peak?
                    [~, closest_peak] = min(abs(loc_peak_m-loc_peaks_m_max(neuron_counter_m + j,1)));

                    % append to all peak locations
                    peak_locations = [peak_locations, closest_peak];
                end
            end
            if count_inside_interval_m / total_counter_m > 0.05
                inside_interval_indices_m = [inside_interval_indices_m; neuron_counter_m + j];
            end
            % check if SD of response latency > 1ms
            if std(peak_locations) < 1
                response_latency_sd_low_m = [response_latency_sd_low_m; neuron_counter_m + j];
            end
        end
        neuron_counter_m = neuron_counter_m + size(stimulus_data_m{i,1},1);
    end
end


% REPEAT FOR Y



% first condition filters out a lot                                                                                                                                                        %TODO


%% check for neurons that comply to all conditions


% inside_interval_indices_m & response_latency_sd_low_m & increased_FR_idx_m

neurons_meet_conditions_m = intersect(intersect(inside_interval_indices_m, response_latency_sd_low_m),increased_FR_idx_m);
neurons_meet_conditions_y = intersect(intersect(inside_interval_indices_y, response_latency_sd_low_y),increased_FR_idx_y);

