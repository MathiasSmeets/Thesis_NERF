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

norm_mean_before_m = mean(normalised_m(:,1:10), 2);
norm_mean_before_y = mean(normalised_y(:,1:10), 2);
norm_mean_after_m = mean(normalised_m(:,11:70), 2);
norm_mean_after_y = mean(normalised_y(:,11:70), 2);

norm_std_before_m = std(normalised_m(:,1:10), [], 2);
norm_std_before_y = std(normalised_y(:,1:10), [], 2);
norm_std_after_m = std(normalised_m(:,11:70), [] , 2);
norm_std_after_y = std(normalised_y(:,11:70), [] , 2);

increased_FR_idx_m = find(norm_mean_after_m > norm_mean_before_m + 5*norm_std_before_m);
increased_FR_idx_y = find(norm_mean_after_y > norm_mean_before_y + 5*norm_std_before_y);        % --> do not really need this?

% problem at the moment: if no spikes in before: will detect after as
% significant even when we only have like 4 spikes

%% detect inhibited neurons using cusum

decreased_FR_idx_m = find(norm_mean_after_m < norm_mean_before_m - 5*norm_std_before_m);
decreased_FR_idx_y = find(norm_mean_after_y < norm_mean_before_y - 5*norm_std_before_y);

inhibitory_m = [];
ilower_m = cell(size(neuron_spikes_m,1),1);
for i = 1:size(neuron_spikes_m,1)
    if std(neuron_spikes_m(i,1:10)) == 0
        [~, ilower_m{i,1}] = cusum_(neuron_spikes_m(i,:),5,1,mean(neuron_spikes_m(i,1:10)), 0.1);
    else
        [~, ilower_m{i,1}] = cusum_(neuron_spikes_m(i,:),5,1,mean(neuron_spikes_m(i,1:10)), std(neuron_spikes_m(i,1:10)));
    end
    if ~isempty(ilower_m{i,1})
        inhibitory_m = [inhibitory_m, i];
    end
end

inhibitory_y = [];
ilower_y = cell(size(neuron_spikes_y,1),1);
for i = 1:size(neuron_spikes_y,1)
    if std(neuron_spikes_y(i,1:10)) == 0
        [~, ilower_y{i,1}] = cusum_(neuron_spikes_y(i,:),5,1,mean(neuron_spikes_y(i,1:10)), 0.1);
    else
        [~, ilower_y{i,1}] = cusum_(neuron_spikes_y(i,:),5,1,mean(neuron_spikes_y(i,1:10)), std(neuron_spikes_y(i,1:10)));
    end
    if ~isempty(ilower_y{i,1})
        inhibitory_y = [inhibitory_y, i];
    end
end

%% find peaks
[peaks_m_max, loc_peaks_m_max] = max(neuron_spikes_m, [], 2);
[peaks_y_max, loc_peaks_y_max] = max(neuron_spikes_y, [], 2);
peak_interval_m = loc_peaks_m_max - 1 : loc_peaks_m_max + 1;
peak_interval_y = loc_peaks_y_max - 1 : loc_peaks_y_max + 1;

mean_before_m = mean(neuron_spikes_m(:,1:10), 2);
mean_before_y = mean(neuron_spikes_y(:,1:10), 2);
mean_after_m = mean(neuron_spikes_m(:,11:70), 2);
mean_after_y = mean(neuron_spikes_y(:,11:70), 2);

std_before_m = std(neuron_spikes_m(:,1:10), [], 2);
std_before_y = std(neuron_spikes_y(:,1:10), [], 2);
std_after_m = std(neuron_spikes_m(:,11:70), [] , 2);
std_after_y = std(neuron_spikes_y(:,11:70), [] , 2);

peaks_locs_m = cell(size(neuron_spikes_m,2),1);
peaks_locs_y = cell(size(neuron_spikes_y,2),1);
inhibs_locs_m = cell(size(neuron_spikes_m,2),1);
inhibs_locs_y = cell(size(neuron_spikes_y,2),1);
inhibited_neurons_m = [];
inhibited_neurons_y = [];
for i = 1:size(neuron_spikes_m,1)
    event_amount = size(stimulus_data_m,2);
    [peaks_locs_m{i,1}, peaks_locs_m{i,2}] = findpeaks(neuron_spikes_m(i,:), 'MinPeakProminence', 5, 'MinPeakHeight', mean_before_m(i,1) + 5*std_before_m(i,1));
    %[peaks_locs_m{i,1}, peaks_locs_m{i,2}] = findpeaks(neuron_spikes_m(i,:),'MinPeakProminence',20, 'MinPeakHeight', 20);

    % detect inhibited neurons
    [inhibs_locs_m{i,1}, inhibs_locs_m{i,2}] = findpeaks(-neuron_spikes_m(i,:), 'MinPeakProminence', 5, 'MinPeakHeight', mean_before_m(i,1) - 5*std_before_m(i,1));
    if ~isempty(inhibs_locs_m{i,1})
        inhibited_neurons_m = [inhibited_neurons_m; i];
    end
end
for i = 1:size(neuron_spikes_y,1)
    [peaks_locs_y{i,1}, peaks_locs_y{i,2}] = findpeaks(neuron_spikes_y(i,:), 'MinPeakProminence', 5, 'MinPeakHeight', mean_before_y(i,1) + 5*std_before_y(i,1));
    %[peaks_locs_y{i,1}, peaks_locs_y{i,2}] = findpeaks(neuron_spikes_y(i,:),'MinPeakProminence',20, 'MinPeakHeight', 20);

    % detect inhibted neurons
    [inhibs_locs_y{i,1}, inhibs_locs_y{i,2}] = findpeaks(-neuron_spikes_y(i,:), 'MinPeakProminence', 5, 'MinPeakHeight', mean_before_y(i,1) - 5*std_before_y(i,1));
    if ~isempty(inhibs_locs_y{i,1})
        inhibited_neurons_y = [inhibited_neurons_y; i];
    end
end

%% secondary neurons:
% loop over each mouse -> neuron -> peaks -> stimulus
% here we loop over all total stimuli for this neuron -> compare with peaks of individual stimulus
% if 5% of stimuli has peak within interval of +-1ms --> keep this peak
% get shortest peak to this peak for each stimulus to calculate std


secondary_neurons_m = [];
others_std_larger = [];
neuron_counter_m = 0;
% loop over all mice
for i = 1:size(stimulus_data_m,1)
    % if this mouse is not empty
    if ~isempty(stimulus_data_m{i,1})
        % loop over all its neurons
        for j = 1:size(stimulus_data_m{i,1},1)
            count_inside_interval_m = 0;
            neuron_spikes_to_keep = [];
            % loop over all peaks
            peaks_to_keep_m = []; 
            for k = 1:size(peaks_locs_m{neuron_counter_m+j,2},2)
                if ~isempty(peaks_locs_m{neuron_counter_m+j,2})
                    peak_counter = 0;
                    peak_interval_m = peaks_locs_m{neuron_counter_m+j,2}(1,k) - 1 : peaks_locs_m{neuron_counter_m+j,2}(1,k) + 1;
                    % loop over all stimuli
                    stimulus_counter_m = 0;
                    for ii = 1:size(stimulus_data_m,2)
                        if ~isempty(stimulus_data_m{i,ii})
                            stimulus_counter_m = stimulus_counter_m + 1;
                            cur_mean_before_m = mean(stimulus_data_m{i,ii}(j,1:10),2);
                            cur_std_before_m = std(stimulus_data_m{i,ii}(j,1:10), [], 2);
                            [cur_peak_m, cur_loc_peak_m] = findpeaks(stimulus_data_m{i,ii}(j,:));
                            members_idx_m = ismember(cur_loc_peak_m, peak_interval_m);
                            if sum(members_idx_m) > 0
                                peak_counter = peak_counter + 1;
                            end
                        end
                    end
                    if peak_counter >= 0.05 * stimulus_counter_m
                        peaks_to_keep_m = [peaks_to_keep_m, peaks_locs_m{neuron_counter_m+j,2}(1,k)];
                    end

                end
            end
            % now we have all peaks that occur frequently
            % calculate psd for this
            % --> do this with peak closest to total peak we kept

            for jj = 1:size(peaks_to_keep_m)
                closest_peaks_m = [];
                if jj ~= 0
                    %loop over stimuli
                    for kk = 1:size(stimulus_data_m,2)
                        if ~isempty(stimulus_data_m{i,kk})
                            if peaks_to_keep_m(jj) + 5 > 70
                                [cur_peak_m, cur_loc_peak_m] = findpeaks(stimulus_data_m{i,kk}(j,peaks_to_keep_m(jj)-5:end),'NPeaks',1);
                                closest_peaks_m = [closest_peaks_m, cur_loc_peak_m+peaks_to_keep_m(jj)-6];
                            elseif peaks_to_keep_m(jj)-5 < 1
                                [cur_peak_m, cur_loc_peak_m] = findpeaks(stimulus_data_m{i,kk}(j,1:peaks_to_keep_m(jj)+5),'NPeaks',1);
                                closest_peaks_m = [closest_peaks_m, cur_loc_peak_m];
                            else
                                [cur_peak_m, cur_loc_peak_m] = findpeaks(stimulus_data_m{i,kk}(j,peaks_to_keep_m(jj)-5:peaks_to_keep_m(jj)+5),'NPeaks',1);
                                closest_peaks_m = [closest_peaks_m, cur_loc_peak_m+peaks_to_keep_m(jj)-6];
                            end
                            %[~, closest_peak_idx] = min(abs(cur_loc_peak_m - peaks_to_keep_m(jj)));
                            %cur_closest_peak_m = cur_loc_peak_m(closest_peak_idx);
                            %if abs(peaks_to_keep_m(jj) - cur_closest_peak_m) < 5     % assume this peak is not related if larger than 5ms
                        
                            %end
                        end
                    end
                end
                if 1 > std(closest_peaks_m)
                %if 1 > std(abs(closest_peaks_m-peaks_to_keep_m(jj)))
                    %this neuron is a secondary neuron
                    secondary_neurons_m = [secondary_neurons_m; neuron_counter_m+j];
                else
                    others_std_larger = [others_std_larger; neuron_counter_m+j];
                end
            end
        end
        neuron_counter_m = neuron_counter_m + size(stimulus_data_m{i,1},1);
    end
end


% REPEAT FOR Y

%% 




