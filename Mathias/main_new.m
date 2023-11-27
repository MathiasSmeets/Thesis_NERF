clear; clc; close all; % make sure to reference paper lab if i end up using this code


%path_to_code = "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\";
path_to_code = "\\nerffs13\takeokalabwip2020\Mathias\data\";

stimulus_data_m = load(path_to_code + "data_after_stimulus_m.mat");
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_y = load(path_to_code + "data_after_stimulus_y.mat");
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

interval_size = 70;
%% filter operation for 0.1ms bins
% for i = 1:1454
%     stimulus_data_m{1,i}(:,101:106) = 0;
% end
%% create 1ms time bins also
% EXTRA_stimulus_data_m = cell(stimulus_data_m);
% EXTRA_stimulus_data_m{1,i}(:,71:700) = [];
% 
% for i = 1:1454
%     for k = 1:15
%         for j = 1:10:700
%             EXTRA_stimulus_data_m{1,i}(k,j) = sum(stimulus_data_m{1,i}(k,j:j+9));
%         end
%     end
%     EXTRA_stimulus_data_m{1,i}(:,71:end) = [];
% end
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

%% create neuron_spikes_m with 1ms timebin
% extra_neuron_spikes_m = zeros(size(neuron_spikes_m,1),70);
% for i = 1:size(neuron_spikes_m,1)
%     for j = 1:69
%         extra_neuron_spikes_m(i,j) = sum(neuron_spikes_m(i,j*10:j*10+9));
%     end
% end
% 
% mean_before_m = mean(extra_neuron_spikes_m(:,1:10), 2);
% std_before_m = std(extra_neuron_spikes_m(:,11:70), [], 2);
%neuron_spikes_m_zscore = (extra_neuron_spikes_m - mean_before_m) ./ std_before_m;


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


%% detect inhibited neurons using cusum

%decreased_FR_idx_m = find(norm_mean_after_m < norm_mean_before_m - 5*norm_std_before_m);
%decreased_FR_idx_y = find(norm_mean_after_y < norm_mean_before_y - 5*norm_std_before_y);

% inhibitory_m = [];
% ilower_m = cell(size(neuron_spikes_m,1),1);
% for i = 1:size(neuron_spikes_m,1)
%     if std(neuron_spikes_m(i,1:10)) == 0
%         [~, ilower_m{i,1}] = cusum_(neuron_spikes_m(i,:),5,1,mean(neuron_spikes_m(i,1:10)), 0.1);
%     else
%         [~, ilower_m{i,1}] = cusum_(neuron_spikes_m(i,:),5,1,mean(neuron_spikes_m(i,1:10)), std(neuron_spikes_m(i,1:10)));
%     end
%     if ~isempty(ilower_m{i,1})
%         inhibitory_m = [inhibitory_m, i];
%     end
% end
% 
% inhibitory_y = [];
% ilower_y = cell(size(neuron_spikes_y,1),1);
% for i = 1:size(neuron_spikes_y,1)
%     if std(neuron_spikes_y(i,1:10)) == 0
%         [~, ilower_y{i,1}] = cusum_(neuron_spikes_y(i,:),5,1,mean(neuron_spikes_y(i,1:10)), 0.1);
%     else
%         [~, ilower_y{i,1}] = cusum_(neuron_spikes_y(i,:),5,1,mean(neuron_spikes_y(i,1:10)), std(neuron_spikes_y(i,1:10)));
%     end
%     if ~isempty(ilower_y{i,1})
%         inhibitory_y = [inhibitory_y, i];
%     end
% end

%% paper lab: page 8 start, page 8 annex paper lab and page 17 extra figure 
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

% compute zscores of psth
std_before_m(std_before_m == 0) = 0.1;
std_before_y(std_before_y == 0) = 0.1;
neuron_spikes_m_zscore = (neuron_spikes_m - mean_before_m) ./ std_before_m;
neuron_spikes_y_zscore = (neuron_spikes_y - mean_before_y) ./ std_before_y;

peaks_locs_m = cell(size(neuron_spikes_m,1),1);
peaks_locs_y = cell(size(neuron_spikes_y,1),1);
inhibs_locs_m = cell(size(neuron_spikes_m,1),1);
inhibs_locs_y = cell(size(neuron_spikes_y,1),1);
for i = 1:size(neuron_spikes_m,1)
    event_amount = size(stimulus_data_m,2);
    %[peaks_locs_m{i,1}, peaks_locs_m{i,2}] = findpeaks(neuron_spikes_m(i,:), 'MinPeakProminence', 5, 'MinPeakHeight', mean_before_m(i,1) + 5*std_before_m(i,1));
    [peaks_locs_m{i,1}, peaks_locs_m{i,2}] = findpeaks(neuron_spikes_m_zscore(i,:), 'MinPeakHeight', 5);

    % detect inhibited neurons
    %[inhibs_locs_m{i,1}, inhibs_locs_m{i,2}] = findpeaks(-neuron_spikes_m(i,:), 'MinPeakProminence', 5, 'MinPeakHeight', mean_before_m(i,1) - 5*std_before_m(i,1));
    %if ~isempty(inhibs_locs_m{i,1})
    %    inhibited_neurons_m = [inhibited_neurons_m; i];
    %end
end
for i = 1:size(neuron_spikes_y,1)
    %[peaks_locs_y{i,1}, peaks_locs_y{i,2}] = findpeaks(neuron_spikes_y(i,:), 'MinPeakProminence', 5, 'MinPeakHeight', mean_before_y(i,1) + 5*std_before_y(i,1));
    [peaks_locs_y{i,1}, peaks_locs_y{i,2}] = findpeaks(neuron_spikes_y_zscore(i,:), 'MinPeakHeight', 5);

    % detect inhibted neurons
    %[inhibs_locs_y{i,1}, inhibs_locs_y{i,2}] = findpeaks(-neuron_spikes_y(i,:), 'MinPeakProminence', 5, 'MinPeakHeight', mean_before_y(i,1) - 5*std_before_y(i,1));
    %if ~isempty(inhibs_locs_y{i,1})
    %    inhibited_neurons_y = [inhibited_neurons_y; i];
    %end
end

%% detect inhibited neurons using cusum

%decreased_FR_idx_m = find(norm_mean_after_m < norm_mean_before_m - 5*norm_std_before_m);
%decreased_FR_idx_y = find(norm_mean_after_y < norm_mean_before_y - 5*norm_std_before_y);

inhibited_m = [];
ilower_m = cell(size(neuron_spikes_m_zscore,1),1);
for i = 1:size(neuron_spikes_m_zscore,1)
    [~, ilower_m{i,1}] = cusum_(neuron_spikes_m_zscore(i,:),5,1,0,2.5);

    if ~isempty(ilower_m{i,1})
        inhibited_m = [inhibited_m, i];
    end
end


%% secondary neurons:
% loop over each mouse -> neuron -> peaks -> stimulus
% here we loop over all total stimuli for this neuron -> compare with peaks of individual stimulus
% if 5% of stimuli has peak within interval of +-1ms --> keep this peak
% get highest peak within 5ms interval of this peak to get the std of the latency

secondary_neurons_m = [];
others_std_larger = [];
test_secondary = [];
neuron_counter_m = 0;
% loop over all mice
for i = 1:size(stimulus_data_m,1)
    % if this mouse is not empty
    if ~isempty(stimulus_data_m{i,1})
        % loop over all its neurons
        for j = 1:size(stimulus_data_m{i,1},1)
            disp(neuron_counter_m+j)
            count_inside_interval_m = 0;
            neuron_spikes_to_keep = [];
            % loop over all peaks (peaks are z_score > 5 from previous section)
            peaks_to_keep_m = []; 
            for k = 1:size(peaks_locs_m{neuron_counter_m+j,2},2)
                if ~isempty(peaks_locs_m{neuron_counter_m+j,2})
                    peak_counter = 0;
                    % get interval around peak: [peak-1ms : peak+1ms]
                    peak_interval_m = peaks_locs_m{neuron_counter_m+j,2}(1,k) - 1 : peaks_locs_m{neuron_counter_m+j,2}(1,k) + 1;
                    % loop over all stimuli for this neuron
                    stimulus_counter_m = 0;
                    for ii = 1:size(stimulus_data_m,2)
                        if ~isempty(stimulus_data_m{i,ii})
                            stimulus_counter_m = stimulus_counter_m + 1;
                            
                            cur_peaks = find(stimulus_data_m{i,ii}(j,peak_interval_m));
                            % check if found peak location is within interval
                            if ismember(2, cur_peaks)
                                peak_counter = peak_counter + 1;
                            end
                        end
                    end
                    if peak_counter >= 0.05 * stimulus_counter_m
                        peaks_to_keep_m = [peaks_to_keep_m, peaks_locs_m{neuron_counter_m+j,2}(1,k)];
                    end

                end
            end
            % now we have all peaks that occur frequently (at least 5% over all stimuli)
            % calculate psd for this
            % --> do this with highest peak in 5ms interval
            max_peak_distance = 5;
            
            if j+neuron_counter_m == 131
                disp("131")
            end

            for jj = 1:size(peaks_to_keep_m,2)
                peak_area_m = [];
                max_peaks = [];
                if jj ~= 0
                    if peaks_to_keep_m(jj) + max_peak_distance > interval_size % edge case
                        peak_area_m = neuron_spikes_m(neuron_counter_m+j,peaks_to_keep_m(jj) - max_peak_distance:end);
                    elseif  peaks_to_keep_m(jj) - max_peak_distance < 1 % edge case
                        peak_area_m = neuron_spikes_m(neuron_counter_m+j,1:peaks_to_keep_m(jj) + max_peak_distance);
                    else % normal case
                        peak_area_m = neuron_spikes_m(neuron_counter_m+j,peaks_to_keep_m(jj)-max_peak_distance:peaks_to_keep_m(jj)+max_peak_distance);
                    end
%                     for kk = 1:size(stimulus_data_m,2)
%                         if ~isempty(stimulus_data_m{i,kk})
%                             if peaks_to_keep_m(jj) + max_peak_distance > interval_size % edge case
%                                 [cur_peak_m, cur_loc_peak_m] = find(stimulus_data_m{i,kk}(j,peaks_to_keep_m(jj) - max_peak_distance:end));
%                                 %[cur_peak_m, cur_loc_peak_m] = findpeaks(cur_psth_zscore(1,peaks_to_keep_m(jj)-5:end),'NPeaks',1,'MinPeakHeight',5);
%                                 closest_peaks_m = [closest_peaks_m, cur_loc_peak_m+peaks_to_keep_m(jj) - max_peak_distance];
%                             elseif peaks_to_keep_m(jj) - max_peak_distance < 1 % edge case
%                                 [cur_peak_m, cur_loc_peak_m] = find(stimulus_data_m{i,kk}(j,1:peaks_to_keep_m(jj) + max_peak_distance));
%                                 %[cur_peak_m, cur_loc_peak_m] = findpeaks(cur_psth_zscore(1,1:peaks_to_keep_m(jj)+5),'NPeaks',1,'MinPeakHeight',5);
%                                 closest_peaks_m = [closest_peaks_m, cur_loc_peak_m];
%                             else % normal case
%                                 cur_peak = find(stimulus_data_m{i,kk}(j,peaks_to_keep_m(jj) - max_peak_distance:peaks_to_keep_m(jj) + max_peak_distance));
%                                 %[cur_peak_m, cur_loc_peak_m] = findpeaks(cur_psth_zscore(1,peaks_to_keep_m(jj)-5:peaks_to_keep_m(jj)+5),'NPeaks',1,'MinPeakHeight',5);
%                                 closest_peaks_m = [closest_peaks_m, cur_loc_peak_m+peaks_to_keep_m(jj) - max_peak_distance];
%                             end
%                             % assume this peak is not related if larger than 5ms
%                             % usually this data is either 0 or 1 so our peak detection does not have to be fancy
%                         end
%                    end
                end
                if neuron_counter_m + j == 967
                    disp("d")
                end

                closest_peaks_m = zeros(1,sum(peak_area_m));
                cur_index = 1;
                peak_area_m = peak_area_m - round(mean_before_m(neuron_counter_m + j,1));
                peak_area_m(peak_area_m<0) = 0;
                for kk = 1:size(peak_area_m,2)
                    if peaks_to_keep_m(jj) - max_peak_distance < 1
                        latency = kk-max_peak_distance-1 - (peaks_to_keep_m(jj) - max_peak_distance-1);
                    else
                        latency = kk-max_peak_distance-1;
                    end
                    closest_peaks_m(cur_index:cur_index+peak_area_m(kk)-1) = latency;
                    cur_index = cur_index + peak_area_m(kk);
                end
                
                if 1 > std(closest_peaks_m) && abs(mean(closest_peaks_m)) < 1
                %if 1 > std(abs(closest_peaks_m-peaks_to_keep_m(jj)))
                    %this neuron is a secondary neuron
                    secondary_neurons_m = [secondary_neurons_m; {neuron_counter_m+j, peaks_to_keep_m(jj)}];
                else
                    others_std_larger = [others_std_larger; {neuron_counter_m+j, peaks_to_keep_m(jj)}];
                end
            end
        end
        neuron_counter_m = neuron_counter_m + size(stimulus_data_m{i,1},1);
    end
end


% REPEAT FOR Y

%% 

output_m = table({secondary_neurons_m}, {others_std_larger}, {inhibited_m}, 'VariableNames',{'Secondary Neurons', 'Neurons with High STD Peaks', 'Inhibited Neurons'});

%% Wilcoxon signed rank test to check if there is a change in activity after stimulus
% assumes independent samples
% does NOT assume normal distribution

greater_mean = [];
lower_mean = [];
for i = 1:size(neuron_spikes_m,1)
    [p_greater, h_greater] = ranksum(neuron_spikes_m(i,11:end), neuron_spikes_m(i,1:10), "Tail", "right");
    [p_lower, h_lower] = ranksum( neuron_spikes_m(i,11:end), neuron_spikes_m(i,1:10), "Tail", "left");
    if h_greater == 1
        greater_mean = [greater_mean, i];
    end
    if h_lower == 1
        lower_mean = [lower_mean, i];
    end
end

% simes' correction?
