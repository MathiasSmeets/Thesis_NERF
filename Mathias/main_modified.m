clear; clc; close all; % make sure to reference paper lab if i end up using this code

%path_to_code = "\\nerffs13\takeokalabwip2020\Mathias\data\";
path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/10kfs/";
other_path = "/mnt/takeokalab/takeokalabwip2023/Mathias/data/";

stimulus_data_m = load(other_path + "data_after_stimulus_m.mat");
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;
stimulus_data_y = load(other_path + "data_after_stimulus_y.mat");
stimulus_data_y = stimulus_data_y.after_stimulus_data_y;

stimulus_data_high_fs_m = load(path_to_code + "after_stimulus_data_m_high_fs.mat");
stimulus_data_high_fs_m = stimulus_data_high_fs_m.after_stimulus_data_m;
stimulus_data_high_fs_y = load(path_to_code + "after_stimulus_data_y_high_fs.mat");
stimulus_data_high_fs_y = stimulus_data_high_fs_y.after_stimulus_data_y;

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


%% make neuron_spikes_m for 0.1ms
total_m_hfs = zeros(1,interval_size*10);
total_y_hfs = zeros(1,interval_size*10);

total_nb_neurons_m_hfs = 0;
total_nb_neurons_y_hfs = 0;
for i = 1:size(stimulus_data_high_fs_m,1)
    total_nb_neurons_m_hfs = total_nb_neurons_m_hfs + size(stimulus_data_high_fs_m{i,1},1);
end
for i = 1:size(stimulus_data_high_fs_y,1)
    total_nb_neurons_y_hfs = total_nb_neurons_y_hfs + size(stimulus_data_high_fs_y{i,1},1);
end
neuron_spikes_m_hfs = zeros(total_nb_neurons_m_hfs, interval_size*10);
neuron_spikes_y_hfs = zeros(total_nb_neurons_y_hfs, interval_size*10);

cur_neuron_nb_m_hfs = 1;
cur_neuron_nb_y_hfs = 1;
for i = 1:size(stimulus_data_high_fs_m,1)
    for j = 1:size(stimulus_data_high_fs_m,2)
        if ~isempty(stimulus_data_high_fs_m{i,j})
            if size(stimulus_data_high_fs_m{i,j},2) < size(total_m_hfs,2)
                total_m_hfs(1:size(stimulus_data_high_fs_m{i,j},2)) = total_m_hfs(1:size(stimulus_data_high_fs_m{i,j},2)) + sum(stimulus_data_high_fs_m{i,j});
                neuron_spikes_m_hfs(cur_neuron_nb_m_hfs:cur_neuron_nb_m_hfs + size(stimulus_data_high_fs_m{i,j},1) - 1, 1:size(stimulus_data_high_fs_m{i,j},2)) = neuron_spikes_m_hfs(cur_neuron_nb_m_hfs:cur_neuron_nb_m_hfs + size(stimulus_data_high_fs_m{i,j}, 1) - 1 ,1:size(stimulus_data_high_fs_m{i,j},2)) + stimulus_data_high_fs_m{i,j};
            else
                total_m_hfs = total_m_hfs + sum(stimulus_data_high_fs_m{i,j});
                neuron_spikes_m_hfs(cur_neuron_nb_m_hfs:cur_neuron_nb_m_hfs + size(stimulus_data_high_fs_m{i,j},1) - 1, :) = neuron_spikes_m_hfs(cur_neuron_nb_m_hfs:cur_neuron_nb_m_hfs + size(stimulus_data_high_fs_m{i,j}, 1) - 1, :) + stimulus_data_high_fs_m{i,j};
            end
        end
    end
    cur_neuron_nb_m_hfs = cur_neuron_nb_m_hfs + size(stimulus_data_high_fs_m{i,1},1);
end
for i = 1:size(stimulus_data_high_fs_y,1)
    for j = 1:size(stimulus_data_high_fs_y,2)
        if ~isempty(stimulus_data_high_fs_y{i,j})
            if size(stimulus_data_high_fs_y{i,j},2) < size(total_y_hfs,2)
                total_y_hfs(1:size(stimulus_data_high_fs_y{i,j},2)) = total_y_hfs(1:size(stimulus_data_high_fs_y{i,j},2)) + sum(stimulus_data_high_fs_y{i,j});
                neuron_spikes_y_hfs(cur_neuron_nb_y_hfs:cur_neuron_nb_y_hfs + size(stimulus_data_high_fs_y{i,j},1) - 1, 1:size(stimulus_data_high_fs_y{i,j},2)) = neuron_spikes_y_hfs(cur_neuron_nb_y_hfs:cur_neuron_nb_y_hfs + size(stimulus_data_high_fs_y{i,j}, 1) - 1 ,1:size(stimulus_data_high_fs_y{i,j},2)) + stimulus_data_high_fs_y{i,j};
            else
                total_y_hfs = total_y_hfs + sum(stimulus_data_high_fs_y{i,j});
                neuron_spikes_y_hfs(cur_neuron_nb_y_hfs:cur_neuron_nb_y_hfs + size(stimulus_data_high_fs_y{i,j}, 1) - 1, :) = neuron_spikes_y_hfs(cur_neuron_nb_y_hfs:cur_neuron_nb_y_hfs + size(stimulus_data_high_fs_y{i,j}, 1) - 1, :) + stimulus_data_high_fs_y{i,j};
            end
        end
    end
    cur_neuron_nb_y_hfs = cur_neuron_nb_y_hfs + size(stimulus_data_high_fs_y{i,1},1);
end

neuron_spikes_m(:,9:11) = 0;
neuron_spikes_y(:,9:11) = 0;
neuron_spikes_m_hfs(:,91:110) = 0;
neuron_spikes_y_hfs(:,91:110) = 0;

save("/scratch/mathiass-takeokalab/01/neuron_spikes_m", "neuron_spikes_m", "-v7.3");
save("/scratch/mathiass-takeokalab/01/neuron_spikes_y", "neuron_spikes_y", "-v7.3");
save("/scratch/mathiass-takeokalab/01/neuron_spikes_m_hfs", "neuron_spikes_m_hfs", "-v7.3");
save("/scratch/mathiass-takeokalab/01/neuron_spikes_y_hfs", "neuron_spikes_y_hfs", "-v7.3");

%% make figure for total spikes
figure
plot(total_m_hfs)
hold on
plot(total_y_hfs)
title("Total Amount of Spikes")
xlabel("Time After Stimulation (ms)")
ylabel("Amount of Spikes (over all tests)")
legend("Learner", "Control")


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
%check_parameter = [];
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
                            
                            % remove activity at stimuli
                            stimulus_data_m{i,ii}(j,9:11) = 0;
                            stimulus_data_high_fs_m{i,ii}(j,90:110) = 0;

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
            
            if j+neuron_counter_m == 831
                disp(peaks_to_keep_m)
            end
            
            for jj = 1:size(peaks_to_keep_m,2)
                [~, adjustment] = max(neuron_spikes_m_hfs(neuron_counter_m+j,peaks_to_keep_m(jj) * 10 - 10:peaks_to_keep_m(jj) * 10 + 10));
                cur_peak = peaks_to_keep_m(jj)*10 + adjustment - 10 - 1;
                latencies = [];
                for kk = 1:size(stimulus_data_high_fs_m,2)
                    if ~isempty(stimulus_data_high_fs_m{i,kk})
                        if cur_peak + 50 > 700
                            peaks_this_interval = find(stimulus_data_high_fs_m{i,kk}(j,cur_peak-50:end));
                        else
                            peaks_this_interval = find(stimulus_data_high_fs_m{i,kk}(j,cur_peak-50:cur_peak+50));
                            latencies = [latencies, peaks_this_interval];
                        end
                    end
                end
                if 16 > std(latencies)
                    secondary_neurons_m = [secondary_neurons_m; {neuron_counter_m+j, peaks_to_keep_m(jj)}];
                else
                    others_std_larger = [others_std_larger; {neuron_counter_m+j, peaks_to_keep_m(jj)}];
                end
                %check_parameter = [neuron_counter_m+j, peaks_to_keep_m(jj), cur_peak, 0, latencies];                          %%%%%%%%%%
                %tmp = neuron_counter_m + j;
                %save("/scratch/mathiass-takeokalab/01/latencies" + tmp, "check_parameter", "-v7.3");
            end


            if 0
                for jj = 1:size(peaks_to_keep_m,2)
                    [~, adjustment] = max(neuron_spikes_m_hfs(neuron_counter_m+j,peaks_to_keep_m(jj) * 10 - 10:peaks_to_keep_m(jj) * 10 + 10));
                    cur_peak = peaks_to_keep_m(jj)*10 + adjustment - 10 - 1;
                    peak_area_m = [];
                    max_peaks = [];
                    if jj ~= 0
                        if cur_peak + max_peak_distance*10 > interval_size*10 % edge case
                            peak_area_m = neuron_spikes_m_hfs(neuron_counter_m+j,cur_peak - max_peak_distance*10:end);                               % "j" should be "neuron_counter + j" in general case
                        elseif  cur_peak - max_peak_distance*10 < 1 % edge case
                            peak_area_m = neuron_spikes_m_hfs(neuron_counter_m+j,1:cur_peak + max_peak_distance*10);                                     % "j" should be "neuron_counter + j" in general case
                        else % normal case
                            peak_area_m = neuron_spikes_m_hfs(neuron_counter_m+j,cur_peak - max_peak_distance*10:cur_peak+max_peak_distance*10);     % "j" should be "neuron_counter + j" in general case
                        end
                    end
 
                    closest_peaks_m = zeros(1,sum(peak_area_m));
                    cur_index = 1;

                    for kk = 1:size(peak_area_m,2)
                        cur_latency = kk - max_peak_distance*10 - 1;
                        closest_peaks_m(cur_index:cur_index + peak_area_m(1,kk)) = cur_latency;
                        cur_index = cur_index + peak_area_m(1,kk);
                    end
                    %peak_area_m = peak_area_m - round(mean_before_m(j,1));                                                         % "j" should be "neuron_counter + j" in general case       didn't do this for now bc this would be wrong (not the same average)
                    %peak_area_m(peak_area_m<0) = 0;
                    % for kk = 1:size(peak_area_m,2)
                    %     if cur_peak - max_peak_distance < 1
                    %         latency = kk-max_peak_distance*10-1 - (cur_peak - max_peak_distance*10-1);
                    %     else
                    %         latency = kk-max_peak_distance*10-1;
                    %     end
                    %     closest_peaks_m(cur_index:cur_index+peak_area_m(kk)-1) = latency;
                    %     cur_index = cur_index + peak_area_m(kk);
                    % end
                    
                    check_parameter = [neuron_counter_m+j, peaks_to_keep_m(jj), closest_peaks_m];                          %%%%%%%%%%
                    tmp = neuron_counter_m + j;
                    save("/scratch/mathiass-takeokalab/01/check_parameter_" + tmp, "check_parameter", "-v7.3");
                    
                    if 10 > std(closest_peaks_m)
                        secondary_neurons_m = [secondary_neurons_m; {neuron_counter_m+j, peaks_to_keep_m(jj)}];
                    else
                        others_std_larger = [others_std_larger; {neuron_counter_m+j, peaks_to_keep_m(jj)}];
                    end

                end
            end
            if 0
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
                        secondary_neurons_m = [secondary_neurons_m; {neuron_counter_m+j, peaks_to_keep_m(jj)}];
                    else
                        others_std_larger = [others_std_larger; {neuron_counter_m+j, peaks_to_keep_m(jj)}];
                    end
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



%% validation
% latencies_29 = load("X:\Mathias\10kfs\01\latencies29");latencies_29 = latencies_29.check_parameter;
% peak=569;figure;subplot(2,1,1); scatter(latencies_29(1,5:end)+peak-50,1:length(latencies_29(1,5:end)));xlim([0 700]);hold on;xline(peak-50);xline(peak+51)
% subplot(2,1,2);plot(neuron_spikes_m(29,:));hold on;plot(neuron_spikes_m_hfs(29,:)


% false positives:
% 462
% 532
% 541
% 561
% 658
% 831




