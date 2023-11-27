
stimulus_data_m_HIGH_FS = load("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\data_after_stimulus_m_HIGHFS_NEW.mat");
stimulus_data_m_HIGH_FS = stimulus_data_m_HIGH_FS.after_stimulus_data_m;

max_peak_distance = 50;
neuron_counter_m = 0;
%124-138

% 124-130: no
% 131: 12    26    29
neuron_spikes_m_HIGH_FS = zeros(15,700);
for i = 1:size(stimulus_data_m_HIGH_FS,1)
    for j = 1:size(stimulus_data_m_HIGH_FS,2)
        if ~isempty(stimulus_data_m_HIGH_FS{i,j})
            neuron_spikes_m_HIGH_FS = neuron_spikes_m_HIGH_FS + stimulus_data_m_HIGH_FS{i,j};
        end
    end
end




for j = 1:15
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