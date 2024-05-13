clear;clc;close all;


for i = setdiff(1:9,[2,8,9])
    disp(i)
    load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_filtered_total_spindle"+i+"m.mat")
    spindle = zeros(size(data));

    %%
    for j = 1:size(data,2)
        accepted = zeros(1,size(data,1));

        % wavelet transform
        [wavelet_coefficients, ~] = cwt(data(:, j), 'amor', 2500, 'FrequencyLimits', [10 16], 'VoicesPerOctave', 8);
        abs_wavelets = abs(wavelet_coefficients);
        mean_abs_wavelets = mean(abs_wavelets);
        
        % determine threshold & find spindle
        accepted(1:250)=1;
        for k = 250:size(data,1)
            useful_values = mean_abs_wavelets(accepted==1);
            previous_data = useful_values(end-249:end);
            cur_threshold = mean(abs(previous_data));
            if mean_abs_wavelets(:,k)>cur_threshold*4.5
                spindle(k,j) = 1;
            else
                accepted(k) = 1;
            end
        end
    end
    %%
    save("/scratch/mathiass-takeokalab/01/spindle_"+i+".mat", "ripple", "-v7.3");
end


[wavelet_coefficients, frequencies] = cwt(EEG_data(i, :), 'amor', 2500, 'FrequencyLimits', [10 16], 'VoicesPerOctave', 8);