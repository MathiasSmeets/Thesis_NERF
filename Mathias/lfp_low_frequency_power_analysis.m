clear;clc;close all;

raw_after_m = load("X:\Mathias\switch_data\correlations\raw_after_m.mat"); raw_after_m = raw_after_m.cur_correlation_after;
threshold_after_m = load("X:\Mathias\switch_data\correlations\threshold_after_m.mat"); threshold_after_m = threshold_after_m.threshold_after;
learner = 1; % change to zero if control

%%
% average_powers_segment = cell(1,9);
% average_replay_segment = cell(1,9);
% for i = 1:9
%     cur_fig = openfig("X:\Mathias\01\01\powerfigure_delta_"+i+".fig");
%     cur_ax_objects = cur_fig.Children;
%     cur_data = cur_ax_objects.ColorData;
% 
%     after_m = raw_after_m{i};
%     threshold_m = threshold_after_m(i);
% 
%     replay_m = after_m;
%     replay_m(replay_m<=threshold_m) = 0;
% 
%     segment_size = fix(numel(replay_m) / size(cur_data,2));           %% should this always be an integer?
%     average_powers_segment{i} = zeros(size(cur_data,2), 1);
%     average_replay_segment{i} = zeros(size(cur_data,2), 1);
% 
%     for j = 1:size(cur_data,2)
%         average_powers_segment{i}(j) = mean(cur_data(:,j));
%         replay_indices = (j-1)*segment_size + 1 : j*segment_size;       %% idem
%         average_replay_segment{i}(j) = mean(replay_m(replay_indices));
%     end
% end

%% or do something with a threshold:
average_replay_segment = cell(1,9);
average_replay_during_power = zeros(1,9);
average_replay_during_no_power = zeros(1,9);

for i = 1:9
    cur_fig = openfig("X:\Mathias\switch_data\LF_signals\powerfigure_delta_"+i+".fig");
    cur_ax_objects = cur_fig.Children;
    cur_data = cur_ax_objects.ColorData;
    meanall = mean(cur_data,'all');
    stdall = std(cur_data(:));
    %zscored_data = zscore(cur_data,[],2);
    zscored_data = (cur_data - meanall) / stdall;
    figure;heatmap(zscored_data,'Colormap',spring, 'GridVisible','off','ColorLimits',[-3 3])

    after_m = raw_after_m{i};
    threshold_m = threshold_after_m(i);
    
    replay_m = after_m;
    replay_m(replay_m<=threshold_m) = 0;
   
    channels_of_interest = get_channels_of_interest(i, 1, learner);

    segment_size = fix(numel(replay_m) / size(cur_data,2));           %% we'll see if this works
    average_replay_segment{i} = zeros(size(cur_data,2), 1);
    % if power is above 30 for at least 50 channels
    [~, highest_segments] = maxk(mean(zscored_data(channels_of_interest,:)),50);
    power = zeros(1,size(cur_data,2));
    for j = 1:size(cur_data,2)
        %if (numel(find(zscored_data(cur_data,j) > 30))) >= numel(channels_of_interest)*0.8
        %if (numel(find(zscored_data(channels_of_interest,j) > 1))) >= numel(channels_of_interest)*0.8
        if ismember(j, highest_segments)
            power(j) = 1;
        end
        replay_indices = (j-1)*segment_size + 1 : j*segment_size;       %% idem
        average_replay_segment{i}(j) = mean(replay_m(replay_indices));
    end
    
    replay_power = 0;
    replay_no_power = 0;
    for j = 1:numel(power)
        if power(j)
            replay_power = replay_power + average_replay_segment{i}(j);
        else
            replay_no_power = replay_no_power + average_replay_segment{i}(j);
        end
    end
    average_replay_during_power(i) = replay_power / numel(find(power));
    average_replay_during_no_power(i) = replay_no_power / numel(find(~power));
    disp("power: " + i)
    disp(numel(find(power)))
    disp("no power: " + i)
    disp(numel(find(~power)))
end
p_m = signrank(average_replay_during_power, average_replay_during_no_power);
disp(p_m)