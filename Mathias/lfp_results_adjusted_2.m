clear;clc;close all;

%% learner
raw_after_m = load("X:\Mathias\switch_data\correlations\raw_after_m.mat"); raw_after_m = raw_after_m.cur_correlation_after;
threshold_after_m = load("X:\Mathias\switch_data\correlations\threshold_after_m.mat"); threshold_after_m = threshold_after_m.threshold_after;

% replay_ripple_m = cell(1,9);
% replay_no_ripple_m = cell(1,9);
% no_replay_ripple_m = cell(1,9);
% no_replay_no_ripple_m = cell(1,9);
mouse_to_exclude_m = []; % m
mouse_to_exclude_y = [2,9]; % y
% compare probability replay occurring during ripple versus non-ripple
% probability_replay_ripple_m = zeros(1,9);
% probability_replay_no_ripple_m = zeros(1,9);
average_replay_during_ripple_m = zeros(1,9);
average_replay_during_no_ripple_m = zeros(1,9);
for i = setdiff(1:9,mouse_to_exclude_m)
    disp(i)
    ripple_m = load("X:/Mathias/switch_data/LF_signals/ripple2_"+i+".mat"); ripple_m = ripple_m.ripple;
    % power_m = load("X:/Mathias/switch_data/LF_signals/power_"+i+"m.mat"); power_m = power_m.power;
    se = strel('line', 10, 90);
    ripple_m = imdilate(ripple_m, se);
    ripple_m = imerode(ripple_m, se);

    after_m = raw_after_m{i};
    threshold_m = threshold_after_m(i);
    
    replay_m = after_m;
    replay_m(replay_m<=threshold_m) = 0;
    disp(mean(replay_m))

    replay_during_ripple_m{i} = 0;
    replay_during_noripple_m{i} = 0;
    counter_ripple = 0;
    counter_non_ripple = 0;
    previous_replay_index = 0;
    for j = 20:numel(ripple_m) % start a bit to the right because replay is shifted
        replay_index = ceil(j/2.5)-7; % replay is shifted (half of template)
        if replay_index > numel(replay_m)
            break
        end
        if replay_index ~= previous_replay_index % because of the sampling difference, avoid taking same replay values multiple times
            if ripple_m(j)
                replay_during_ripple_m{i} = replay_during_ripple_m{i} + replay_m(replay_index);
                counter_ripple = counter_ripple+1;
            else
                replay_during_noripple_m{i} = replay_during_noripple_m{i} + replay_m(replay_index);
                counter_non_ripple = counter_non_ripple + 1;
            end
        end
        previous_replay_index = replay_index;
    end

    average_replay_during_ripple_m(i) = replay_during_ripple_m{i} / counter_ripple;
    average_replay_during_no_ripple_m(i) = replay_during_noripple_m{i} / counter_non_ripple;
end

%% control
raw_after_y = load("X:\Mathias\switch_data\correlations\raw_after_y.mat"); raw_after_y = raw_after_y.cur_correlation_after;
threshold_after_y = load("X:\Mathias\switch_data\correlations\threshold_after_y.mat"); threshold_after_y = threshold_after_y.threshold_after;

% replay_ripple_y = cell(1,9);
% replay_no_ripple_y = cell(1,9);
% no_replay_ripple_y = cell(1,9);
% no_replay_no_ripple_y = cell(1,9);
mouse_to_exclude_m = []; % m
mouse_to_exclude_y = [2,9]; % y
% compare probability replay occurring during ripple versus non-ripple
% probability_replay_ripple_y = zeros(1,9);
% probability_replay_no_ripple_y = zeros(1,9);
average_replay_during_ripple_y = zeros(1,9);
average_replay_during_no_ripple_y = zeros(1,9);
for i = setdiff(1:9,mouse_to_exclude_y)
    disp(i)
    ripple_y = load("X:/Mathias/switch_data/LF_signals/ripple2_"+i+"y.mat"); ripple_y = ripple_y.ripple;
    % power_y = load("X:/Mathias/switch_data/LF_signals/power_"+i+"m.mat"); power_y = power_y.power;
    se = strel('line', 5, 90);
    ripple_y = imdilate(ripple_y, se);
    ripple_y = imerode(ripple_y, se);

    after_y = raw_after_y{i};
    threshold_y = threshold_after_y(i);
    
    replay_y = after_y;
    replay_y(replay_y<=threshold_y) = 0;

    replay_during_ripple_y{i} = 0;
    replay_during_noripple_y{i} = 0;
    counter_ripple = 0;
    counter_non_ripple = 0;
    previous_replay_index = 0;
    for j = 20:numel(ripple_y) % start a bit to the right because replay is shifted
        replay_index = ceil(j/2.5)-7; % replay is shifted (half of template)
        if replay_index > numel(replay_y)
            break
        end
        if replay_index ~= previous_replay_index % because of the sampling difference, avoid taking same replay values multiple times
            if ripple_y(j)
                replay_during_ripple_y{i} = replay_during_ripple_y{i} + replay_y(replay_index);
                counter_ripple = counter_ripple+1;
            else
                replay_during_noripple_y{i} = replay_during_noripple_y{i} + replay_y(replay_index);
                counter_non_ripple = counter_non_ripple + 1;
            end
        end
        previous_replay_index = replay_index;
    end

    average_replay_during_ripple_y(i) = replay_during_ripple_y{i} / counter_ripple;
    average_replay_during_no_ripple_y(i) = replay_during_noripple_y{i} / counter_non_ripple;
end

%% figures

figure
boxplot([average_replay_during_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))',average_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(average_replay_during_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1)),average_replay_during_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
scatter(ones(size(average_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1))*2,average_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
line([ones(size(average_replay_during_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1),1), ones(size(average_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1),1)*2]',[average_replay_during_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))', average_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))']','Color','green')
p_m = signrank(average_replay_during_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))', average_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))');
ylabel('Average Replay Strength')
disp("p-value learner: " + p_m)

figure
boxplot([average_replay_during_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))',average_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(average_replay_during_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1)),average_replay_during_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))', 'filled', 'blue')
scatter(ones(size(average_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1))*2,average_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))', 'filled', 'blue')
line([ones(size(average_replay_during_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1),1), ones(size(average_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1),1)*2]',[average_replay_during_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))', average_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))']','Color','green')
p_y = signrank(average_replay_during_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))', average_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))');
ylabel('Average Replay Strength')
disp("p-value control: " + p_y)
