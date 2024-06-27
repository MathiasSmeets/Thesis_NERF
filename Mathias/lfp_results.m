clear;clc;close all; %% OLD FUNCTION, GO TO lfp_results_adjusted.m FOR UP TO DATE FUNCTION

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
average_power_during_replay_m = zeros(1,9);
average_power_during_no_replay_m = zeros(1,9);
for i = setdiff(1:9,mouse_to_exclude_m)
    disp(i)
    % ripple_m = load("X:/Mathias/switch_data/LF_signals/ripple2_"+i+"m.mat"); ripple_m = ripple_m.ripple;
    power_m = load("X:/Mathias/switch_data/LF_signals/power_"+i+"m.mat"); power_m = power_m.power;
    % se = strel('line', 5, 90);
    % ripple_m = imdilate(ripple_m, se);
    % ripple_m = imerode(ripple_m, se);

    after_m = raw_after_m{i};
    threshold_m = threshold_after_m(i);
    
    replay_m = (after_m>threshold_m);

    
    % replay_ripple_m{i} = 0;
    % replay_no_ripple_m{i} = 0;
    % no_replay_ripple_m{i} = 0;
    % no_replay_no_ripple_m{i} = 0;

    power_during_replay_m = zeros(floor(numel(replay_m)*2.5),1);
    power_during_no_replay_m = zeros(floor(numel(replay_m)*2.5),1);
    replay_counter = 0;
    no_replay_counter = 0;
    for j = 1:numel(replay_m)
        % j is in ms
        % translate this to the sampling of LPF (2.5kHz)
        % ripple_index_m = ceil(j/2.5);
        % if replay_m(j) && ripple_m(ripple_index_m)
        %     replay_ripple_m{i} = replay_ripple_m{i} + 1;
        % elseif replay_m(j) && ~ripple_m(ripple_index_m)
        %     replay_no_ripple_m{i} = replay_no_ripple_m{i} + 1;
        % elseif ~replay_m(j) && ripple_m(ripple_index_m)
        %     no_replay_ripple_m{i} = no_replay_ripple_m{i} + 1;
        % elseif ~replay_m(j) && ~ripple_m(ripple_index_m)
        %     no_replay_no_ripple_m{i} = no_replay_no_ripple_m{i} + 1;
        % end
        if replay_m(j)
            % because of the difference in sampling frequency:
            % 1 --> 1*2.5 = 2.5 =   [1, 2, 3]
            % 2 --> 2*2.5 = 5 =     [4, 5]
            % 3 --> 3*2.5 = 7.5 =   [6, 7, 8]
            % 4 --> 4*2.5 = 10 =    [9, 10]
            % ...
            if mod(j,2) == 1
                power_during_replay_m(replay_counter+1) = power_m(ceil(j*2.5)-2);
                power_during_replay_m(replay_counter+2) = power_m(ceil(j*2.5)-1);
                power_during_replay_m(replay_counter+3) = power_m(ceil(j*2.5));
                replay_counter = replay_counter+3;
            else
                power_during_replay_m(replay_counter+1) = power_m(ceil(j*2.5)-2);
                power_during_replay_m(replay_counter+2) = power_m(ceil(j*2.5)-1);
                replay_counter = replay_counter+2;
            end
        else
            if mod(j,2) == 1
                power_during_no_replay_m(no_replay_counter+1) = power_m(ceil(j*2.5)-2);
                power_during_no_replay_m(no_replay_counter+2) = power_m(ceil(j*2.5)-1);
                power_during_no_replay_m(no_replay_counter+2) = power_m(ceil(j*2.5)-1);
                no_replay_counter = no_replay_counter+3;
            else
                power_during_no_replay_m(no_replay_counter+1) = power_m(ceil(j*2.5)-2);
                power_during_no_replay_m(no_replay_counter+2) = power_m(ceil(j*2.5)-1);
                no_replay_counter = no_replay_counter+2;
            end
        end
    end
    power_during_replay_m = power_during_replay_m(1:replay_counter);
    power_during_no_replay_m = power_during_no_replay_m(1:no_replay_counter);
    average_power_during_replay_m(i) = mean(power_during_replay_m)/385;
    average_power_during_no_replay_m(i) = mean(power_during_no_replay_m)/385;
    
    % probability_replay_ripple_m(i) = replay_ripple_m{i} / (numel(find(ripple_m))*10);
    % probability_replay_no_ripple_m(i) = replay_no_ripple_m{i} / (numel(find(~ripple_m))*10);
end

%% control
raw_after_y = load("X:\Mathias\switch_data\correlations\raw_after_y.mat"); raw_after_y = raw_after_y.cur_correlation_after;
threshold_after_y = load("X:\Mathias\switch_data\correlations\threshold_after_y.mat"); threshold_after_y = threshold_after_y.threshold_after;

replay_ripple_y = cell(1,9);
replay_no_ripple_y = cell(1,9);
no_replay_ripple_y = cell(1,9);
no_replay_no_ripple_y = cell(1,9);
mouse_to_exclude_m = []; % m
mouse_to_exclude_y = [2,9]; % y
% compare probability replay occurring during ripple versus non-ripple
probability_replay_ripple_y = zeros(1,9);
probability_replay_no_ripple_y = zeros(1,9);
for i = setdiff(1:9,mouse_to_exclude_y)
    disp(i)
    ripple_y = load("X:/Mathias/switch_data/LF_signals/ripple2_"+i+"y.mat"); ripple_y = ripple_y.ripple;

    se = strel('line', 5, 90);
    ripple_y = imdilate(ripple_y, se);
    ripple_y = imerode(ripple_y, se);

    after_y = raw_after_y{i};
    threshold_y = threshold_after_y(i);
    
    replay_y = (after_y>threshold_y);

    
    replay_ripple_y{i} = 0;
    replay_no_ripple_y{i} = 0;
    no_replay_ripple_y{i} = 0;
    no_replay_no_ripple_y{i} = 0;

    for j = 1:numel(replay_y)
        % j is in ms
        % translate this to the sampling of LPF (2.5kHz)
        ripple_index_y = ceil(j/2.5);
        if replay_y(j) && ripple_y(ripple_index_y)
            replay_ripple_y{i} = replay_ripple_y{i} + 1;
        elseif replay_y(j) && ~ripple_y(ripple_index_y)
            replay_no_ripple_y{i} = replay_no_ripple_y{i} + 1;
        elseif ~replay_y(j) && ripple_y(ripple_index_y)
            no_replay_ripple_y{i} = no_replay_ripple_y{i} + 1;
        elseif ~replay_y(j) && ~ripple_y(ripple_index_y)
            no_replay_no_ripple_y{i} = no_replay_no_ripple_y{i} + 1;
        end
    end
    probability_replay_ripple_y(i) = replay_ripple_y{i} / (numel(find(ripple_y))*10);
    probability_replay_no_ripple_y(i) = replay_no_ripple_y{i} / (numel(find(~ripple_y))*10);
end

%% figures
figure
boxplot([probability_replay_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))',probability_replay_ripple_m(setdiff(1:9,mouse_to_exclude_m))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(probability_replay_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1)),probability_replay_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
scatter(ones(size(probability_replay_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1))*2,probability_replay_ripple_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
line([ones(size(probability_replay_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1),1), ones(size(probability_replay_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1),1)*2]',[probability_replay_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))', probability_replay_ripple_m(setdiff(1:9,mouse_to_exclude_m))']','Color','green')
p_m = signrank(probability_replay_no_ripple_m(setdiff(1:9,mouse_to_exclude_m))', probability_replay_ripple_m(setdiff(1:9,mouse_to_exclude_m))');
ylabel('Probability Replay Occurring')
figure
boxplot([probability_replay_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))',probability_replay_ripple_y(setdiff(1:9,mouse_to_exclude_y))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(probability_replay_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1)),probability_replay_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))', 'filled', 'blue')
scatter(ones(size(probability_replay_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1))*2,probability_replay_ripple_y(setdiff(1:9,mouse_to_exclude_y))', 'filled', 'blue')
line([ones(size(probability_replay_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1),1), ones(size(probability_replay_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1),1)*2]',[probability_replay_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))', probability_replay_ripple_y(setdiff(1:9,mouse_to_exclude_y))']','Color','green')
p_y = signrank(probability_replay_no_ripple_y(setdiff(1:9,mouse_to_exclude_y))', probability_replay_ripple_y(setdiff(1:9,mouse_to_exclude_y))');
ylabel('Probability Replay Occurring')
disp("p-value learner: " + p_m)
disp("p-value control: " + p_y)
disp("OLD FUNCTION, GO TO lfp_results_adjusted.m FOR UP TO DATE FUNCTION")

