clear;clc;close all;

raw_after = load("X:\Mathias\switch_data\correlations\raw_after_m.mat"); raw_after = raw_after.cur_correlation_after;
load("X:\Mathias\switch_data\correlations\threshold_after_m.mat");%threshold_after

replay_ripple = cell(1,9);
replay_no_ripple = cell(1,9);
no_replay_ripple = cell(1,9);
no_replay_no_ripple = cell(1,9);
mouse_to_exclude = [];
% compare probability replay occurring during ripple versus non-ripple
probability_ripple_during_replay = zeros(1,9);
probability_noripple_during_replay = zeros(1,9);
for i = setdiff(1:9,mouse_to_exclude)
    disp(i)
    load("X:/Mathias/switch_data/LF_signals/ripple2_"+i+".mat")%ripple

    se = strel('line', 5, 90);
    ripple = imdilate(ripple, se);
    ripple = imerode(ripple, se);

    after = raw_after{i};
    threshold = threshold_after(i);
    
    replay = (after>threshold);

    transitions = diff([0, replay, 0]);
    % Find the starting indices of blocks of ones
    starts = find(transitions == 1);
    % Find the ending indices of blocks of ones
    ends = find(transitions == -1) - 1;
    % Combine the start and end indices into a 2xY array
    blocks = [starts; ends];

    ripple_during_replay{i} = 0;
    noripple_during_replay{i} = 0;
    replay_counter = 0;
    for j = 1:size(blocks,2)

        starttime = ceil(blocks(1,j)*2.5);
        endtime = ceil(blocks(2,j)*2.5);
        if ~isempty(find(ripple(starttime:endtime), 1))
            ripple_during_replay{i} = ripple_during_replay{i} + 1;
        elseif isempty(find(ripple(starttime:endtime), 1))
            noripple_during_replay{i} = noripple_during_replay{i} + 1;
        end
        replay_counter = replay_counter + blocks(2,j)-blocks(1,j)+1;
    end

    probability_ripple_during_replay(i) = ripple_during_replay{i} / numel(find(ripple));
    probability_noripple_during_replay(i) = noripple_during_replay{i} / numel(find(~ripple));



    transitions_ripple = diff([0, ripple', 0]);
    % Find the starting indices of blocks of ones
    starts_ripple = find(transitions_ripple == 1);
    % Find the ending indices of blocks of ones
    ends_ripple = find(transitions_ripple == -1) - 1;
    % Combine the start and end indices into a 2xY array
    blocks_ripple = [starts_ripple; ends_ripple];

    replay_during_ripple{i} = 0;
    noreplay_during_ripple{i} = 0;
    ripple_counter = 0;
    for j = 1:size(blocks_ripple,2)
        starttime = ceil(blocks_ripple(1,j)/2.5);
        endtime = min(ceil(blocks_ripple(2,j)/2.5), numel(replay));
        if ~isempty(find(replay(starttime:endtime), 1))
            replay_during_ripple{i} = replay_during_ripple{i} + 1;
        elseif isempty(find(replay(starttime:endtime), 1))
            noreplay_during_ripple{i} = noreplay_during_ripple{i} + 1;
        end
        ripple_counter = ripple_counter + blocks_ripple(2,j)-blocks_ripple(1,j)+1;
    end

    probability_replay_during_ripple(i) = replay_during_ripple{i} / numel(find(replay));
    probability_noreplay_during_ripple(i) = noreplay_during_ripple{i} / numel(find(~replay));

end

figure
boxplot([probability_noripple_during_replay(setdiff(1:9,mouse_to_exclude))',probability_ripple_during_replay(setdiff(1:9,mouse_to_exclude))'],'Labels',{'Not During Replay','During Replay'})
hold on
scatter(ones(size(probability_noripple_during_replay(setdiff(1:9,mouse_to_exclude))',1)),probability_noripple_during_replay(setdiff(1:9,mouse_to_exclude))', 'filled', 'blue')
scatter(ones(size(probability_ripple_during_replay(setdiff(1:9,mouse_to_exclude))',1))*2,probability_ripple_during_replay(setdiff(1:9,mouse_to_exclude))', 'filled', 'blue')
line([ones(size(probability_noripple_during_replay(setdiff(1:9,mouse_to_exclude))',1),1), ones(size(probability_ripple_during_replay(setdiff(1:9,mouse_to_exclude))',1),1)*2]',[probability_noripple_during_replay(setdiff(1:9,mouse_to_exclude))', probability_ripple_during_replay(setdiff(1:9,mouse_to_exclude))']','Color','green')
p_m = signrank(probability_noripple_during_replay(setdiff(1:9,mouse_to_exclude))', probability_ripple_during_replay(setdiff(1:9,mouse_to_exclude))');
ylabel('Probability Ripple Occurring')
disp("p-value: " + p_m)

figure
boxplot([probability_noreplay_during_ripple(setdiff(1:9,mouse_to_exclude))',probability_replay_during_ripple(setdiff(1:9,mouse_to_exclude))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(probability_noreplay_during_ripple(setdiff(1:9,mouse_to_exclude))',1)),probability_noreplay_during_ripple(setdiff(1:9,mouse_to_exclude))', 'filled', 'blue')
scatter(ones(size(probability_replay_during_ripple(setdiff(1:9,mouse_to_exclude))',1))*2,probability_replay_during_ripple(setdiff(1:9,mouse_to_exclude))', 'filled', 'blue')
line([ones(size(probability_noreplay_during_ripple(setdiff(1:9,mouse_to_exclude))',1),1), ones(size(probability_replay_during_ripple(setdiff(1:9,mouse_to_exclude))',1),1)*2]',[probability_noreplay_during_ripple(setdiff(1:9,mouse_to_exclude))', probability_replay_during_ripple(setdiff(1:9,mouse_to_exclude))']','Color','green')
p_m = signrank(probability_noreplay_during_ripple(setdiff(1:9,mouse_to_exclude))', probability_replay_during_ripple(setdiff(1:9,mouse_to_exclude))');
ylabel('Probability Replay Occurring')
disp("p-value: " + p_m)


