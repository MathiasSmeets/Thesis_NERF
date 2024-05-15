clear;clc;close all;

raw_after = load("X:\Mathias\switch_data\correlations\raw_after_m.mat"); raw_after = raw_after.cur_correlation_after;
load("X:\Mathias\switch_data\correlations\threshold_after_m.mat");%threshold_after

replay_ripple = cell(1,9);
replay_no_ripple = cell(1,9);
no_replay_ripple = cell(1,9);
no_replay_no_ripple = cell(1,9);
mouse_to_exclude = [];
% compare probability replay occurring during ripple versus non-ripple
probability_replay_ripple = zeros(1,9);
probability_replay_no_ripple = zeros(1,9);
for i = setdiff(1:9,mouse_to_exclude)
    disp(i)
    load("X:/Mathias/switch_data/LF_signals/ripple2_"+i+".mat")%ripple

    se = strel('line', 5, 90);
    ripple = imdilate(ripple, se);
    ripple = imerode(ripple, se);

    after = raw_after{i};
    threshold = threshold_after(i);
    
    replay = (after>threshold);

    
    replay_ripple{i} = 0;
    replay_no_ripple{i} = 0;
    no_replay_ripple{i} = 0;
    no_replay_no_ripple{i} = 0;

    for j = 1:numel(replay)
        % j is in ms
        % translate this to the sampling of LPF (2.5kHz)
        ripple_index = ceil(j/2.5);
        if replay(j) && ripple(ripple_index)
            replay_ripple{i} = replay_ripple{i} + 1;
        elseif replay(j) && ~ripple(ripple_index)
            replay_no_ripple{i} = replay_no_ripple{i} + 1;
        elseif ~replay(j) && ripple(ripple_index)
            no_replay_ripple{i} = no_replay_ripple{i} + 1;
        elseif ~replay(j) && ~ripple(ripple_index)
            no_replay_no_ripple{i} = no_replay_no_ripple{i} + 1;
        end
    end
    probability_replay_ripple(i) = replay_ripple{i} / (numel(find(ripple))*10);
    probability_replay_no_ripple(i) = replay_no_ripple{i} / (numel(find(~ripple))*10);
end

figure
boxplot([probability_replay_no_ripple(setdiff(1:9,mouse_to_exclude))',probability_replay_ripple(setdiff(1:9,mouse_to_exclude))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(probability_replay_no_ripple(setdiff(1:9,mouse_to_exclude))',1)),probability_replay_no_ripple(setdiff(1:9,mouse_to_exclude))', 'filled', 'blue')
scatter(ones(size(probability_replay_ripple(setdiff(1:9,mouse_to_exclude))',1))*2,probability_replay_ripple(setdiff(1:9,mouse_to_exclude))', 'filled', 'blue')
line([ones(size(probability_replay_no_ripple(setdiff(1:9,mouse_to_exclude))',1),1), ones(size(probability_replay_ripple(setdiff(1:9,mouse_to_exclude))',1),1)*2]',[probability_replay_no_ripple(setdiff(1:9,mouse_to_exclude))', probability_replay_ripple(setdiff(1:9,mouse_to_exclude))']','Color','green')
p_m = signrank(probability_replay_no_ripple(setdiff(1:9,mouse_to_exclude))', probability_replay_ripple(setdiff(1:9,mouse_to_exclude))');
ylabel('Probability Replay Occurring')
disp("p-value: " + p_m)

