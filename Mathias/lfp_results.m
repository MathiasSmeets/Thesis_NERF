clear;clc;close all;

raw_after = load("X:\Mathias\switch_data\correlations\raw_after_m.mat"); raw_after = raw_after.cur_correlation_after;
load("X:\Mathias\switch_data\correlations\threshold_after_m.mat");%threshold_after

replay_ripple = cell(1,9);
replay_no_ripple = cell(1,9);
no_replay_ripple = cell(1,9);
no_replay_no_ripple = cell(1,9);
for i = setdiff(1:9,[2,8,9])
    disp(i)
    load("X:/Mathias/switch_data/LF_signals/ripple_"+i+".mat")%ripple
    after = raw_after{i};
    threshold = threshold_after(i);
    
    replay = (after>threshold);

    
    replay_ripple{i} = 0;
    replay_no_ripple{i} = 0;
    no_replay_ripple{i} = 0;
    no_replay_no_ripple{i} = 0;

    for j = 1:numel(replay)
        ripple_index = ceil(j/10);
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
end
% compare probability replay occurring during ripple versus non-ripple
probability_replay_ripple = zeros(1,9);
probability_replay_no_ripple = zeros(1,9);
for i = setdiff(1:9,[2,8,9])
    probability_replay_ripple(i) = replay_ripple{i} / (numel(find(ripple))*10);
    probability_replay_no_ripple(i) = replay_no_ripple{i} / (numel(find(~ripple))*10);
end

figure
boxplot([probability_replay_no_ripple(setdiff(1:9,[2,8,9]))',probability_replay_ripple(setdiff(1:9,[2,8,9]))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(probability_replay_no_ripple(setdiff(1:9,[2,8,9]))',1)),probability_replay_no_ripple(setdiff(1:9,[2,8,9]))', 'filled', 'blue')
scatter(ones(size(probability_replay_ripple(setdiff(1:9,[2,8,9]))',1))*2,probability_replay_ripple(setdiff(1:9,[2,8,9]))', 'filled', 'blue')
line([ones(size(probability_replay_no_ripple(setdiff(1:9,[2,8,9]))',1),1), ones(size(probability_replay_ripple(setdiff(1:9,[2,8,9]))',1),1)*2]',[probability_replay_no_ripple(setdiff(1:9,[2,8,9]))', probability_replay_ripple(setdiff(1:9,[2,8,9]))']','Color','green')
p_m = signrank(probability_replay_no_ripple(setdiff(1:9,[2,8,9]))', probability_replay_ripple(setdiff(1:9,[2,8,9]))');
ylabel('Probability Replay Occurring')


