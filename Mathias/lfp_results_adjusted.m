clear;clc;close all;


%% learner
raw_after_m = load("X:\Mathias\switch_data\correlations\raw_after_m.mat"); raw_after_m = raw_after_m.cur_correlation_after;
threshold_after_m = load("X:\Mathias\switch_data\correlations\threshold_after_m.mat"); threshold_after_m = threshold_after_m.threshold_after;
mouse_to_exclude_m = [];
% compare probability replay occurring during ripple versus non-ripple
probability_ripple_during_replay_m = zeros(1,9);
probability_noripple_during_replay_m = zeros(1,9);
for i = setdiff(1:9,mouse_to_exclude_m)
    disp(i)
    ripple_m = load("X:/Mathias/switch_data/LF_signals/ripple2_"+i+".mat"); ripple_m = ripple_m.ripple;

    se = strel('line', 5, 90);
    ripple_m = imdilate(ripple_m, se);
    ripple_m = imerode(ripple_m, se);

    after_m = raw_after_m{i};
    threshold_m = threshold_after_m(i);
    
    replay_m = (after_m>threshold_m);

    transitions_m = diff([0, replay_m, 0]);
    % Find the starting indices of blocks of ones
    starts_m = find(transitions_m == 1);
    % Find the ending indices of blocks of ones
    ends_m = find(transitions_m == -1) - 1;
    % Combine the start and end indices into a 2xY array
    blocks_m = [starts_m; ends_m];
    blocks_m = blocks_m - 7.5;

    ripple_during_replay_m{i} = 0;
    noripple_during_replay_m{i} = 0;
    for j = 1:size(blocks_m,2)
        starttime = ceil(blocks_m(1,j)*2.5);
        endtime = ceil(blocks_m(2,j)*2.5);
        if ~isempty(find(ripple_m(starttime:endtime), 1))
            ripple_during_replay_m{i} = ripple_during_replay_m{i} + 1;
        elseif isempty(find(ripple_m(starttime:endtime), 1))
            noripple_during_replay_m{i} = noripple_during_replay_m{i} + 1;
        end
    end

    probability_ripple_during_replay_m(i) = (ripple_during_replay_m{i} / numel(find(ripple_m)));
    probability_noripple_during_replay_m(i) = (noripple_during_replay_m{i} / numel(find(~ripple_m)));



    transitions_ripple_m = diff([0, ripple_m', 0]);
    % Find the starting indices of blocks of ones
    starts_ripple_m = find(transitions_ripple_m == 1);
    % Find the ending indices of blocks of ones
    ends_ripple_m = find(transitions_ripple_m == -1) - 1;
    % Combine the start and end indices into a 2xY array
    blocks_ripple_m = [starts_ripple_m; ends_ripple_m];

    replay_during_ripple_m{i} = 0;
    noreplay_during_ripple_m{i} = 0;
    for j = 1:size(blocks_ripple_m,2)
        starttime = ceil(blocks_ripple_m(1,j)/2.5);
        endtime = min(ceil(blocks_ripple_m(2,j)/2.5), numel(replay_m));
        if ~isempty(find(replay_m(starttime:endtime), 1))
            replay_during_ripple_m{i} = replay_during_ripple_m{i} + 1;
        elseif isempty(find(replay_m(starttime:endtime), 1))
            noreplay_during_ripple_m{i} = noreplay_during_ripple_m{i} + 1;
        end
    end

    probability_replay_during_ripple_m(i) = (replay_during_ripple_m{i} / numel(find(replay_m)));
    probability_noreplay_during_ripple_m(i) = (noreplay_during_ripple_m{i} / numel(find(~replay_m)));
end


%% control

raw_after_y = load("X:\Mathias\switch_data\correlations\raw_after_y.mat"); raw_after_y = raw_after_y.cur_correlation_after;
threshold_after_y = load("X:\Mathias\switch_data\correlations\threshold_after_y.mat"); threshold_after_y = threshold_after_y.threshold_after;

mouse_to_exclude_y = 9;
% compare probability replay occurring during ripple versus non-ripple
probability_ripple_during_replay_y = zeros(1,9);
probability_noripple_during_replay_y = zeros(1,9);
for i = setdiff(1:9,mouse_to_exclude_y)
    disp(i)
    ripple_y = load("X:/Mathias/switch_data/LF_signals/ripple2_"+i+"y.mat"); ripple_y = ripple_y.ripple;

    se = strel('line', 5, 90);
    ripple_y = imdilate(ripple_y, se);
    ripple_y = imerode(ripple_y, se);

    after_y = raw_after_y{i};
    threshold_y = threshold_after_y(i);
    
    replay_y = (after_y>threshold_y);

    transitions_y = diff([0, replay_y, 0]);
    % Find the starting indices of blocks of ones
    starts_y = find(transitions_y == 1);
    % Find the ending indices of blocks of ones
    ends_y = find(transitions_y == -1) - 1;
    % Combine the start and end indices into a 2xY array
    blocks_y = [starts_y; ends_y];
    blocks_y = blocks_y - 7.5;

    ripple_during_replay_y{i} = 0;
    noripple_during_replay_y{i} = 0;
    for j = 1:size(blocks_y,2)
        starttime = ceil(blocks_y(1,j)*2.5);
        endtime = ceil(blocks_y(2,j)*2.5);
        if ~isempty(find(ripple_y(starttime:endtime), 1))
            ripple_during_replay_y{i} = ripple_during_replay_y{i} + 1;
        elseif isempty(find(ripple_y(starttime:endtime), 1))
            noripple_during_replay_y{i} = noripple_during_replay_y{i} + 1;
        end
    end

    probability_ripple_during_replay_y(i) = (ripple_during_replay_y{i} / numel(find(ripple_y)));
    probability_noripple_during_replay_y(i) = (noripple_during_replay_y{i} / numel(find(~ripple_y)));



    transitions_ripple_y = diff([0, ripple_y', 0]);
    % Find the starting indices of blocks of ones
    starts_ripple_y = find(transitions_ripple_y == 1);
    % Find the ending indices of blocks of ones
    ends_ripple_y = find(transitions_ripple_y == -1) - 1;
    % Combine the start and end indices into a 2xY array
    blocks_ripple_y = [starts_ripple_y; ends_ripple_y];

    replay_during_ripple_y{i} = 0;
    noreplay_during_ripple_y{i} = 0;
    for j = 1:size(blocks_ripple_y,2)
        starttime = ceil(blocks_ripple_y(1,j)/2.5);
        endtime = min(ceil(blocks_ripple_y(2,j)/2.5), numel(replay_y));
        if ~isempty(find(replay_y(starttime:endtime), 1))
            replay_during_ripple_y{i} = replay_during_ripple_y{i} + 1;
        elseif isempty(find(replay_y(starttime:endtime), 1))
            noreplay_during_ripple_y{i} = noreplay_during_ripple_y{i} + 1;
        end
    end

    probability_replay_during_ripple_y(i) = (replay_during_ripple_y{i} / numel(find(replay_y)));
    probability_noreplay_during_ripple_y(i) = (noreplay_during_ripple_y{i} / numel(find(~replay_y)));

end

%% figures
figure
boxplot([probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))',probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))',1)),probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
scatter(ones(size(probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))',1))*2,probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
line([ones(size(probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))',1),1), ones(size(probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))',1),1)*2]',[probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))', probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))']','Color','green')
p_m = signrank(probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))', probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))');
ylabel('Probability Replay Occurrs')
disp("p-value: " + p_m)

figure
boxplot([probability_noreplay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))',probability_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))'],'Labels',{'Not During Replay','During Replay'})
hold on
scatter(ones(size(probability_noreplay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1)),probability_noreplay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
scatter(ones(size(probability_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1))*2,probability_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
line([ones(size(probability_noreplay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1),1), ones(size(probability_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))',1),1)*2]',[probability_noreplay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))', probability_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))']','Color','green')
p_m = signrank(probability_noreplay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))', probability_replay_during_ripple_m(setdiff(1:9,mouse_to_exclude_m))');
ylabel('Probability Ripple Occurring')
disp("p-value: " + p_m)

figure
boxplot([probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))',probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))'],'Labels',{'Not During Ripple','During Ripple'})
hold on
scatter(ones(size(probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))',1)),probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))', 'filled', 'blue')
scatter(ones(size(probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))',1))*2,probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))', 'filled', 'blue')
line([ones(size(probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))',1),1), ones(size(probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))',1),1)*2]',[probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))', probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))']','Color','green')
p_y = signrank(probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))', probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))');
ylabel('Probability Replay Occurrs')
disp("p-value: " + p_y)

figure
boxplot([probability_noreplay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))',probability_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))'],'Labels',{'Not During Replay','During Replay'})
hold on
scatter(ones(size(probability_noreplay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1)),probability_noreplay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))', 'filled', 'blue')
scatter(ones(size(probability_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1))*2,probability_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))', 'filled', 'blue')
line([ones(size(probability_noreplay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1),1), ones(size(probability_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))',1),1)*2]',[probability_noreplay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))', probability_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))']','Color','green')
p_y = signrank(probability_noreplay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))', probability_replay_during_ripple_y(setdiff(1:9,mouse_to_exclude_y))');
ylabel('Probability Ripple Occurring')
disp("p-value: " + p_y)

figure
data2 = [probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))'; probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))'; probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))';probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))'];
g2 = [zeros(length(probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))'),1); ones(length(probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))'),1); 2*ones(length(probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))'),1); 3*ones(length(probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))'),1)];
boxplot(data2, g2, 'Labels',{'Not During Ripple', 'During Ripple','Not During Ripple','During Ripple'})
hold on
scatter(ones(size(probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m)),1)), probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m)), 'filled', 'blue')
scatter(ones(size(probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m)),1))*2, probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))', 'filled', 'blue')
line([ones(numel(probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))),1), ones(numel(probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))),1)*2]',[probability_noripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))', probability_ripple_during_replay_m(setdiff(1:9,mouse_to_exclude_m))']','Color','green')

scatter(ones(size(probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y)),1))*3, probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y)), 'filled', 'blue')
scatter(ones(size(probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y)),1))*4, probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y)), 'filled', 'blue')
line([ones(numel(probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))),1)*3, ones(numel(probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))),1)*4]',[probability_noripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))', probability_ripple_during_replay_y(setdiff(1:9,mouse_to_exclude_y))']','Color','green')
ylabel('Probability Replay Occurrs')


%% create more figures
startt = 350101 -90;
endd = 350101 +100;
ripple_m = load("X:/Mathias/switch_data/LF_signals/ripple2_1.mat"); ripple_m = ripple_m.ripple;
%begin = find(ripple_m(ceil(startt*2.5):ceil(endd*2.5)),1);
%einde = find(ripple_m(ceil(startt*2.5):ceil(endd*2.5)),1,'last');
begin = ceil(236);
einde = ceil(266);
figure
subplot(3,1,1)
cur_file = matfile("X:\Mathias\switch_data\LF_signals\LF_filtered_rec1_part1m_ripple2.mat");
raw_data = cur_file.cur_filtered_data(floor((startt)*2.5):floor((endd)*2.5),:);
plot((1:numel(raw_data(:,1)))/2500,raw_data(:,50)) % let's show raw signal of random channel
xline(begin/2500);xline(einde/2500)
xlim([0.0075 0.18])
subplot(3,1,2)
power_ = sum(raw_data.^2,2);
plot((1:numel(power_))/2500,power_)
xline(begin/2500);xline(einde/2500)
xlim([0.0075 0.18])
subplot(3,1,3)
cur_correlation = raw_after_m{1};
plot((1:numel(cur_correlation(startt:endd)))/1000 + 0.0075,cur_correlation(startt:endd)) % shift to the right because of how correlation is calculated
hold on
plot((1:numel(cur_correlation(startt:endd)))/1000+0.0075,ones(size(1:numel(cur_correlation(startt:endd))))*threshold_after_m(1),'red')
ylim([0 2*threshold_after_m(1)])
xline(begin/2500);xline(einde/2500)
xlim([0.0075 0.18])


