stimulus_data_m = load("X:\Mathias\switch_data\data_after_stimulus\after_stimulus_data_m_horridge.mat");
stimulus_data_m = stimulus_data_m.after_stimulus_data_m;

output_m = load("X:\Mathias\switch_data\neurons_of_interest\neurons_of_interest_horridge_m.mat");
output_m = output_m.output_m;

inhibited_m = load("X:\Mathias\switch_data\neurons_of_interest\inhibited_horridge_m.mat");
inhibited_m = inhibited_m.inhibited_m;

rest_data_m = load("X:\Mathias\switch_data\data_after_stimulus\after_data_m.mat");
rest_data_m = rest_data_m.after_data;
%% get template

amount_intervals = 60;
last_interval_data = zeros(1,size(stimulus_data_m,1));
template = cell(1,10);
for i = size(stimulus_data_m,1)
    % get last interval
    % get last interval
    counter = size(stimulus_data_m,2);
    last_interval_data(i) = counter;
    while isempty(stimulus_data_m{i,counter})
        counter = counter - 1;
        last_interval_data(i) = counter;
    end
    cur_avg = zeros(size(stimulus_data_m{i,last_interval_data(i)}));
    % get average of last 60 intervals
    for j = last_interval_data(i):-1:last_interval_data(i)-amount_intervals
        cur_avg = cur_avg + stimulus_data_m{i,j};
    end
    template{i} = cur_avg/amount_intervals;
end


%% check for replay in rest data
for i = size(stimulus_data_m)
    cur_correlation = zeros(1,(i));
    %for j = 1:last_interval_data(i)-1
    for j = 1:size(rest_data_m,2)
        cur_correlation(j) = sum(template{i}.*stimulus_data_m{i,j},'all') / (size(stimulus_data_m{i,j},1) * size(stimulus_data_m{i,j},2));
    end
end
