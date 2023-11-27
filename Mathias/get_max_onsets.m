fields = fieldnames(event_main);
numbers = setdiff(7:34, [20, 32, 33, 34]);

cur_max = 0;

for i = 1:numel(fields)
    cur_onsets = event_main.(fields{i}).events.onsets;
    
    if length(cur_onsets) > cur_max
        cur_max = length(cur_onsets);
    end
end