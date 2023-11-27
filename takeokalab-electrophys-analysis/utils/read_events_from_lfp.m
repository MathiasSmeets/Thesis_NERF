[~, idx] = find(lfp>1000); %lfp 1 channel
fs = 30000;

onsets = [];
offsets = [];
savethefirst = [];
flag = 0;
for i = 1:length(idx)-1
if (idx(i+1) - idx(i))<1000
    if flag == 0
        flag = 1;
        onsets = [onsets, idx(i)];
    end 
else
    if flag == 0
        onsets = [onsets idx(i)];
    else
        flag = 0;
        offsets = [offsets, idx(i)];
    end
end
end
offsets = [offsets, idx(end)];

duration = [];
for i = 1:length(onsets)
    duration = [duration offsets(i)-onsets(i)];
end

%transform in seconds

events.onsets = onsets/fs;
events.offsets = offsets/fs;
events.duration = duration/fs;

