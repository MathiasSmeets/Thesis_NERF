function events = read_events_from_npy(folder, names, fs)

% Read TTL events saved as npy files from openephys (to be used when
% recording events with the digital input). To work it needs the npy-matlab
% library https://github.com/kwikteam/npy-matlab
%
% Input: folder cotaining the npy TTL files
%
% Mattia D'Andola, Oct 2020

if nargin<2, names = {}; end
if nargin<3, fs = 30000; end

    
ts = double(readNPY(fullfile(folder,'timestamps.npy')));
ch = double(readNPY(fullfile(folder,'channels.npy')));
fw = double(readNPY(fullfile(folder,'full_words.npy')));

nchs = unique(ch);
if isempty(names)
    for i = 1:length(nchs)
        names{i} = sprintf('TTL%d',i);
    end
else
    if ~iscell(names)
        names = {names};
    end
end
for i = 1:length(nchs)
    events.onsets = (ts(ch==i & fw==1))/fs';
    events.offsets = (ts(ch==i & fw==0))/fs';
    events.name = names{i};
end


    
    