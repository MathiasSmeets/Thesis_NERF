function preprocess_openephys_recording(dfile, sfile, rec_chs, eve_chs, eve_names,sep_chs, cut, filter_params)


% Preprocess openephys recordings to create:
% - a data file with no header for spike sorting
% - an event file for event analog channels (if any)
%
% Inputs: 
% - dfile = complete path for output file (with no extension!!)
% - sfile = complete path to structure json file (usually
% 'structure.obein'). If sep_chs = 1, sfile = folder that contains the
% separate channel recordings
% - rec_chs = index of recording channels
% - eve_chs = index of event channels or path to folder containing npy
% files for digital TTL
% - sep_chs = binary. set it 1 if recording has separated channel files and
% no json. In that case, put the folder as sfile
% - cut = parts to remove from the recording [start sec , end sec];
%       use '0' to indicate boundaries -- examples:
%       [0, 100] -> upload only the first 100 s
%       [30, 100] -> upload from t=30s to t=100s
%       [100, 0] -> upload from t=100 to the end of the recording
% - filter_params = if specified, data are filtered following the specified
% parameters (see function 'data_filter' for parameters specifications)
%
% Mattia D'Andola, Jan 2020

% use this channel map if recording with intan 32 ch headstage and using
% CN 32 channel 1b electrode (and its cooresponding channelmap for sorting)

CHMAP = [8, 11, 13, 2, 9, 12, 14, 1, 16, 15, 3, 4, 5, 6, 10, 7, 25, 24, 22, 21, 20, 19, 31, 32, 17, 18, 30, 29, 28, 27, 23, 26];

if nargin < 4
    eve_chs = [];
end
if nargin < 5
    sep_chs = 0;
end
if nargin < 6
    cut = [0, 0];
end
if nargin < 7
    filter_params = [];
end

% upload data
disp('Uploading data...')
if sep_chs
    for i = rec_chs
        namefile = [sfile, '/101_CH', num2str(i), '.continuous'];
        [Dtemp, ~, info] = load_open_ephys_data_faster(namefile);
        D.Data(i,:) = single(Dtemp);
    end
    D.Header.sample_rate = info.header.sampleRate;
else
    D = load_open_ephys_binary(sfile, 'continuous',1,'mmap');
end

%select data of electrophys recordings, and save .dat file
if ~isempty(rec_chs)
    %calculate start and end bin of data save
    if cut(1) == 0
        startbin = 1;
    else
        startbin = cut(1)*D.Header.sample_rate;
    end
    if cut(2) == 0
        endbin = size(D.Data.Data(1).mapped,2);
    else
        endbin = cut(2)*D.Header.sample_rate;
    end
    if isstruct(D)
        for i = rec_chs
            lfp(i,:) = D.Data.Data(1).mapped(i,startbin:endbin);
        end
    else
        lfp = D.Data(rec_chs, startbin:endbin);
    end
    %reorder channels
    lfp = lfp(CHMAP, :);
    if ~isempty(filter_params)
        disp('Filtering data. This process could take a while...')
        lfp = data_filter(lfp,filter_params);
    end
    lfpnew = reshape(lfp,[1, size(lfp,1)*size(lfp,2)]);

    disp('Writing new data file for spike sorting...')
    fid = fopen([dfile '.dat'], 'w');
    fwrite(fid, lfpnew, 'int16');
    fclose(fid);
end
    
%select data of analog events
if ~isempty(eve_chs)
    disp('Writing an events file for analysis...')
    if isnumeric(eve_chs)
        for e = length(eve_chs)
            if sep_chs 
                %to be implemented
                disp('implement it!')
            else
                stimdata = D.Data.Data(1).mapped(eve_chs(e),:);
            end

            %clean data
            stimdata = stimdata - median(stimdata); %subtract baseline
            [~, idx] = find(stimdata>max(stimdata)*0.005); % find points in which the signal is higher than the baseline

            %create a binary vector for stim on/off
            trsign = zeros(1, size(stimdata,2));
            trsign(idx) = 1;

            %find transitions
            trpoints = diff(trsign);
            [~, onsets] = find(trpoints == 1);
            [~, offsets] = find(trpoints == -1);
        % 
        %     % cut the last one to be sure t have the same size
        %     onsets = onsets(1:end-1);

            % transform to seconds
            onsets = (onsets*(1/D.Header.sample_rate))';
            offsets = (offsets*(1/D.Header.sample_rate))';

            events(e).stim = trsign;
            events(e).onsets = onsets;
            events(e).offsets = offsets;
            events(e).type = eve_names{e};
        end
    elseif ischar(eve_chs)
        events = read_events_from_npy(eve_chs);
    end
    save([dfile '_events.mat'], 'events', '-v7.3')
end
