function KwikImport(folder, filename, fs)

% Import data from kwik/kwx files (output of spikesorting with
% klusta2+phy2).
%
% Mattia D'Andola,
% Jan 2020


[~, outname, ~] = fileparts(filename);
path = fullfile(folder, filename);
maskPath = fullfile(folder, [outname, '.kwx']);

info = h5info(path, '/channel_groups');
shanks = {info.Groups.Name};
MUA = [];
Good = [];
maxClus = 0;

%iterate on electrode's shanks
for i = 1:length(shanks)
    info = h5info(path, [shanks{i}, '/clusters/main']);
    clusters = {info.Groups.Name};

    chMask = squeeze(h5read(maskPath, [shanks{i}, '/features_masks'],...
        [2,1,1],[1,Inf,Inf]));
    shankCh = double(h5readatt(path, shanks{1},'channel_order') + 1);
    
    valid = cellfun(@(x) uint8(h5readatt(path, x, 'cluster_group')), clusters);
    clusters = cellfun(@(x) strsplit(x, '/'), clusters, 'UniformOutput', false);
    clusters = cellfun(@(x) str2num(x{end}), clusters);
    clMUA = sort(clusters(valid == 1));
    clGood = sort(clusters(valid == 2));
    
    clusters = hdf5read(path, [shanks{i}, '/spikes/clusters/main']);
    channels = zeros(size(clusters));
    for j = unique(clusters)'
        clusterIdx = clusters == j;
        nSpk = sum(clusterIdx);
        cl = cellfun(@(x) (x>0), num2cell(chMask(1:3:end, clusterIdx), 1),...
            'UniformOutput', false);
        cl = cell2mat(cl);
        channels(clusterIdx) = sum(sum(cl,2).*shankCh/sum(cl(:)));  
    end
    
    times = hdf5read(path, [shanks{i}, '/spikes/time_samples']);
    times = double(times) / fs;
    muaIdx = ismember(clusters, clMUA);
    goodIdx = ismember(clusters, clGood);
    
    tmpNeus = clusters + maxClus;
    
    maxClus = maxClus + max(tmpNeus) + 1;
    
    MUA = [MUA; [times(muaIdx), double(tmpNeus(muaIdx)), channels(muaIdx), double(clusters(muaIdx))]];
    Good = [Good; [times(goodIdx), double(tmpNeus(goodIdx)), channels(goodIdx), double(clusters(goodIdx))]];
    
end

[~, Idx] = sort(MUA(:,1));

dest = [outname, '_MUA.mat'];
output = matfile(dest, 'Writable', true);
output.header = {'time', 'neuId', 'channel', 'clusterId'};
output.spikes = MUA(Idx, :);

[~, Idx] = sort(Good(:,1));

dest = [outname, '_Good.mat'];
output = matfile(dest, 'Writable', true);
output.header = {'time', 'neuId', 'channel', 'clusterId'};
output.spikes = Good(Idx, :);

end