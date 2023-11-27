function TemplateImport(folder, outname)

% Import data from phy template files (output of spikesorting with
% Kilosort/spyking circus + phy2).
% Make use of python helpers of the spike-master repository of CortexLab
% (https://github.com/cortex-lab/spikes)
%
% Inputs:
%   - folder: folder where spike sorting data are stored. 
%   - outname: name for saving results
%
% Mattia D'Andola,
% Apr 2020


%get spike infos from phy2 files
spk = loadKSdir(folder);

%Good clusters
goodN = spk.cids(spk.cgs == 2);
goodIdx = ismember(spk.clu, goodN);

%MUA clusters
muaN = find(spk.cgs == 1);
muaIdx = ismember(spk.clu, muaN);

%remember:
% - spikeTemplates are the templates initially assigned by the automatic
% sorting
% - spikeCluters are the ones after merging/splitting

% shanks = {}; %to be implemented once known how templateGUI processes multiple shanks

%save general infos
infos.data_path = folder; %path of the data
infos.data_type = spk.dtype; %type of data in bin file
infos.offset = 0; %data offset
infos.Fs = spk.sample_rate; %sample rate
infos.n_channels = spk.n_channels_dat; %number of channels in the electrode
infos.xcoords = spk.xcoords; %coordintes of the channels in the electrode
infos.ycoords = spk.ycoords;
infos.temps_shape = spk.temps; %template shape for each channel
infos.winv = spk.winv; %inverse of whitening mask

%separate MUA spikes from SU spikes
MUA = [spk.st(muaIdx), double(spk.clu(muaIdx)), nan(length(find(muaIdx==1)),1), spk.tempScalingAmps(muaIdx)];
Good = [spk.st(goodIdx), double(spk.clu(goodIdx)), nan(length(find(goodIdx==1)),1), spk.tempScalingAmps(goodIdx)];
    
%sort data by time and save
[~, Idx] = sort(MUA(:,1));
dest = [outname, '_MUA.mat'];
output = matfile(dest, 'Writable', true);
output.header = {'time', 'neuId', 'channel', 'tempScalingAmps'};
output.spikes = MUA(Idx, :);
output.infos = infos;

[~, Idx] = sort(Good(:,1));
dest = [outname, '_Good.mat'];
output = matfile(dest, 'Writable', true);
output.header = {'time', 'neuId', 'channel', 'tempScalingAmps'};
output.spikes = Good(Idx, :);
output.infos = infos;

end