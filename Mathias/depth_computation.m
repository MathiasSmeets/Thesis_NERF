clear;close all; clc

%% compute depth of neurons
recording_database_updated_20200805_JCfreshStart_1;
depth_m = cell(1,9);
channel_m = cell(1,9);

for j = 1:9
    frY = load("X:\Mathias\switch_data\tables\frY_switched_"+j+".mat");frY = frY.frM;
    depths_M = zeros(size(frY,1),1);
    channels_M = zeros(size(depths_M));
    for i = 1:size(frY,1)
        if isnan(frY.Id_spont(i))
            load(fullfile(sort_yokes_horridge{frY.Recording(i)},'chanMap.mat'));
            clu = frY.Id_horr(i);
            idepth = get_cluster_depth(sort_yokes_horridge{frY.Recording(i)}, clu);
        else
            load(fullfile(sort_yokes_before{frY.Recording(i)},'chanMap.mat'));
            clu = frY.Id_spont(i);
            idepth = get_cluster_depth(sort_yokes_before{frY.Recording(i)}, clu);
        end
        depths_M(i) = abs(idepth.depth' - ycoords(masters_depth(frY.Recording(i))));
        channels_M(i) = idepth.channel;
    end
    depth_m{j} = depths_M';
    channel_m{j} = channels_M';
    clearvars frY
end