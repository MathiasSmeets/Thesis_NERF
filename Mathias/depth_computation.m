clear;close all; clc


recording_database_updated_20200805_JCfreshStart_1;
depth_np2 = cell(1,9);

for j = 1:3
    frM = load("X:\Mathias\switch_data\tables\frM_np2_"+j+".mat");frM = frM.frM;
    depths_M = zeros(size(frM,1),1);
for i = 1:size(frM,1)
    if isnan(frM.Id_spont(i))
        load(fullfile(sort_masters_horridge{frM.Recording(i)},'chanMap.mat'));
        clu = frM.Id_horr(i);
        idepth = get_cluster_depth(sort_masters_horridge{frM.Recording(i)}, clu);
    else
        load(fullfile(sort_masters_before{frM.Recording(i)},'chanMap.mat'));
        clu = frM.Id_spont(i);
        idepth = get_cluster_depth(sort_masters_before{frM.Recording(i)}, clu);
    end
    depths_M(i) = abs(idepth.depth' - ycoords(masters_depth(frM.Recording(i))));
end
depth_np2{j} = depths_M';
clearvars frM
end
