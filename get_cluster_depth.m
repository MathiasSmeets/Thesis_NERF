function depths = get_cluster_depth(spkpath, clusters)

%calculate the depth of a specific cluster

[cldata, clheader, xlraw] = tsvread(fullfile(spkpath,'cluster_info.tsv'));
idx = []; 
depths = table();
depths.cluster = clusters';
for i = 1:length(clusters)
    idx(i) = find(cldata(:,1) == clusters(i));
end
depths.depth = cldata(idx,7);
depths.channel = cldata(idx,6);

