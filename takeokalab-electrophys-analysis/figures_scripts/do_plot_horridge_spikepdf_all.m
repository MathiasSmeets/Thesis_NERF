clear all
recordings_database;

depth_binning = 10;

figure

[~, cdfs_M, ~, depthBins_M] = spike_pdf_depth(sort_masters, 'plot_results', 0, 'depth_binning', depth_binning, 'depth_offset', master_ch_limit);

%align results for different depth, sum and smooth 
N = length(depthBins_M);
depthlim = zeros(N,2);
ypl_M = [];
for i = 1:N
    depthX_M{i} = depthBins_M{i}(1:end-1)+mean(diff(depthBins_M{i}))/2;
    p1 = [depthX_M{i}];
    depthlim(i,:) = [min(p1) max(p1)];
end
[~,idmaxM] = max(depthlim(:,2));
ypl_M = zeros(1,size(depthX_M{idmaxM},2));

for i = 1:N
    ypl_M(i,1:size(depthX_M{i},2)) = smoothdata(sum(cdfs_M{i},2),'gaussian',10);
end

[~, cdfs_Y, ~, depthBins_Y] = spike_pdf_depth(sort_yokes, 'plot_results', 0, 'depth_binning', depth_binning,'depth_offset', yokes_ch_limit);
N = length(depthBins_Y);
depthlim = zeros(N,2);
ypl_Y = [];
for i = 1:N
    depthX_Y{i} = depthBins_Y{i}(1:end-1)+mean(diff(depthBins_Y{i}))/2;
    p1 = [depthX_Y{i}];
    depthlim(i,:) = [min(p1) max(p1)];
end
[~,idmaxY] = max(depthlim(:,2));
ypl_Y = zeros(1,size(depthX_Y{idmaxY},2));

for i = 1:N
    ypl_Y(i,1:size(depthX_Y{i},2)) = smoothdata(sum(cdfs_Y{i},2),'gaussian',10);
end

absmax = max([depthX_M{idmaxM}, depthX_Y{idmaxY}]);
subplot(1,2,1)
fill([nanmean(ypl_M)+nanstd(ypl_M) fliplr(nanmean(ypl_M)-nanstd(ypl_M))], [depthX_M{idmaxM} fliplr(depthX_M{idmaxM})], [.5 .5 .5], 'FaceAlpha', .5,'linestyle','none')
hold on
plot(nanmean(ypl_M), depthX_M{idmaxM}, '-k', 'linewidth', 1.5)
ylabel('depth (µm)')
xlabel('summed FR (spk/s)')
title('Master mice')
load(fullfile(sort_masters{idmaxM},'chanMap.mat'), 'ycoords');
gt = get(gca, 'ytick');
ylim([max(depthX_M{idmaxM})-absmax max(depthX_M{idmaxM})])
set(gca, 'ytick', gt, 'yticklabel', abs(gt-ycoords(master_ch_limit(idmaxM))))

subplot(1,2,2)
fill([nanmean(ypl_Y)+nanstd(ypl_Y) fliplr(nanmean(ypl_Y)-nanstd(ypl_Y))], [depthX_Y{idmaxY} fliplr(depthX_Y{idmaxY})], [.5 .5 .5], 'FaceAlpha', .5,'linestyle','none')
hold on
plot(nanmean(ypl_Y), depthX_Y{idmaxY}, '-k', 'linewidth', 1.5)
ylabel('depth (µm)')
xlabel('summed FR (spk/s)')
title('Yocked mice')
load(fullfile(sort_yokes{idmaxY},'chanMap.mat'), 'ycoords');
gt = get(gca, 'ytick');
ylim([max(depthX_Y{idmaxY})-absmax max(depthX_Y{idmaxY})])

set(gca, 'ytick', gt, 'yticklabel', abs(gt-ycoords(yokes_ch_limit(idmaxY))))


