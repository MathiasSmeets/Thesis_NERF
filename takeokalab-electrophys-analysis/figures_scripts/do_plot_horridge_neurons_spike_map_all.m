clear all

%upload recordings list and measured depths
recordings_database;
nch = 385;
% tbin = [0:60:540; 60:60:600]';

%% compute spikes occurances time and depth for each neuron
allmasterneu = [];
allmasterdepth = [];
for i = 1:length(sort_masters)
%     %create time binning
    load(fullfile(raw_masters{i},'events.mat'),'events');
    evnames = {events.type};
    evtid = find(contains(evnames,'horridge'));
    
    % bin1 = first minute (adaptation) + bin2 = remaining part of stimulation (consolidation) + bin3 = remaining part of recording (retaining)
    file = dir(fullfile(raw_masters{i},'*.ap.bin'));
    maxtime = file.bytes/2/nch;
    tbin = [0 60; 60 events(evtid).offsets(end); events(evtid).offsets(end) maxtime]; 
    
    %compute spiking activity 
    [~, depth, neu_resp_ratio] = neurons_spike_map(sort_masters{i},'add_events_range',1,'events_name','horridge','plot_results',0,'time_binning',tbin);
    
    %adjust depths with real coordinates based on last channel inside
    load(fullfile(sort_masters{i},'chanMap.mat'), 'ycoords');
    depth = abs(depth - ycoords(master_ch_limit(i)));
    allmasterneu = [allmasterneu; neu_resp_ratio];
    allmasterdepth = [allmasterdepth depth];
end

% repeat all for Yokes
allyokneu = [];
allyokdepth = [];
for i = 1:length(sort_yokes)
%     %create time binning
    load(fullfile(raw_yokes{i},'events.mat'),'events');
    evnames = {events.type};
    evtid = find(contains(evnames,'horridge'));
    
    file = dir(fullfile(raw_yokes{i},'*.ap.bin'));
    maxtime = file.bytes/2/nch;
    % bin1 = first minute (adaptation) + bin2 = remaining part of stimulation (consolidation) + bin3 = remaining part of recording (retaining)
    tbin = [0 60; 60 events(evtid).offsets(end); events(evtid).offsets(end) maxtime]; 
    [~, depth, neu_resp_ratio] = neurons_spike_map(sort_yokes{i},'add_events_range',1,'events_name','horridge','plot_results',0,'time_binning',tbin);
    depth = abs(depth - ycoords(yokes_ch_limit(i)));

    allyokneu = [allyokneu; neu_resp_ratio];
    allyokdepth = [allyokdepth depth];
end

%average for depth
dbin = 0:50:1200;
dplot = 5:50:1200;
cont = 0;
neuMmean = zeros(length(dbin)-1,size(tbin,1));
neuYmean = zeros(length(dbin)-1,size(tbin,1));

for d = 2:length(dbin)
    cont = cont+1;
    idx = find(allmasterdepth>dbin(d-1) & allmasterdepth<=dbin(d));
    if ~isempty(idx)
        neuMmean(cont, :) = mean(allmasterneu(idx,:));
    end
    idx = find(allyokdepth>dbin(d-1) & allyokdepth<=dbin(d));
    if ~isempty(idx)
        neuYmean(cont, :) = mean(allyokneu(idx,:));
    end
end

%group for responses' time bin prevalency
%1: >= 75% during first minute
%2: >=75% during horridge, but more distributed
%3: inside horridge ~= outside horridge
%4: >=75% after horridge
group1M = [];
group2M = [];
group3M = [];
group4M = [];
for i = 1:size(allmasterneu, 1)
    if allmasterneu(i,1)>=.75
        group1M = [group1M i];
    elseif allmasterneu(i,1)+allmasterneu(i,2)>=.75
        group2M = [group2M i];
    elseif allmasterneu(i,1)+allmasterneu(i,2)>=.25 && allmasterneu(i,3)>=.25
        group3M = [group3M i];
    elseif allmasterneu(i,3)>=.75
        group4M = [group4M i];
    end
end

group1Y = [];
group2Y = [];
group3Y = [];
group4Y = [];
for i = 1:size(allyokneu, 1)
    if allmasterneu(i,1)>=.75
        group1Y = [group1Y i];
    elseif allyokneu(i,1)+allyokneu(i,2)>=.75
        group2Y = [group2Y i];
    elseif allyokneu(i,1)+allyokneu(i,2)>=.25 && allyokneu(i,3)>=.25
        group3Y = [group3Y i];
    elseif allyokneu(i,3)>=.75
        group4Y = [group3Y i];
    end
end


%% plot results

%scatter plots alla neurons
figure
cols = [1 0 0; 0 1 0; 0 0 1];
lab = {'adaptation', 'consolidation', 'retaining'};
for i = 1:3
    subplot(2,3,i)
    plot(allmasterneu(:,i), allmasterdepth, 'o', 'color', cols(i,:))
    xlabel('% of spiking during time period')
    ylabel('depth [µm]')
    title(lab{i})
end

for i = 4:6
    subplot(2,3,i)
    plot(allyokneu(:,i-3), allyokdepth, 'o', 'color', cols(i-3,:))
    xlabel('% of spiking during time period')
    ylabel('depth [µm]')
    title(lab{i-3})
end

% color scale plot averages
figure
subplot(1,2,1)
imagesc(1:size(tbin,1), dplot, neuMmean)
set(gca, 'Ydir', 'normal')

subplot(1,2,2)
imagesc(1:size(tbin,1), dplot, neuYmean)
set(gca, 'Ydir', 'normal')

% histograms groups
f = figure;
dmax = max([allmasterdepth allyokdepth]);

subplot(2,4,1)
h1 = histfit(allmasterdepth(group1M),10,'kernel');
h1(1).Horizontal = 'on';
Hy = h1(2).XData;
Hx = h1(2).YData;
h1(2).Visible = 'off';
hold on
plot(Hx, Hy, 'LineWidth',2)
hold off
set(gca,'ylim',[0 dmax],'ytick',0:200:dmax,'ydir','reverse')
ylabel('Depth [\mu m]')
xlabel('Neurons')
title('Master, Group 1')
text(max(h1(1).YData)-1,1200,sprintf('%d neurons',length(group1M)));

subplot(2,4,2)
h2 = histfit(allmasterdepth(group2M),10,'kernel');
h2(1).Horizontal = 'on';
Hy = h2(2).XData;
Hx = h2(2).YData;
h2(2).Visible = 'off';
hold on
plot(Hx, Hy, 'LineWidth',2)
hold off
set(gca,'ylim',[0 dmax],'ytick',0:200:dmax,'ydir','reverse')
ylabel('Depth [\mu m]')
xlabel('Neurons')
title('Master, Group 2')
text(max(h2(1).YData)-1,1200,sprintf('%d neurons',length(group2M)));


subplot(2,4,3)
h3 = histfit(allmasterdepth(group3M),10,'kernel');
h3(1).Horizontal = 'on';
Hy = h3(2).XData;
Hx = h3(2).YData;
h3(2).Visible = 'off';
hold on
plot(Hx, Hy, 'LineWidth',2)
hold off
set(gca,'ylim',[0 dmax],'ytick',0:200:dmax,'ydir','reverse')
ylabel('Depth [\mu m]')
xlabel('Neurons')
title('Master, Group 3')
text(max(h3(1).YData)-1,1200,sprintf('%d neurons',length(group3M)));


subplot(2,4,4)
h4 = histfit(allmasterdepth(group4M),10,'kernel');
h4(1).Horizontal = 'on';
Hy = h4(2).XData;
Hx = h4(2).YData;
h4(2).Visible = 'off';
hold on
plot(Hx, Hy, 'LineWidth',2)
hold off
set(gca,'ylim',[0 dmax],'ytick',0:200:dmax,'ydir','reverse')
ylabel('Depth [\mu m]')
xlabel('Neurons')
title('Master, Group 4')
text(max(h4(1).YData)-1,1200,sprintf('%d neurons',length(group4M)));


subplot(2,4,5)
h5 = histfit(allyokdepth(group1Y),10,'kernel');
h5(1).Horizontal = 'on';
Hy = h5(2).XData;
Hx = h5(2).YData;
h5(2).Visible = 'off';
hold on
plot(Hx, Hy, 'LineWidth',2)
hold off
set(gca,'ylim',[0 dmax],'ytick',0:200:dmax,'ydir','reverse')
ylabel('Depth [\mu m]')
xlabel('Neurons')
title('Yoked, Group 1')
text(max(h5(1).YData)-1,1200,sprintf('%d neurons',length(group1Y)));

subplot(2,4,6)
h6 = histfit(allyokdepth(group2Y),10,'kernel');
h6(1).Horizontal = 'on';
Hy = h6(2).XData;
Hx = h6(2).YData;
h6(2).Visible = 'off';
hold on
plot(Hx, Hy, 'LineWidth',2)
hold off
set(gca,'ylim',[0 dmax],'ytick',0:200:dmax,'ydir','reverse')
ylabel('Depth [\mu m]')
xlabel('Neurons')
title('Yoked, Group 2')
text(max(h6(1).YData)-1,1200,sprintf('%d neurons',length(group2Y)));

subplot(2,4,7)
h7 = histfit(allyokdepth(group3Y),10,'kernel');
h7(1).Horizontal = 'on';
Hy = h7(2).XData;
Hx = h7(2).YData;
h7(2).Visible = 'off';
hold on
plot(Hx, Hy, 'LineWidth',2)
hold off
set(gca,'ylim',[0 dmax],'ytick',0:200:dmax,'ydir','reverse')
ylabel('Depth [\mu m]')
xlabel('Neurons')
title('Yoked, Group 3')
text(max(h7(1).YData)-1,1200,sprintf('%d neurons',length(group3Y)));


subplot(2,4,8)
h8 = histfit(allyokdepth(group4Y),10,'kernel');
h8(1).Horizontal = 'on';
Hy = h8(2).XData;
Hx = h8(2).YData;
h8(2).Visible = 'off';
hold on
plot(Hx, Hy, 'LineWidth',2)
hold off
set(gca,'ylim',[0 dmax],'ytick',0:200:dmax,'ydir','reverse')
ylabel('Depth [\mu m]')
xlabel('Neurons')
title('Yoked, Group 4')
text(max(h8(1).YData)-1,1200,sprintf('%d neurons',length(group4Y)));

hs = [h1, h2, h3, h4, h5, h6, h7, h8];
%compare fittings
for i = 1:4
    xM{i} = get(hs(2,i), 'XData');
    yM{i} = get(hs(2,i), 'YData');
    
    xY{i} = get(hs(2,i+4), 'XData');
    yY{i} = get(hs(2,i+4), 'YData');
end

figure
colors = jet;
colors = colors(1:size(colors,1)/4:end,:);
subplot(1,2,1), hold on
for i = 1:4
    plot(xM{i}, yM{i}, 'color', colors(i,:))
end
xlabel('depth [\mu m]')
ylabel('# of neurons')
legend('group 1', 'group2', 'group3', 'group 4')
title('masters')

subplot(1,2,2), hold on
for i = 1:4
    plot(xY{i}, yY{i}, 'color', colors(i,:))
end
legend('group 1', 'group2', 'group3', 'group 4')
xlabel('depth [\mu m]')
ylabel('# of neurons')
legend('group 1', 'group2', 'group3', 'group 4')
title('yokes')


% maxC = max([max(h1(1).YData) max(h2(1).YData) max(h3(1).YData) max(h4(1).YData) max(h5(1).YData) max(h6(1).YData) max(h7(1).YData) max(h8(1).YData)]);
% fc = get(f, 'children');
% for i = 1:length(fc)
%     fc(i).XLim = [0 maxC];
% end
% hp2 = get(subplot(1,2,2), 'Position');
% cbar = colorbar('Position',  [hp2(1)+hp2(2)+0.01  hp2(2)  0.03  hp2(1)+hp2(2)*2.1]);
% cbar.Label.String = '% of spikes during phase';
% subplot(1,2,1)
% plot(allmasterneu, allmasterdepth, 'o')
% xlabel('% of spikes during learning')
% ylabel('depth')
% set(gca, 'ydir', 'reverse')
% title('master mice')
% 
% subplot(1,2,2)
% plot(allyokneu, allyokdepth, 'o')
% xlabel('% of spikes during learning')
% ylabel('depth')
% set(gca, 'ydir', 'reverse')
% title('yokes mice')

