% mouse mouse
% intervals 271-intervals_together0
clear;close all;clc

interval = 271; %301 2
mouse = 2;
intervals_to_plot = 80;

%% image 1: example

load("X:\Mathias\switch_data\clusters\assemblies_horridge_m.mat")
load("X:\Mathias\switch_data\data_after_stimulus\after_stimulus_data_m_horridge.mat")


figure
heatmap(after_stimulus_data_m{mouse,interval})
grid('off')
saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_1.png")
saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_1.svg")


%% image 2: neurons of interest

load("X:\Mathias\switch_data\neurons_of_interest\neurons_of_interest_horridge_m.mat")
load("X:\Mathias\switch_data\neurons_of_interest\inhibited_horridge_m.mat")
%neuron_counter = size(after_stimulus_data_m{1,1},1)+size(after_stimulus_data_m{2,1},1)+size(after_stimulus_data_m{3,1},1)+1;
neuron_counter = 1;
for i = 1:mouse-1
    neuron_counter = neuron_counter + size(after_stimulus_data_m{i,1},1);
end

noi = get_neurons_of_interest(after_stimulus_data_m{mouse,interval}, output_m, inhibited_m, neuron_counter);
figure
heatmap(after_stimulus_data_m{mouse,interval}(noi,:))
grid('off')
saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_2.png")
saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_2.svg")

%% image 3: take multiple bins together 

new_data = zeros(size(after_stimulus_data_m{mouse,interval}(noi,:),1),4);
for i = 1:size(new_data,1)
    for j = 1:size(new_data,2)
        new_data(i,j) = sum(after_stimulus_data_m{mouse,interval}(noi(i),11+(j-1)*15:10+j*15));
    end
end
figure
h=heatmap(new_data,'CellLabelColor','none');
grid('off')

saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_3.png")
saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_3.svg")

%% image mouse: multiple intervals combined

new_data = zeros(size(after_stimulus_data_m{mouse,interval}(noi,:),1),4*30);
counter = 1;
for i = 1:size(new_data,1)
    for j = interval:interval+29
        for k = 1:4
            new_data(i,counter) = sum(after_stimulus_data_m{mouse,j}(noi(i),11+(k-1)*15:10+k*15));
            counter = counter + 1;
        end
    end
    counter = 1;
end
figure
hm=heatmap(new_data(:,1:intervals_to_plot),'CellLabelColor','none');
% grid('off')
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(hm); % Undocumented
ax = S.Axes;    % Undocumented
clear('cleanup')
% Remove grids
hm.GridVisible = 'off';
% Place lines around selected columns and row
% Assumes columns and rows are 1 unit in size!
col = [0.5,4.5];    
row = [0.5, 10.5];
xline(ax, [col, col], 'k-', 'LineWidth', 1.5); % see footnotes [1,2]
%yline(ax, [row, row], 'k-', 'LineWidth', 1.5); % see footnotes [1,2]

saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_4.png")
saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_4.svg")


%% image 5: clusters ica

load("X:\Mathias\switch_data\clusters\assemblies_horridge_m.mat")
index = ceil(interval/30);
assembly = total_assemblies{mouse,index}{1};

assembly_data = new_data(assembly,:);
figure
heatmap(assembly_data(:,1:intervals_to_plot))
grid('off')

saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_5.png")
saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_5.svg")

%% image 5b: activity

load("X:\Mathias\switch_data\clusters\activity_horridge_m.mat")
activity = total_activity{mouse,index}(1:intervals_to_plot);
figure
plot(abs(activity),'LineWidth',2)
ylabel("Activity")

[pks, locs] = findpeaks(abs(activity), "MinPeakHeight", 0.4*max(abs(activity)));
hold on
scatter(locs, pks, "*", "LineWidth",2)

saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_5b.png")
saveas(gcf, "C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\images\overview_figure\fig_5b.svg")

%% image 6

for jj = 1:length(locs)
    raw_data_index = (j-1)*30 + ceil((locs(jj))/ceil(60/4));
    position_in_data = mod(locs(jj)-1,ceil(60/15))*15+1;
    cur_raw_data = stimulus_data_m{mouse,raw_data_index};
    cur_assembly_data = cur_raw_data(template_cluster{mouse}, position_in_data:min(position_in_data+15-1, 70));
    cur_template(:,1:size(cur_assembly_data,2)) = cur_template(:,1:size(cur_assembly_data,2)) + double(cur_assembly_data);
    counter = counter + 1;
end







