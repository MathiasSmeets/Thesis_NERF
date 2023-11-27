clear;clc;close all;
recording_database_updated_20220718;

for i = setdiff(7:34, [20, 32, 33, 34])

    load(fullfile(sort_masters_horridge{i},'events.mat'),'events');
    save("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\eventM_" + i + ".mat", "events")
    load(fullfile(sort_yokes_horridge{i},'events.mat'),'events');
    save("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\eventY_" + i + ".mat", "events")


end