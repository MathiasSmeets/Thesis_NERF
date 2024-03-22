clear;clc;close all;

%% get data

if ispc
    volume_base = '\\nerffs13';
    volume_base2 = '\\nerffs17';
elseif ismac
    volume_base = '/volumes/';
else  % cluster
    volume_base = '/mnt/takeokalab/';
    volume_base2 = '/mnt/takeokalab/';
end

path_to_code = "\takeokalabwip2023\Mathias\switch_data\tables";
write_path = fullfile(volume_base2, "/takeokalabwip2023/Mathias/switch_data/sos_data");


%%
sos_results_m = cell(3,2);
sos_results_y = cell(3,2);
counter = 1;
for i = 1:3
    cur_frM = load(fullfile(volume_base2, path_to_code, "frM_np2_" + i));
    cur_frM = cur_frM.frM;

    % remove rows that contain zero in the recording column
    cur_frM(~cur_frM.Recording,:) = [];

    baseline_data_m = cur_frM.Fr_spont;
    cur_average_m = mean(cell2mat(baseline_data_m),2);
    cur_std_m = std(cell2mat(baseline_data_m), [], 2);
    
    sos_results_m{counter, 1} = cur_average_m;
    sos_results_m{counter, 2} = cur_std_m;

    clearvars cur_frM

    % cur_frY = load(fullfile(volume_base2, path_to_code, "frY_np2_" + i));
    % cur_frY = cur_frY.frM;
    % 
    % % remove rows that contain zero in the recording column
    % cur_frY(~cur_frY.Recording,:) = [];
    % 
    % baseline_data_y = cur_frY.Fr_spont;
    % cur_average_y = mean(cell2mat(baseline_data_y),2);
    % cur_std_y = std(cell2mat(baseline_data_y), [], 2);
    % 
    % sos_results_y{counter, 1} = cur_average_y;
    % sos_results_y{counter, 2} = cur_std_y;
    % 
    % clearvars cur_frY
    % 
    counter = counter + 1;
end

save(fullfile(write_path, "sos_results_np2"), "sos_results_m")
%save(fullfile(write_path, "sos_results_y"), "sos_results_y")