% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Creation of frY and frM table for analysis
%
% 29/04/22
% Simon Lavaud, Modified by Charlotte Bichara
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear;clc;close all;
path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/";
path_to_jeremy = "/mnt/takeokalab/takeokalabwip2023/Jeremy/";
path_1 = genpath(path_to_code);
path_2 = genpath(path_to_jeremy);
addpath(path_1)
addpath(path_2)
% stimulation: Y:\past lab members\Mattia\Spinal_Cord\Sorted\20210511\horridge

for mathias_param = setdiff(7:34, [20, 32, 33, 34])
%for mathias_param = 10

% - - - Help on how to use it:
%
% You will probably already have a frM frY table from preivous recording
% The best is to load your previous table, rename it as frM_old frY_old.
% Then you can run for a specific recording the frM frY algorithm.
% It will return a frM frY only for one recording.
% You can copy and past this new frM frY in the old version.
% Then delete the new frM and frY, and rename the frM_old and frY_old
% as frM and frY. You can now save it in the right folder.



% Select 1 for bad masters/yokes and 0 for good ones
% BAD = 0;
% Select the database of interest
recording_database_updated_20220718;

% Initialize the variable
nch = 385;
Fs = 30000;
analysis_type = 1;

% Column of the table for analysis
Recording = [];
Id_spont = [];
Fr_spont = [];
Id_horr = [];
Fr_early = [];
Fr_middle = [];
Fr_late = [];
class = [];


%%

% % % % % % % % % % % % % % % % % % % % % % % %
% - - -
% - - - - - frM
% - - -
% % % % % % % % % % % % % % % % % % % % % % % %

clearvars frM

% Initialize inner variable
cont = 0;

% Iteration on the recording
% If only one recording to add : i = recording number
% If all the recording: [7 8 10 11 12 13 14 15 16 17 18 19 20 21 22]

% i = number of the recording
for i = mathias_param
    QQ=i;
    
    % Display the recording that is analysed
    disp(sprintf('Master %d',i))
    
    spk = loadKSdir(sort_masters_horridge{i});
    clusters = spk.cids(spk.cgs == 2);
    
    %crete time binning
    load(fullfile(sort_masters_horridge{i},'events.mat'),'events');
    save("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\eventM_" + mathias_param + ".mat", "events")
    load(fullfile(sort_yokes_horridge{i},'events.mat'),'events');
    save("C:\Users\Mathi\OneDrive\Documenten\Master_3\Thesis\code\Mathias\data\eventY_" + mathias_param + ".mat", "events")

    evnames = {events.type};
    evtid = find(contains(evnames,'horridge'));
    if isempty(evtid)==1
        evtid = find(contains(evnames,'opto'));
    end
    
    %Load first and last event for Eearly and Late
    %For it: look at the event file.

    if i==8 % WT
        %End of master horridge number 8 at 32min
        lim_start=1; %length(find(events.onsets<(22*60)));
        lim_event=156; %length(find(events.onsets<(32*60) & events.onsets>(22*60)));
    elseif i==9 % Tlx3
        %Start horridge 22 End of master horridge number 9 at 32min
        lim_start=259;
        lim_event=1036;
    elseif i==10 %09/11/21 // Tlx3
        lim_start=1;
        lim_event=1394;
        %         elseif i==12 %22/11/21 // Ptf1a
        %             lim_start=length(find(events.onsets<(20*60)));
        %             lim_event=length(find(events.onsets<(60*30) & events.onsets>(20*60)));
    elseif i==11 %17/11/21 // Ptf1a
        lim_start=1; %length(find(events.onsets<(15*60)));
        lim_event=2035; %length(find(events.onsets<(60*23) & events.onsets>(15*60)));
    elseif i==12 % Ptf1a
        lim_start=21; %OK checked CB;
        lim_event=578; %length(find(events.onsets<(60*23) & events.onsets>(15*60)));
    elseif i==13 %06/12/21 // Ptf1a
        lim_start=length(find(events.onsets<(15*60)));
        lim_event=length(find(events.onsets<(60*23) & events.onsets>(15*60)));
    elseif i==14 %04/01/22 // Ptf1a
        lim_start=find(events.onsets>(15*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=length(find(events.onsets<(60*23) & events.onsets>(15*60)));
    elseif i==15 %11/01/22 // Ptf1a
        lim_start=find(events.onsets>(15*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=length(find(events.onsets<(60*23) & events.onsets>(15*60)));
    elseif i==16 %18/01/22 // Ptf1a
        lim_start=find(events.onsets>(14*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=length(find(events.onsets<(60*23) & events.onsets>(15*60)));
    elseif i==17 %25/01/22 // Ptf1a
        lim_start=133; % OK checked CB;;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=1345; %length(find(events.onsets<(60*25) & events.onsets>(18*60)));
    elseif i==18 %07/02/22 // Ptf1a
        lim_start=3;%find(events.onsets>(17*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=1228;%length(find(events.onsets<(60*25) & events.onsets>(17*60)));
    elseif i==19 %14/02/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=286; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==20 %21/02/22 // Ptf1a
        lim_start=4; %find(events.onsets>(15*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=885; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==21 %10/03/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=569; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==22 %14/03/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=293; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==23 %24/03/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        lim_event=1787; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==24 %19/04/22 // Ptf1a
        lim_start=2; %find(events.onsets>(15*60),1)+1;
        lim_event=788; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==25 %28/02/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=239; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==26 %25/05/22 // Ptf1a
        lim_start=2; %find(events.onsets>(15*60),1)+1;
        %lim_start=length(find(events.onsets<(16*60) & events.onsets>(10*60)));
        lim_event=3372; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==27 %22/06/22 // Ptf1a
        lim_start=1; 
        lim_event=610;
    elseif i==28 %25/07/22 // Ptf1a
        lim_start=1; 
        lim_event=324;
    elseif i==29 %19/08/22 // Ptf1a
        lim_start=1; 
        lim_event=336;
    elseif i==30 %21/09/22 // Ptf1a
        lim_start=2; 
        lim_event=1070;
    elseif i==31 %23/09/22 // Ptf1a
        lim_start=1; 
        lim_event=128;
%     elseif i==32 %05/10/22 // Ptf1a  /Removed /!\ done Using "events old". If you use the corected eventfile ("events"), change the lim.
%         lim_start=131; 
%         lim_event=354;   
    elseif i==32 %11/10/22 // Ptf1a 
        lim_start=1; 
        lim_event=103;         
    elseif i==33 %15/11/22 // Ptf1a 
        lim_start=1; 
        lim_event=183; %corrected
    elseif i==34 %16/11/22 // Ptf1a   (Ran as a Yoked the Copy pasted)
        lim_start=3;
        lim_event=419;
    

        
        % Initialization parameter in case of problem
    else
        lim_start=2;
        lim_event=length(events(evtid).onsets);
    end
    
    if lim_start<1
        lim_start=2;
    end
    
    % Computation of the timeline Early, mid, late
    % Mid not use anymore
    clearvars dstim tdist
    a=1;
    %Frequency
    for j = (lim_start+1):lim_event
        tdist(j) = events(evtid).onsets(j) - events(evtid).onsets(j-1);
        %a=a+1;
    end
    
    dstim = smoothdata(tdist);
    idmax = find(dstim>1,1);
    if isempty(idmax)
        idmax=lim_start+10;
    end
    
    if i==14 | i==16
        filter_non_onset=find(events.onsets>800 & events.onsets<1000,1);
    else
        filter_non_onset=1;
    end
    
    % Tbin
    % 0 to 600 - spont
    % 600 to time early (early phase)
    % time early to time mid (mid phase)
    % time mid to end horridge (late phase)
    % end horridge - end+600 (after phase)
    
    if ~isempty(idmax)
        tbin = [events(evtid).onsets(lim_start)-600, events(evtid).onsets(lim_start)-1;... %Spont
            events(evtid).onsets(lim_start), events(evtid).onsets(idmax);... %Early
            events(evtid).onsets(idmax), events(evtid).onsets(lim_event);... %Mid
            events(evtid).onsets(lim_event), events(evtid).onsets(lim_start)+600]; %Late
    else
        tbin = [events(evtid).onsets(1)-600, events(evtid).onsets(1)-1;...
            events(evtid).onsets(1), 60;...
            60, events(evtid).onsets(end);...
            events(evtid).onsets(end), events(evtid).onsets(1)+600];
    end
    
    % Create columm :
    % Recording number
    % Id horridge / Id spont / Id after (same on long recording)
    L = 0.0001;
    step = 0.1;
    for j = 1:length(clusters)
        cont = cont+1;
        disp(sprintf('master %d, cluster long %d',i,j))
        Recording(cont) = i;
        Id_horr(cont) = clusters(j);
        Id_spont(cont) = clusters(j);
        Id_after(cont) = clusters(j);
        
        % Create the FR column for each phases depending on Tbin
        for k = 1:size(tbin,1)
            win = [tbin(k,1) tbin(k,1)+L];
            frmov = [];
            edges = tbin(k,1):L:tbin(k,2);
            frmov = histcounts(spk.st(spk.clu == clusters(j)),edges);
            switch k
                case 1
                    Fr_spont{cont} = frmov;
                case 2
                    Fr_early{cont} = frmov;
                case 3
                    Fr_middle{cont} = frmov;
                case 4
                    Fr_late{cont} = frmov;
            end
        end
        
        k=4;
        win = [1500 1500+L];
        frmov = [];
        edges = win(1):L:win(2);
        frmov = histcounts(spk.st(spk.clu == clusters(j)),edges);
        Fr_after{cont} = frmov;
        
        Fr_average_after(cont)=mean(frmov);
        
        Fr_average_spont(cont)=mean(Fr_spont{cont});
        Tbin_tot{cont}=tbin;
    end

end


% Create the frM table

frM = table();
frM.Recording = Recording';
frM.Tbin = Tbin_tot';
frM.Id_spont = Id_spont';
frM.Fr_spont = Fr_spont';
frM.Fr_average_spont = Fr_average_spont';
frM.Id_horr = Id_horr';
frM.Fr_early = Fr_early';
frM.Fr_middle = Fr_middle';
frM.Fr_late = Fr_late';
frM.Id_after = Id_after';
frM.Fr_after = Fr_after';
frM.Fr_average_after = Fr_average_after';




%%

% % % % % % % % % % % % % % % % % % % % % % % %
% - - -
% - - - - - frY
% - - -
% % % % % % % % % % % % % % % % % % % % % % % %

% Work exactly as the Master code

RecordingY = [];
Id_spontY = [];
Fr_spontY = [];
Id_horrY = [];
Fr_earlyY = [];
Fr_middleY = [];
Fr_lateY = [];
classY = [];

cont = 0;
for i = mathias_param % [7 8 9 10 11 12 13 14 15 16 17 18] %1:length(sort_yokes_horridge) 7 8 9 10 11 12 13
    QQ=i;
    
    disp(sprintf('Yoke %d',i))
    
    spk = loadKSdir(sort_yokes_horridge{i});
    clusters = spk.cids(spk.cgs == 2);
    
    %crete time binning
    load(fullfile(sort_yokes_horridge{i},'events.mat'),'events');
    evnames = {events.type};
    evtid = find(contains(evnames,'horridge'));
    if isempty(evtid)==1
        evtid = find(contains(evnames,'opto'));
    end
    
    if i==7 %31/08/21       
        lim_start=1; %1; horridge fail first time %length(find(events.onsets<(17*60)));
        lim_event=1603; %length(find(events.onsets<(60*25) & events.onsets>(60*17)));
    elseif i==8
        %End of yoked horridge number 8 at 27min / start at 17
        lim_start=97; %1; horridge fail first time %length(find(events.onsets<(17*60)));
        lim_event=252; %length(find(events.onsets<(60*25) & events.onsets>(60*17)));
    elseif i==9 %16/11/21
        lim_start=1;
        lim_event=1382;
    elseif i==10 %22/11/21
        lim_start=1; %length(find(events.onsets<(15*60)));
        lim_event=1826; %length(find(events.onsets<(60*25) & events.onsets>(60*15)));
    elseif i==11 %02/12/21
        lim_start=1; %11;%length(find(events.onsets<(15.1*60)));
        lim_event=558;%length(find(events.onsets<(60*25) & events.onsets>(60*15.1)));
    elseif i==12 %07/12/21
        lim_start=length(find(events.onsets<(15*60)));
        lim_event=length(find(events.onsets<(60*23) & events.onsets>(60*15)));
    elseif i==13 %05/01/22
        lim_start=find(events.onsets<(18*60) & events.onsets>(16*60),1)+1;
        lim_event=204; %197; %length(find(events.onsets<(60*25) & events.onsets>(60*17)));
    elseif i==14 %14/01/22
        test=events.onsets;
        lim_start=find(test>(15*60) & test<(16*60),1)+1;
        lim_event=length(find(test<(60*23) & test>(60*15)));
    elseif i==15 %21/01/22
        test=events.onsets;
        lim_start=find(test>(14*60) & test<(16*60),1)+1;
        lim_event=length(find(test<(60*23) & test>(60*15)));
    elseif i==16 %01/02/22
        %test=events.onsets;
        lim_start=1;%find(test>(14*60) & test<(16*60),1)+1;
        lim_event=1213;%length(find(test<(60*23) & test>(60*15)));
    elseif i==17 %08/02/22
        test=events.onsets;
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=1224;%length(find(test<(60*25) & test>(60*17)));
    elseif i==18 %17/02/22
        test=events.onsets;
        lim_start=3; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=287;%length(find(test<(60*25) & test>(60*17)));
    elseif i==19 %01/03/22
        test=events.onsets;
        lim_start=76; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=957; %length(find(test<(60*25) & test>(60*17)));
    elseif i==20 %11/03/22
        test=events.onsets;
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=569; %length(find(test<(60*25) & test>(60*17)));
    elseif i==21 %23/03/22
        test=events.onsets;
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=293; %length(find(test<(60*25) & test>(60*17)));
    elseif i==22 %13/04/22
        test=events.onsets;
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=1786; %length(find(test<(60*25) & test>(60*17)));
    elseif i==23 %27/04/22
        test=events.onsets;
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=787; %length(find(test<(60*25) & test>(60*17)));
    elseif i==24 %29/04/22
        test=events.onsets;
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=239; %length(find(test<(60*25) & test>(60*17)));     
    elseif i==25 %07/06/22
        test=events.onsets;
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=2648; %length(find(test<(60*25) & test>(60*17)));     
    elseif i==26 %05/07/22 // Ptf1a
        lim_start=1;
        lim_event=610;
    elseif i==27 %26/07/22 // Ptf1a
        lim_start=1;
        lim_event=568;
    elseif i==28 %05/09/22 // Ptf1a
        lim_start=1;
        lim_event=279;
    elseif i==29 %22/09/22 // Ptf1a
        lim_start=1;
        lim_event=1035;
    elseif i==30 %27/09/22 // Ptf1a   
        lim_start=1;
        lim_event=127;
    elseif i==31 %06/10/22 // Ptf1a   
        lim_start=2;
        lim_event=221;
    elseif i==32 %18/10/22 // Ptf1a   
        lim_start=1;
        lim_event=183; %corrected
    elseif i==33 %07/12/22 // Ptf1a
        lim_start=107;
        lim_event=287;
    elseif i==34 %19/01/23 // Ptf1a
        lim_start=0;
        lim_event=416;    
        
    else
        lim_start=2;
        lim_event=length(events(evtid).onsets);
    end
    
    if lim_start<1
        lim_start=2
    end
    
    clearvars dstim tdist a
    a=1;
    for j = (lim_start+1):lim_event
        tdist(a) = events(evtid).onsets(j) - events(evtid).onsets(j-1);
        a=a+1;
    end
    
    if i==13
        dstim = smoothdata(tdist,'gaussian', 2);
    else
        dstim = smoothdata(tdist);
    end
    
    idmax = find(dstim>1,1);
    if isempty(idmax)
        idmax=lim_start+10;
    end
    
    %if i==8
    %filter_non_onset=find(events.onsets>800 & events.onsets<1000,1);
    if i==13
        filter_non_onset=find(events.onsets>800 & events.onsets<1400,1);
    else
        filter_non_onset=1;
    end
    
    if i==15
        idmax=779;
    end
    
    
    
    if ~isempty(idmax)
        if i == 33 % Faux depart noise before early
            tbin = [300, 899;... %Spont
            events(evtid).onsets(lim_start), events(evtid).onsets(lim_start+1);... %Early
            events(evtid).onsets(lim_start+1), events(evtid).onsets(lim_event);... %Mid
            events(evtid).onsets(lim_event), events(evtid).onsets(lim_start)+600]; %Late

        else
        tbin = [events(evtid).onsets(lim_start)-600, events(evtid).onsets(lim_start)-1;... %Spont
            events(evtid).onsets(lim_start), events(evtid).onsets(idmax);... %Early
            events(evtid).onsets(idmax), events(evtid).onsets(lim_event);... %Mid
            events(evtid).onsets(lim_event), events(evtid).onsets(lim_start)+600]; %Late
        end

    else
        tbin = [events(evtid).onsets(1)-600, events(evtid).onsets(1)-1;...
            events(evtid).onsets(1), 60;...
            60, events(evtid).onsets(end);...
            events(evtid).onsets(end), events(evtid).onsets(1)+600];
    end
    
    L = 0.0001; % timebin in seconds
    step = 0.1;
    for j = 1:length(clusters)
        disp(sprintf('master %d, yoke long %d',i,j))
        cont = cont+1;
        RecordingY(cont) = i;
        Id_horrY(cont) = clusters(j);
        Id_spontY(cont) = clusters(j);
        Id_afterY(cont) = clusters(j);
        
        for k = 1:size(tbin,1)
            win = [tbin(k,1) tbin(k,1)+L];
            frmov = [];
            edges = tbin(k,1):L:tbin(k,2);
            frmov = histcounts(spk.st(spk.clu == clusters(j)),edges);
            switch k
                case 1
                    Fr_spontY{cont} = frmov;
                case 2
                    Fr_earlyY{cont} = frmov;
                case 3
                    Fr_middleY{cont} = frmov;
                case 4
                    Fr_lateY{cont} = frmov;
            end
        end
        k=4;
        win = [1500 1500+L];
        frmov = [];
        edges = win(1):L:win(2);
        frmov = histcounts(spk.st(spk.clu == clusters(j)),edges);
        Fr_afterY{cont} = frmov;
        
        Fr_average_afterY(cont)=mean(frmov);
        Fr_average_spontY(cont)=mean(Fr_spontY{cont});
        Tbin_totY{cont}=tbin;
        
    end
end


frY = table();
frY.Recording = RecordingY';
frY.Tbin = Tbin_totY';
frY.Fr_average_spont = Fr_average_spontY';
frY.Id_spont = Id_spontY';
frY.Fr_spont = Fr_spontY';
frY.Id_horr = Id_horrY';
frY.Fr_early = Fr_earlyY';
frY.Fr_middle = Fr_middleY';
frY.Fr_late = Fr_lateY';
frY.Id_after = Id_afterY';
frY.Fr_after = Fr_afterY';
frY.Fr_average_after = Fr_average_afterY';

%%
% % % % % % % % % % % % % % % % % % % % % % % %
% - - -
% - - - - - Up-regulation / Down-regulation frM frY
% - - -
% - - -
% - - -
% - - - - - 0=non participatiig; 1:increased; 2=decreased; 3=disappeared
% - - -
% - - -
% % % % % % % % % % % % % % % % % % % % % % % %

% - - - - - - - - - -
% - - - - - Initialize all the variable
% - - - - - - - - - -
clearvars compa_after zdata_early zdata_mid zdata_late zdata
clearvars evolution class2 evolution_cdf_pre_post evolution_fi evolution_cdf_pre_post_fi
clearvars compa_after zdata_early zdata_early_up zdata_early_down
clearvars zdata_mid zdata_mid_up zdata_mid_down zdata_late zdata
clearvars zdata_late zdata_late_up zdata_late_down
clearvars zdata zdata2 zdata3 zgroup_early zgroup_mid zgroup_late
clearvars mod_early_up_zscore mod_early_down_zscore
clearvars mod_mid_up_zscore mod_mid_down_zscore
clearvars mod_late_up_zscore mod_late_down_zscore
clearvars zdata_EandMmerge_up zdata_EandMmerge_down zdata_EandMmerge
clearvars mod_EandMmerge_up_zscore mod_EandMmerge_down_zscore zgroup_EandMmerge
clearvars EandMmerge


% - - - Parameter for up/down regulation
% Threshold of z-score detection
threshold_zscore=2;
% Threshold for p-value comparison (old methods)
pvalue_threshold=0.001;
% Smoothing factor of the Fr (0.001 is no smoothing)
smoothing_factor=0.001;
% Threshold for persistent/intermitent modulation (0 = 0%, 1=100% of Fr modulated)
threshold_per_zscore=0.5;
threshold_weak=0.2;

%%

% - - - - - - - -
% - - - Start calculating FRM
% - - - - - - -

clc
clearvars zgroup_EandMmerge zdata_EandMmerge zdata_EandMmerge_up zdata_EandMmerge_down
clearvars zgroup_late zdata_late zdata_late_up zdata_late_down



% Iterating on the size of previously calculated frM
for i = 1:size(frM,1)
    
    %Simon methods - initialization
    evolution=0;
    evolution_cdf_pre_post=0;
    %compa_after=0;
    
    % Compute all the Fr of all the phase for futur calculation
    FrM_spont_s{i}=smoothdata(frM.Fr_spont{i},'gaussian', smoothing_factor);
    FrM_early_s{i}=smoothdata(frM.Fr_early{i},'gaussian', smoothing_factor);
    FrM_middle_s{i}=smoothdata(frM.Fr_middle{i},'gaussian', smoothing_factor);
    FrM_late_s{i}=smoothdata(frM.Fr_late{i},'gaussian', smoothing_factor);
    FrM_after_s{i}=smoothdata(frM.Fr_after{i},'gaussian', smoothing_factor);
    
    % Mean of the FR spont
    base = mean(frM.Fr_spont{i}(300:599));
    % Std of Fr spont
    basestd = std(frM.Fr_spont{i}(300:599));
    
    % If std spont = 0 then assign std=0.1
    if basestd==0
        basestd=0.1;
    end
    
    % - - -
    % Z-score calculation for each phase : EARLY MID / LATE
    % - - -
    
    % intialization of the variable
    b=0;
    c=0;
    d=0;
    clearvars zdata
    mod_EandMmerge_up=[];
    mod_EandMmerge_down=[];
    
    clearvars zdata
    zdata=[];
    for j=1:size(FrM_early_s{i},2)
        zdata(j)= (FrM_early_s{i}(j) - base) / basestd;
    end
    zdata_early_full=zdata;
    
	clearvars zdata
    zdata=[];
    for j=1:size(FrM_middle_s{i},2)
        zdata(j)= (FrM_middle_s{i}(j) - base) / basestd;
    end
    zdata_mid_full=zdata;
    
    
    % Merge EARLY and MID
    EandMmerge=[FrM_early_s{i} FrM_middle_s{i}];
    
    
    for j=1:size(EandMmerge,2)
        
        % Z-score calculation
        zdata(j)= (EandMmerge(j) - base) / basestd;
        
        % If crossing the z-score threshold (either + or -) then b value increase
        if zdata(j)>threshold_zscore | zdata(j)<-threshold_zscore
            b=b+1;
        end
        
        % If crossing the z-score threshold (+) then c value increase
        if zdata(j)>threshold_zscore
            c=c+1;
            % If extraction of the modulation of the unit during the phase
            % Here is the calculation (here juste Fr phase - Fr spont)
            mod_EandMmerge_up=[(EandMmerge(j))-(base) mod_EandMmerge_up];
        end
        
        % If crossing the z-score threshold (-) then c value increase
        if zdata(j)<-threshold_zscore
            d=d+1;
            %mod_early_down=[(frY.Fr_early{i}(a)+1)/(base+1) mod_early_down];
            mod_EandMmerge_down=[(EandMmerge(j))-(base) mod_EandMmerge_down];
        end
        
    end
    
    % Percentage of zscore significant for + and - the threshold
    zdata_EandMmerge(i)=b/size(EandMmerge,2);
    % Percentage of zscore significant for + the threshold
    zdata_EandMmerge_up(i)=c/size(EandMmerge,2);
    % Percentage of zscore significant for - the threshold
    zdata_EandMmerge_down(i)=d/size(EandMmerge,2);
    % Raw zscore
    zdata_EandMmerge_full=zdata;
    
    % - - - - - - -
    % - - - Classify up/down regulation
    % - - - - - - -
    
    
    % Normal increase /// CAT = 1
    % example : z-score up > percentage threshold persistent % AND z-score down < percentage threshold intermitent
    if zdata_EandMmerge_up(i)>threshold_per_zscore & zdata_EandMmerge_down(i)<threshold_weak % ==0 %& mean(frM.Fr_early{i})>base
        zgroup_EandMmerge(i)=1;
        
        % Weak increase /// CAT = 10
    elseif zdata_EandMmerge_up(i)>threshold_weak & zdata_EandMmerge_down(i)<threshold_weak %==0 %& mean(frM.Fr_early{i})>base
        zgroup_EandMmerge(i)=10;
        
        % Not participating /// CAT = 0
    elseif zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)<threshold_weak
        zgroup_EandMmerge(i)=0;
        
        % Normal decrease /// CAT = 2
    elseif  zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)>threshold_per_zscore %& mean(frM.Fr_early{i})<base
        zgroup_EandMmerge(i)=2;
        
        % Weak decrease /// CAT = 20
    elseif zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)>threshold_weak
        zgroup_EandMmerge(i)=20;
        
        % Increase AND Decrease /// CAT = 3 or 4, depending on the order of increase/decrease
    elseif zdata_EandMmerge_up(i)>threshold_per_zscore & zdata_EandMmerge_down(i)>threshold_per_zscore
        if zdata_EandMmerge_up(i)>zdata_EandMmerge_down(i)
            zgroup_EandMmerge(i)=3;
        else
            zgroup_EandMmerge(i)=4;
        end
        
        % Not able to categorise /// CAT = 5, neurons that behave wierdly
    else
        zgroup_EandMmerge(i)=5;
    end
    
    
    % LATE
    % /!\ same methods as for the merge EARLY MID
    a=1;
    b=0;
    c=0;
    d=0;
    zdata=[];
    mod_late_up=[];
    mod_late_down=[];
    for j=1:size(FrM_late_s{i},2)
        zdata(j)= (FrM_late_s{i}(j) - base) / basestd;
        
        if zdata(j)>threshold_zscore | zdata(j)<-threshold_zscore
            b=b+1;
        end
        if zdata(j)>threshold_zscore
            c=c+1;
            mod_late_up=[(frM.Fr_late{i}(j))-(base) mod_late_up];
        end
        if zdata(j)<-threshold_zscore
            d=d+1;
            mod_late_down=[(frM.Fr_late{i}(j))-(base) mod_late_down];
        end
        
    end
    
    zdata_late(i)=b/size(FrM_late_s{i},2);
    zdata_late_up(i)=c/size(FrM_late_s{i},2);
    zdata_late_down(i)=d/size(FrM_late_s{i},2);
    zdata_late_full=zdata;
    
    % Normal increase
    if zdata_late_up(i)>threshold_per_zscore & zdata_late_down(i)<threshold_weak % ==0 %& mean(frM.Fr_early{i})>base
        zgroup_late(i)=1;
        
        % Weak increase
    elseif zdata_late_up(i)>threshold_weak & zdata_late_down(i)<threshold_weak %==0 %& mean(frM.Fr_early{i})>base
        zgroup_late(i)=10;
        
        % Not participating
    elseif zdata_late_up(i)<threshold_weak & zdata_late_down(i)<threshold_weak
        zgroup_late(i)=0;
        
        % Normal decrease
    elseif  zdata_late_up(i)<threshold_weak & zdata_late_down(i)>threshold_per_zscore %& mean(frM.Fr_early{i})<base
        zgroup_late(i)=2;
        
        % Weak decrease
    elseif zdata_late_up(i)<threshold_weak & zdata_late_down(i)>threshold_weak
        zgroup_late(i)=20;
        
        % Increase AND Decrease
    elseif zdata_late_up(i)>threshold_per_zscore & zdata_late_down(i)>threshold_per_zscore
        if zdata_late_up(i)>zdata_late_down(i)
            zgroup_late(i)=3;
        else
            zgroup_late(i)=4;
        end
        
        % Not able to categorise
    else
        zgroup_late(i)=5;
    end
    
    % Raw zscore data
    %frM.zscore_EandMmerge_full{i}=zdata_EandMmerge_full;
  	frM.zscore_early_full{i}=zdata_early_full;
    frM.zscore_mid_full{i}=zdata_mid_full;
    frM.zscore_late_full{i}=zdata_late_full;

    
    clearvars zdata_EandMmerge_full zdata_late_full
end
% 
% frM.zscore_early_full=zdata_early_full';
% frM.zscore_mid_full=zdata_mid_full';

frM.zscore_group_EandMmerge=zgroup_EandMmerge';
%frM.zscore_EandMmerge=zdata_EandMmerge';
frM.zscore_EandMmerge_up=zdata_EandMmerge_up';
frM.zscore_EandMmerge_down=zdata_EandMmerge_down';

frM.zscore_group_late=zgroup_late';
frM.zscore_late=zdata_late';
frM.zscore_late_up=zdata_late_up';
frM.zscore_late_down=zdata_late_down';

%%


clearvars zgroup_EandMmerge zdata_EandMmerge zdata_EandMmerge_up zdata_EandMmerge_down
clearvars zgroup_late zdata_late zdata_late_up zdata_late_down
% - - - - - - - -
% - - - Start calculating FrY
% - - - - - - - -

% Iterating on the size of previously calculated frM
for i = 1:size(frY,1)
    
    % Compute all the Fr of all the phase for futur calculation
    FrY_spont_s{i}=smoothdata(frY.Fr_spont{i},'gaussian', smoothing_factor);
    FrY_early_s{i}=smoothdata(frY.Fr_early{i},'gaussian', smoothing_factor);
    FrY_middle_s{i}=smoothdata(frY.Fr_middle{i},'gaussian', smoothing_factor);
    FrY_late_s{i}=smoothdata(frY.Fr_late{i},'gaussian', smoothing_factor);
    FrY_after_s{i}=smoothdata(frY.Fr_after{i},'gaussian', smoothing_factor);
    
    % Mean of the FR spont
    base = mean(frY.Fr_spont{i}(300:599));
    % Std of Fr spont
    basestd = std(frY.Fr_spont{i}(300:599));
    
    % If std spont = 0 then assign std=0.1
    if basestd==0
        basestd=0.1;
    end
    
    % - - -
    % Z-score calculation for each phase : EARLY / LATE
    % - - -
    
    % Early and Mid MERGE = EARLY
    
    clearvars zdata
    zdata=[];
    for j=1:size(FrY_early_s{i},2)
        zdata(j)= (FrY_early_s{i}(j) - base) / basestd;
    end
    zdata_early_full=zdata;
    
	clearvars zdata
    zdata=[];
    for j=1:size(FrY_middle_s{i},2)
        zdata(j)= (FrY_middle_s{i}(j) - base) / basestd;
    end
    zdata_mid_full=zdata;
    
    
    % intialization of the variable
    b=0;
    c=0;
    d=0;
    clearvars zdata
    mod_EandMmerge_up=[];
    mod_EandMmerge_down=[];
    
    % Merge EARLY and MID
    EandMmerge=[FrY_early_s{i} FrY_middle_s{i}];
    
    
    for j=1:size(EandMmerge,2)
        
        % Z-score calculation
        zdata(j)= (EandMmerge(j) - base) / basestd;
        
        % If crossing the z-score threshold (either + or -) then b value increase
        if zdata(j)>threshold_zscore | zdata(j)<-threshold_zscore
            b=b+1;
        end
        
        % If crossing the z-score threshold (+) then c value increase
        if zdata(j)>threshold_zscore
            c=c+1;
            % If extraction of the modulation of the unit during the phase
            % Here is the calculation (here juste Fr phase - Fr spont)
            mod_EandMmerge_up=[(EandMmerge(j))-(base) mod_EandMmerge_up];
        end
        
        % If crossing the z-score threshold (-) then c value increase
        if zdata(j)<-threshold_zscore
            d=d+1;
            %mod_early_down=[(frY.Fr_early{i}(a)+1)/(base+1) mod_early_down];
            mod_EandMmerge_down=[(EandMmerge(j))-(base) mod_EandMmerge_down];
        end
        
    end
    
    % Percentage of zscore significant for + and - the threshold
    zdata_EandMmerge(i)=b/size(EandMmerge,2);
    % Percentage of zscore significant for + the threshold
    zdata_EandMmerge_up(i)=c/size(EandMmerge,2);
    % Percentage of zscore significant for - the threshold
    zdata_EandMmerge_down(i)=d/size(EandMmerge,2);
    % Raw zscore
    zdata_EandMmerge_full=zdata;
    
    % - - - - - - -
    % - - - Classify up/down regulation
    % - - - - - - -
    
    
    % Normal increase /// CAT = 1
    % example : z-score up > percentage threshold persistent % AND z-score down < percentage threshold intermitent
    if zdata_EandMmerge_up(i)>threshold_per_zscore & zdata_EandMmerge_down(i)<threshold_weak % ==0 %& mean(frM.Fr_early{i})>base
        zgroup_EandMmerge(i)=1;
        
        % Weak increase /// CAT = 10
    elseif zdata_EandMmerge_up(i)>threshold_weak & zdata_EandMmerge_down(i)<threshold_weak %==0 %& mean(frM.Fr_early{i})>base
        zgroup_EandMmerge(i)=10;
        
        % Not participating /// CAT = 0
    elseif zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)<threshold_weak
        zgroup_EandMmerge(i)=0;
        
        % Normal decrease /// CAT = 2
    elseif  zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)>threshold_per_zscore %& mean(frM.Fr_early{i})<base
        zgroup_EandMmerge(i)=2;
        
        % Weak decrease /// CAT = 20
    elseif zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)>threshold_weak
        zgroup_EandMmerge(i)=20;
        
        % Increase AND Decrease /// CAT = 3 or 4, depending on the order of increase/decrease
    elseif zdata_EandMmerge_up(i)>threshold_per_zscore & zdata_EandMmerge_down(i)>threshold_per_zscore
        if zdata_EandMmerge_up(i)>zdata_EandMmerge_down(i)
            zgroup_EandMmerge(i)=3;
        else
            zgroup_EandMmerge(i)=4;
        end
        
        % Not able to categorise /// CAT = 5, neurons that behave wierdly
    else
        zgroup_EandMmerge(i)=5;
    end
    
    
    % LATE
    % /!\ same methods as for the merge EARLY MID
    a=1;
    b=0;
    c=0;
    d=0;
    zdata=[];
    mod_late_up=[];
    mod_late_down=[];
    for j=1:size(FrY_late_s{i},2)
        zdata(j)= (FrY_late_s{i}(j) - base) / basestd;
        
        if zdata(j)>threshold_zscore | zdata(j)<-threshold_zscore
            b=b+1;
        end
        if zdata(j)>threshold_zscore
            c=c+1;
            mod_late_up=[(frY.Fr_late{i}(j))-(base) mod_late_up];
        end
        if zdata(j)<-threshold_zscore
            d=d+1;
            mod_late_down=[(frY.Fr_late{i}(j))-(base) mod_late_down];
        end
        
    end
    
    zdata_late(i)=b/size(FrY_late_s{i},2);
    zdata_late_up(i)=c/size(FrY_late_s{i},2);
    zdata_late_down(i)=d/size(FrY_late_s{i},2);
    zdata_late_full=zdata;
    
    % Normal increase
    if zdata_late_up(i)>threshold_per_zscore & zdata_late_down(i)<threshold_weak % ==0 %& mean(frM.Fr_early{i})>base
        zgroup_late(i)=1;
        
        % Weak increase
    elseif zdata_late_up(i)>threshold_weak & zdata_late_down(i)<threshold_weak %==0 %& mean(frM.Fr_early{i})>base
        zgroup_late(i)=10;
        
        % Not participating
    elseif zdata_late_up(i)<threshold_weak & zdata_late_down(i)<threshold_weak
        zgroup_late(i)=0;
        
        % Normal decrease
    elseif  zdata_late_up(i)<threshold_weak & zdata_late_down(i)>threshold_per_zscore %& mean(frM.Fr_early{i})<base
        zgroup_late(i)=2;
        
        % Weak decrease
    elseif zdata_late_up(i)<threshold_weak & zdata_late_down(i)>threshold_weak
        zgroup_late(i)=20;
        
        % Increase AND Decrease
    elseif zdata_late_up(i)>threshold_per_zscore & zdata_late_down(i)>threshold_per_zscore
        if zdata_late_up(i)>zdata_late_down(i)
            zgroup_late(i)=3;
        else
            zgroup_late(i)=4;
        end
        
        % Not able to categorise
    else
        zgroup_late(i)=5;
    end
    
    % Raw zscore data
    %frY.zscore_EandMmerge_full{i}=zdata_EandMmerge_full;
  	frY.zscore_early_full{i}=zdata_early_full;
    frY.zscore_mid_full{i}=zdata_mid_full;
    frY.zscore_late_full{i}=zdata_late_full;
    
    clearvars zdata_EandMmerge_full zdata_late_full
end

frY.zscore_group_EandMmerge=zgroup_EandMmerge';
%frY.zscore_EandMmerge=zdata_EandMmerge';
frY.zscore_EandMmerge_up=zdata_EandMmerge_up';
frY.zscore_EandMmerge_down=zdata_EandMmerge_down';

frY.zscore_group_late=zgroup_late';
frY.zscore_late=zdata_late';
frY.zscore_late_up=zdata_late_up';
frY.zscore_late_down=zdata_late_down';

%% frM

% - - - - - - - -
% - - - Computation of the Depth
% - - - - - - -

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
frM.depth = depths_M';

%% frY

for i = 1:size(frY,1)
    if isnan(frY.Id_spont(i))
        load(fullfile(sort_yokes_horridge{frY.Recording(i)},'chanMap.mat'));
        clu = frY.Id_horr(i);
        idepth = get_cluster_depth(sort_yokes_horridge{frY.Recording(i)}, clu);
    else
        clu = frY.Id_spont(i);
        load(fullfile(sort_yokes_before{frY.Recording(i)},'chanMap.mat'));
        idepth = get_cluster_depth(sort_yokes_before{frY.Recording(i)}, clu);
    end
        depths_Y(i) = abs(idepth.depth' - ycoords(yokes_depth(frY.Recording(i))));

    
end
frY.depth = depths_Y';

%%
% - - - - - - - - 
% Identify Dead neurons frM
% - - - - - - - - - 

% Number of the recording studied
% recordingsM=find(frM.Recording==QQ);

% for i=1:length(recordingsM)
%     j=recordingsM(i);
%     
%     if frM.Fr_average_after(j) < 1
%         
%             figure(j)
%             plot([frM.Fr_spont{j} frM.Fr_early{j} frM.Fr_middle{j} frM.Fr_late{j} frM.Fr_after{j}])
%             hold on
%             line([600 600],[0 30],'Color','red')
%             iter=600+length(frM.Fr_early{j});
%             line([iter iter],[0 30],'Color','black')
%             iter=iter+length(frM.Fr_middle{j});
%             line([iter iter],[0 30],'Color','black')
%             iter=iter+length(frM.Fr_late{j});
%             line([iter iter],[0 30],'Color','red')
%             iter=iter+600;
%             line([iter iter],[0 30],'Color','red')
%             hold off
%         
%     end
%     
% end
%%
for j= 38%  %Line of the dying neurons
    frM.zscore_group_late(j)=3;
    frM.zscore_group_EandMmerge(j)=3;
end


%% frY 

% Number of the recording studied
% recordingsY=find(frY.Recording==33);
% 
% for i=1:length(recordingsY)
%     j=recordingsY(i);
%     
%     if frY.Fr_average_after(j) < 1
%         
%             figure(j)
%             plot([frY.Fr_spont{j} frY.Fr_early{j} frY.Fr_middle{j} frY.Fr_late{j} frY.Fr_after{j}])
%             hold on
%             line([600 600],[0 30],'Color','red')
%             iter=600+length(frY.Fr_early{j});
%             line([iter iter],[0 30],'Color','black')
%             iter=iter+length(frY.Fr_middle{j});
%             line([iter iter],[0 30],'Color','black')
%             iter=iter+length(frY.Fr_late{j});
%             line([iter iter],[0 30],'Color','red')
%             iter=iter+600;
%             line([iter iter],[0 30],'Color','red')
%             hold off
%         
%     end
%     
% end
%%
for j= 4 %Line of the dying neurons
    frY.zscore_group_late(j)=3;
    frY.zscore_group_EandMmerge(j)=3;
end




%% cleanup mathias
path_to_code = "/mnt/takeokalab/takeokalabwip2023/Mathias/";
save(path_to_code + "/code/Mathias/data/frM_" + mathias_param + ".mat", "frM", "-v7.3")
save(path_to_code + "/code/Mathias/data/frY_" + mathias_param + ".mat", "frY", "-v7.3")
clearvars; close all;
clear
end

