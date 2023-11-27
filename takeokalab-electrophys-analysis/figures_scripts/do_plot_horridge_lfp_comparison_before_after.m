clear all
%% Parameters and data
fignum = 3; % 1: example average PSD; 2: summed by freq band (barplot); 3: summed by fre band (scatter plot)

normband = 1;
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
bands = [1:6];
allch = 1;

%Master mouse #1 --> good master (1)
%Master mouse #2 --> bad master (1)
%Master mouse #3 --> good master (2)
%Master mouse #4 --> bad master (2)
%Master mouse #5 --> bad master (3)
%Master mouse #6 --> noisy recording, excluded
%Master mouse #7 --> good master (3)
%Master mouse #8 --> noisy recording, excluded

master{1} = 'D:\Experiments\Data\SC\20200907\psdmean_before_after_allch.mat';
master{2} = 'D:\Experiments\Data\SC\20200921\Master_3\psdmean_before_after_allch.mat';
master{3} = 'D:\Experiments\Data\SC\20201002\Master 7\psdmean_before_after_allch.mat';
master{4} = 'D:\Experiments\Data\SC\20201012\Master9\psdmean_before_after_allch.mat';
% master{5} = 'D:\Experiments\Data\SC\20201014\Master10\psdmean_before_after_allch.mat';

%Master2 recording did not have the recovery period, if re-insterted
%decomment the if cycle in the bad_master section below
% master_bad{1} = 'D:\Experiments\Data\SC\20200916\Master 2\psdmean_before_after_allch.mat';
% master_bad{1} = 'D:\Experiments\Data\SC\20200921\Master_4\psdmean_before_after_recovery_allch.mat';
% master_bad{2} = 'D:\Experiments\Data\SC\20200923\Master_5\psdmean_before_after_recovery_allch.mat';

master_bad = {};

yoked{1} = 'D:\Experiments\Data\SC\20200914\psdmean_before_after_allch.mat';
yoked{2} = 'D:\Experiments\Data\SC\20200916\Yoked 3\psdmean_before_after_allch.mat';
yoked{3} = 'D:\Experiments\Data\SC\20200923\Yoked_4\psdmean_before_after_allch.mat';
yoked{4} = 'D:\Experiments\Data\SC\20200928\Yoked 5\psdmean_before_after_allch.mat';
yoked{5} = 'D:\Experiments\Data\SC\20201005\psdmean_before_after_allch.mat';

%% Organize data of master mice
disp('loading a processing master')

sizemean = [];
psdmean_all = {};
fall = {};
for i = 1:length(master)
    load(master{i})
    clear psdmean_temp
    if allch
        for j = 1:length(psdmean)
            psdmean_temp{1}(j,:) = psdmean{j}(:,1)';
            psdmean_temp{2}(j,:) = psdmean{j}(:,2)';
        end
        psdmean_all{i} = [mean(psdmean_temp{1}); mean(psdmean_temp{2})];
        fall{i} = f';
        sizemean = [sizemean size(psdmean_all{i},2)];
    else
        psdmean_all{i} = psdmean{6}';
        fall{i} = f';
        sizemean = [sizemean size(psdmean_all{i},2)];
    end
end
clear psdmean
[maxL, indmax] = max(sizemean);
for i = 1:length(psdmean_all)
    if sizemean(i)<maxL
        for j = 1:2
            psdmean{i}(j,:) = interp1(fall{i}, psdmean_all{i}(j,:), fall{indmax});
        end
    else
        for j = 1:2
            psdmean{i}(j,:) = psdmean_all{i}(j,:);
        end
    end
end
f = fall{indmax};

for i = 1:length(master)
    switch fignum
        case 1
            P_master_B(i,:) = psdmean{i}(1,:)./max(psdmean{i}(1,:));
            P_master_A(i,:) = psdmean{i}(2,:)./max(psdmean{i}(2,:));
        case {2,3}
            P_master_B(i,:) = psdmean{i}(1,:);
            P_master_A(i,:) = psdmean{i}(2,:);
    end
    f_master{i} = f;
end

%% Organize data of bad masters
if ~isempty(master_bad)
    disp('loading and processing bad master')
    sizemean = [];
    clear psdmean_all fall f 
    for i = 1:length(master_bad)
        load(master_bad{i})
       clear psdmean_temp
        if allch
            for j = 1:length(psdmean)
                psdmean_temp{1}(j,:) = psdmean{j}(:,1)';
                psdmean_temp{2}(j,:) = psdmean{j}(:,2)';
            end
            psdmean_all{i} = [mean(psdmean_temp{1}); mean(psdmean_temp{2})];
            fall{i} = f';
            sizemean = [sizemean size(psdmean_all{i},2)];
        else
            psdmean_all{i} = psdmean{6}';
            fall{i} = f';
            sizemean = [sizemean size(psdmean_all{i},2)];
        end
    end
    clear psdmean
    [maxL, indmax] = max(sizemean);
    for i = 1:length(psdmean_all)
        if sizemean(i)<maxL
            for j = 1:2
                psdmean{i}(j,:) = interp1(fall{i}, psdmean_all{i}(j,:), fall{indmax});
            end
            if i ~= 1 %only if Master2 recording included
                psdmean{i}(3,:) = interp1(fall{i}, psdmean_all{i}(3,:), fall{indmax});
            end
        else
            for j = 1:2
                psdmean{i}(j,:) = psdmean_all{i}(j,:);
            end
            if i ~= 1 %only if Master2 recording included
                psdmean{i}(3,:) = psdmean_all{i}(3,:);
            end
        end
    end
    f = fall{indmax};
    for i = 1:length(master_bad)
        switch fignum
            case 1
                P_master_bad_B(i,:) = psdmean{i}(1,:)./max(psdmen{i}(1,:));
                P_master_bad_A(i,:) = psdmean{i}(2,:)./max(psdmean{i}(2,:));
                if i ~= 1
                    P_master_bad_R(i-1,:) = psdmean{i}(3,:)./max(psdmean{i}(3,:));
                end
                P_master_bad_R(i,:) = psdmean{i}(3,:)./max(psdmean{i}(3,:));
            case {2,3}
                P_master_bad_B(i,:) = psdmean{i}(1,:);
                P_master_bad_A(i,:) = psdmean{i}(2,:);
                if i ~= 1 %only if Master2 recording included
                    P_master_bad_R(i-1,:) = psdmean{i}(3,:);
                end
    %                 P_master_bad_R(i,:) = psdmean{i}(3,:);
        end
        f_master_bad{i} = f;
    end
end

%% Organize data of yoked mice
disp('loading and processing yoked')
clear psdmean_all fall f
sizemean = [];
psdmean_all = {};
fall = {};
for i = 1:length(yoked)
    load(yoked{i})
    clear psdmean_temp
    if allch
        for j = 1:length(psdmean)
            psdmean_temp{1}(j,:) = psdmean{j}(:,1)';
            psdmean_temp{2}(j,:) = psdmean{j}(:,2)';
        end
        psdmean_all{i} = [nanmean(psdmean_temp{1}); nanmean(psdmean_temp{2})];
        fall{i} = f';
        sizemean = [sizemean size(psdmean_all{i},2)];
    else
        psdmean_all{i} = psdmean{6}';
        fall{i} = f';
        sizemean = [sizemean size(psdmean_all{i},2)];
    end
end
clear psdmean
[maxL, indmax] = max(sizemean);
for i = 1:length(psdmean_all)
    if sizemean(i)<maxL
        for j = 1:2
            psdmean{i}(j,:) = interp1(fall{i}, psdmean_all{i}(j,:), fall{indmax});
        end
    else
        for j = 1:2
            psdmean{i}(j,:) = psdmean_all{i}(j,:);
        end
    end
end
f = fall{indmax};
for i = 1:length(psdmean)
    switch fignum
        case 1
            P_yoked_B(i,:) = psdmean{i}(1,:)./max(psdmean{i}(1,:));
            P_yoked_A(i,:) = psdmean{i}(2,:)./max(psdmean{i}(2,:));
        case {2,3}
            P_yoked_B(i,:) = psdmean{i}(1,:);
            P_yoked_A(i,:) = psdmean{i}(2,:);
    end
    f_yoked{i} = f;
end


% P_Master_before = mean(P_master_B);
% S_Master_before = std(P_master_B);
% P_Master_after = mean(P_master_A);
% S_Master_after = std(P_master_A);
% 
% P_Master_bad_before = mean(P_master_bad_B);
% S_Master_bad_before = std(P_master_bad_B);
% P_Master_bad_after = mean(P_master_bad_A);
% S_Master_bad_after = std(P_master_bad_A);

P_Master_before = mean(P_master_B);
S_Master_before = std(P_master_B);
P_Master_after = mean(P_master_A);
S_Master_after = std(P_master_A);

if ~isempty(master_bad)
    P_Master_bad_before = mean(P_master_bad_B);
    S_Master_bad_before = std(P_master_bad_B);
    P_Master_bad_after = mean(P_master_bad_A);
    S_Master_bad_after = std(P_master_bad_A);
    P_Master_bad_recovery = mean(P_master_bad_R);
    S_Master_bad_recovery = std(P_master_bad_R);
end

P_Yoked_before = mean(P_yoked_B);
S_Yoked_before = std(P_yoked_B);
P_Yoked_after = mean(P_yoked_A);
S_Yoked_after = std(P_yoked_A);




switch fignum
    %% PLOT TYPE 1
    case 1
        figure
        subplot(1,3,1)
        hold on
            %normalize for representation
        plot(f_master{1},10*log10(P_Master_before), 'color', colors(1,:));
        plot(f_master{1},10*log10(P_Master_after), 'color', colors(2,:));
        xlabel('f [Hz]')
        ylabel('PSD [db/Hz]')
        xlim([0 100])
        legend({'before','after'})
        title('Master')


        subplot(1,3,2)
        hold on
        plot(f_master_bad{1},10*log10(P_Master_bad_before), 'color', colors(1,:));
        plot(f_master_bad{1},10*log10(P_Master_bad_after), 'color', colors(2,:));
        xlabel('f [Hz]')
        ylabel('PSD [db/Hz]')
        xlim([0 100])
        legend({'before','after'})
        title('Bad Master')

        subplot(1,3,3)
        hold on

        x = f_yoked{1}';
        y = 10*log10(P_Yoked_before);
        curve1 = y + S_Yoked_before;
        curve2 = y - S_Yoked_before;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        h = fill(x2, inBetween, colors(1,:));
        set(h,'facealpha',.2)
        plot(x, y, 'color',colors(1,:), 'LineWidth', 2); 

        x = f_yoked{1}';
        y = 10*log10(P_Yoked_after);
        curve1 = y + S_Yoked_after;
        curve2 = y - S_Yoked_after;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        h = fill(x2, inBetween, colors(2,:));
        set(h,'facealpha',.2)
        plot(x, y, 'color',colors(2,:), 'LineWidth', 2); 


        xlabel('f [Hz]')
        ylabel('PSD [db/Hz]')
        xlim([0 100])
        legend({'before','after'})
        title('Yoked')
        
        
    case {2, 3}
        labs = {'0.5-4 Hz','4-8 Hz','8-13 Hz','13-30 Hz','30-60 Hz', '60-100 Hz'};
        cm = hot;

      %% PLOT TYPE 2
        %%%%%%%%% MASTER %%%%%%%%%%%
        
        B_master_B = lfp_bandsum_comparison(P_master_B,f_master{1},'plt',0);
        B_master_A = lfp_bandsum_comparison(P_master_A,f_master{1},'plt',0);
        if normband
            for i = 1:size(B_master_A,2)
                B_master_A(:,i) = B_master_A(:,i)./B_master_B(:,i);
                B_master_B(:,i) = B_master_B(:,i)./B_master_B(:,i);
            end
        end
        BM(:,1) = mean(B_master_B,2);
        BMstd(:,1) = std(B_master_B')/sqrt(size(B_master_B,2));
        BM(:,2) = mean(B_master_A,2);
        BMstd(:,2) = std(B_master_A'/sqrt(size(B_master_A,2)));  
        
        BM = BM(bands,:);
        BMstd = BMstd(bands,:);
        
        figure, hold on
        if fignum == 2
            bar(1:length(bands),BM);
            set(gca, 'XTickLabel',labs(bands))
            ngroups = size(BM, 1);
            nbars = size(BM, 2);
            % Calculating the width for each bar group
            groupwidth = min(0.8, nbars/(nbars + 1.5));
            for i = 1:nbars
                x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                errorbar(x, BM(:,i), BMstd(:,i), '.');
            end
             hold off
            legend('before','after')
        elseif fignum == 3
            bandcolors = cm(1:floor(size(cm,1)/length(bands)):end,:);
            for i = 1:length(bands)
                errorbar(1:2, BM(i,:), BMstd(i,:), 'o-', 'color',bandcolors(i,:))
            end
            hold off
            set(gca, 'XTick', [1 2], 'XtickLabel', {'before','after'})
            legend(labs)
        end
       
        ylabel('$\sum$ (PSD) [au]','Interpreter','Latex')
        title('Master')
        
        if ~isempty(master_bad)
            %%%%%%%% MASTER BAD %%%%%%%%%%%
            B_master_bad_B = lfp_bandsum_comparison(P_master_bad_B,f_master_bad{1},'plt',0);
            B_master_bad_A = lfp_bandsum_comparison(P_master_bad_A,f_master_bad{1},'plt',0);
            B_master_bad_R = lfp_bandsum_comparison(P_master_bad_R,f_master_bad{1},'plt',0);
            if normband
                for i = 1:size(B_master_bad_B,2)
                    B_master_bad_A(:,i) = B_master_bad_A(:,i)./B_master_bad_B(:,i);
                    if i ~= 1 %only if bad_master n.1 is included
                        B_master_bad_R(:,i-1) = B_master_bad_R(:,i-1)./B_master_bad_B(:,i);
                    end
                    B_master_bad_B(:,i) = B_master_bad_B(:,i)./B_master_bad_B(:,i);
                end

            end
            BMB(:,1) = mean(B_master_bad_B,2);
            BMBstd(:,1) = std(B_master_bad_B')/sqrt(size(B_master_bad_B,2));
            BMB(:,2) = mean(B_master_bad_A,2);
            BMBstd(:,2) = std(B_master_bad_A')/sqrt(size(B_master_bad_A,2));  
            BMB(:,3) = mean(B_master_bad_R,2);
            BMBstd(:,3) = std(B_master_bad_R')/sqrt(size(B_master_bad_R,2));

            BMB = BMB(bands,:);
            BMBstd = BMBstd(bands,:);
        
            figure, hold on
            bar(1:length(bands),BMB);
            set(gca, 'XTickLabel',labs(bands))

            ngroups = size(BMB, 1);
            nbars = size(BMB, 2);
            % Calculating the width for each bar group
            groupwidth = min(0.8, nbars/(nbars + 1.5));
            for i = 1:nbars
                x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                errorbar(x, BMB(:,i), BMBstd(:,i), '.');
            end
            hold off
            legend('before','after')
            ylabel('$\sum$ (PSD) [au]','Interpreter','Latex')
            title('Bad master')

        
        end
        
        % YOKED
        
        clear PP
        B_yoked_B = lfp_bandsum_comparison(P_yoked_B,f_yoked{1},'label_names',{'before','after'},'plt',0);
        B_yoked_A = lfp_bandsum_comparison(P_yoked_A,f_yoked{1},'label_names',{'before','after'},'plt',0);
        if normband
            for i = 1:size(B_yoked_A,2)
                B_yoked_A(:,i) = B_yoked_A(:,i)./B_yoked_B(:,i);
                B_yoked_B(:,i) = B_yoked_B(:,i)./B_yoked_B(:,i);
            end
        end
        BY(:,1) = mean(B_yoked_B,2);
        BYstd(:,1) = std(B_yoked_B')/sqrt(size(B_yoked_B,2));
        BY(:,2) = mean(B_yoked_A,2);
        BYstd(:,2) = std(B_yoked_A')/sqrt(size(B_yoked_A,2));
        
        BY = BY(bands,:);
        BYstd = BYstd(bands,:);
        
        figure, hold on
        if fignum == 2
            bar(1:length(bands),BY);
            set(gca, 'XTickLabel',labs(bands))
            ngroups = size(BY, 1);
            nbars = size(BY, 2);
            % Calculating the width for each bar group
            groupwidth = min(0.8, nbars/(nbars + 1.5));
            for i = 1:nbars
                x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                errorbar(x, BY(:,i), BYstd(:,i), '.');
            end
            hold off
            legend('before','after')
        elseif fignum == 3
            bandcolors = cm(1:floor(size(cm,1)/length(bands)):end,:);
            for i = 1:length(bands)
                errorbar(1:2, BY(i,:), BYstd(i,:), 'o-', 'color',bandcolors(i,:))
            end
            hold off
            set(gca, 'XTick', [1 2], 'XtickLabel', {'before','after'})
            legend(labs)
        end
        ylabel('$\sum$ (PSD) [au]','Interpreter','Latex')
        title('Yoked')
    end
    

        