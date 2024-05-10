load("plots_Jeremy/psds/masters/20220525.mat")
%%
bound = 2;
low_bound = -bound;
high_bound = bound;
 %%
close all
fig = figure(1);
%fig.WindowState = 'maximized';
sgtitle("Power in corresponding frequency bands")
subplot(5,1,1)
if cond < 2
    psd_freq_no_training = psd_freq_1(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
    time_to_analyse = -timeframe(3)*60+timeframe(6)*60;
else
    time_to_analyse = timeframe(end)*60;
    psd_freq_no_training = psd_freq_1;
end
imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');
% xticks(1:endrest/5/binSize:endrest/binSize)
% xticklabels(1:endrest/60/5:endrest/60)
caxis([low_bound high_bound])
a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
title("0.5-4Hz")
xlabel("time [s]")
ylabel("Channel depth [µm]")

subplot(5,1,2)
if cond < 2
    psd_freq_no_training = psd_freq_2(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
else
    psd_freq_no_training = psd_freq_2;
end
imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');
caxis([low_bound high_bound])
a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
title("4-8Hz")
xlabel("time [s]")
ylabel("Channel depth [µm]")

subplot(5,1,3)
if cond < 2
    psd_freq_no_training = psd_freq_3(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
else
    psd_freq_no_training = psd_freq_3;
end
imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');
caxis([low_bound high_bound])
a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
title("8-12Hz")
xlabel("time [s]")
ylabel("Channel depth [µm]")

subplot(5,1,4)
if cond < 2
    psd_freq_no_training = psd_freq_4(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
else
    psd_freq_no_training = psd_freq_4;
end
imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');
caxis([low_bound high_bound])
a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
title("12-35Hz")
xlabel("time [s]")
ylabel("Channel depth [µm]")

subplot(5,1,5)
if cond < 2
    psd_freq_no_training = psd_freq_5(ceil(timeframe(3)*60/25)+1:floor(timeframe(6)*60/25),:);
else
    psd_freq_no_training = psd_freq_5;
end
imagesc(0:time_to_analyse/size(psd_freq_no_training,1):time_to_analyse-1,1:size(psd_freq_no_training,2)*10,flip(psd_freq_no_training,2)');


caxis([low_bound high_bound])
a = colorbar; ylabel(a,"Power Z-score",'Rotation',270); a.Label.Position(1) = 3;
title("35-200")
xlabel("time [s]")
ylabel("Channel depth [µm]")

 %% Save image
 
if cond == 2
    saveas(fig,"plots_Jeremy/psds/masters/"+lfpFilename(43:50)+".fig")
    saveas(fig,"plots_Jeremy/psds/masters/"+lfpFilename(43:50)+".emf")
    delete(fig)
    save("plots_Jeremy/psds/masters/"+lfpFilename(43:50)+".mat")
elseif cond == 3
    saveas(fig,"plots_Jeremy/psds/yokes/"+lfpFilename(43:50)+".fig")
    saveas(fig,"plots_Jeremy/psds/yokes/"+lfpFilename(43:50)+".emf")
    delete(fig)
    save("plots_Jeremy/psds/yokes/"+lfpFilename(43:50)+".mat")
end