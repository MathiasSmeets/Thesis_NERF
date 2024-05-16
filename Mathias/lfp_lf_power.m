clear;clc;close all;
%load and store neuropixel lfp data
load("/mnt/takeokalab/takeokalabwip2023/Mathias/switch_data/LF_signals/LF_1m.mat");


seg = [];
%%downsample frm 2500Hz to 500HZ
reduction_factor = 5;
for i = 1:fix(size(data,2)/reduction_factor)
    for j = 1:385
        seg(j,i) = mean(data(j,1+i*reduction_factor-reduction_factor:i*reduction_factor));
    end
end


%smooth data
Fs = 500;
seg_smoothed = smoothdata(seg')';


seg_o = seg;
seg = seg_smoothed;

               
T = 1/Fs;             % Sampling period       
L = 8000;             % Length of signal in ms
t = (0:L-1)*T;        % Time vector
max_freq = 50;         % max frequency to keep for plotting, in Hz
L_ = L*Fs/1000;        %number of value to take for L size timeBin at Fs Hz
%initialize power spectrum matrix

 % - - - - - Initialization of the frequency bands of interest

PSD_over_time_Delta = [];
    PSD_over_time_Theta = [];
    PSD_over_time_Alpha = [];
    PSD_over_time_Beta = [];
    PSD_over_time_Gamma = [];
    PSD_over_time_more = [];
for j = 1:size(seg,1)


    for i = 1:size(seg,2)/(L_)
    
        
        X=seg(j,1+((L_)*i)-(L_):((L_)*i));
        
        %%plot(1000*t(1:4000),X(1:4000))
        %%title("EEG trace")
        %%xlabel("t (milliseconds)")
        %%ylabel("X(t)")
        P1=[];
        P2=[];
        Y=[];
        f=[];
        Y = fft(X);
        P2 = abs(Y/(L_));
        P1 = P2(1:(L_)/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        
        f = Fs*(0:((L_)/2))/(L_);
        %plot(f,P1)
        
        %title("Single-Sided Amplitude Spectrum of X(t)")
        %xlabel("f (Hz)")
        %ylabel("|P1(f)|")
        
        Y = fft(X);
        P2 = abs(Y/(L_));
        P1 = P2(1:(L_));
        P1(2:end-1) = 2*P1(2:end-1);
        
        
        
        %plot(f(1,1:1500),P1) 
        %title("Single-Sided Amplitude Spectrum of S(t)")
        %xlabel("f (Hz)")
        %ylabel("|P1(f)|")
        
    PSD_over_time_Delta(j,i) = mean(P1(1,1:4*4));
    PSD_over_time_Theta(j,i) = mean(P1(1,4*4:8*4));
    PSD_over_time_Alpha(j,i) = mean(P1(1,8*4:12*4));
    PSD_over_time_Beta(j,i) = mean(P1(1,12*4:30*4));
    PSD_over_time_Gamma(j,i) = mean(P1(1,30*4:50*4));
    PSD_over_time_more(j,i) = mean(P1(1,50*4:300*4));
    end
end
figure
h = heatmap(PSD_over_time_Delta(:,:),'Colormap',spring, 'GridVisible','off','ColorLimits',[0 45]);

figure
h = heatmap(PSD_over_time_Theta(:,:),'Colormap',spring, 'GridVisible','off','ColorLimits',[0 30]);
figure
h = heatmap(PSD_over_time_Alpha(:,:),'Colormap',spring, 'GridVisible','off','ColorLimits',[0 30]);
figure
h = heatmap(PSD_over_time_Beta(:,:),'Colormap',spring, 'GridVisible','off','ColorLimits',[0 30]);
figure
h = heatmap(PSD_over_time_Gamma(:,:),'Colormap',spring, 'GridVisible','off','ColorLimits',[0 30]);
figure
h = heatmap(PSD_over_time_more(:,:),'Colormap',spring, 'GridVisible','off','ColorLimits',[0 30]);