function [datafil, t, MV, dataf] = moving_average(data, time, WindowSize, StepSize, chs, fs)

%@LFP/MOVING_AVERAGE compute moving average on lfp objects, with
%'WindowSize' and 'StepSize' parameters
%
%Mattia D'Andola, Jun 2013

if nargin<5
    fs = 30000;
end
dt = 1/fs;
SampleNum = floor((size(data,2) - WindowSize)/StepSize) + 1;
t = (0:SampleNum-1) * dt * StepSize + mean(time(1:WindowSize));

for ch = chs   
    dataf = zeros(WindowSize, SampleNum);
    for k = 1:WindowSize
        dataf(k,:) = data(ch,k:StepSize:end-WindowSize+k);
    end
    datafil(ch,:) = mean(dataf);
end
% l(i).value = tempvalue;
% l(i).time = t;
MV = std(dataf).^2;
