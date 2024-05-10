function Spikes_raster_plot(data,Neurons, neuSel, interval)    

% Plot spike raster plot for the units contained in the Neurons DB (or a
% subset) in a specified interval. 
% Use interval = [0 0] to get an aligned raster plot for the whole
% recording

% if nargin<3 || isempty(neuSel)
%     neuSel = 1:length(Neurons.NeuId);
% end
% if nargin<4
%     interval = [1 20];
% end
chgap = -0.25;
Fs = 40000;
figure
hold on

% if prova == [0 0]
%     for i = neuSel
%         neuDB = Neurons(i, :);
%         
%         idx = find(neuDB.SpikeTime{:}
% else
time = data.continuous.time; %time in [s]
stimbin = zeros(1,length(time));
for i = 1:length(data.events.onsets)
    idxbin = find(time>=data.events(1).onsets(i) & time<=data.events(1).offsets(i));
    stimbin(idxbin) = ones(1,length(idxbin));
end
    
idtime = find(time>=interval(1) & time<=interval(2));
plot(time(idtime), stimbin(idtime), 'b', 'linewidth', 2)
ynames{1} = 'Optostim';

cont = 1;
for i = neuSel
    neuDB = Neurons(i, :);

    idx = find(neuDB.SpikeTime{:}>interval(1) & neuDB.SpikeTime{:}<interval(2));
    numspikes = length(idx);
    xx = ones(3*numspikes,1)*nan;
    yy = ones(3*numspikes,1)*nan; 

    xx(1:3:3*numspikes)=neuDB.SpikeTime{1}(idx);
    xx(2:3:3*numspikes)=neuDB.SpikeTime{1}(idx);

    yy(1:3:3*numspikes) = (i)*chgap;
    yy(2:3:3*numspikes) = yy(1:3:3*numspikes)+0.2;
    plot(xx, yy, 'k', 'linewidth', 1);
    cont = cont+1;
    ynames{cont} = sprintf('Neu %d', neuDB.NeuId);
end
set(gca, 'ytick', [(chgap*(length(neuSel))):-chgap:chgap,0.5], 'yticklabel', fliplr(ynames)) 
%end