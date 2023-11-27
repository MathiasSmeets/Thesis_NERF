function Neu_summary_plot(Neurons, selNeu, stim)

% Plot PSTH and info on latencies and SNR of the selected neurons
%
% Mattia D'Andola
% Jan 2020

if nargin < 2
    selNeu = 1:length(Neurons.NeuId);
end

if nargin < 3
    stim = [100 500]; %period of stimulation
end

for i = selNeu
    
    
    neuDB = Neurons(i, :);
    fig = figure;
    plot(neuDB.Psth, 'b', 'LineWidth', 2);
    hold on
    plot(neuDB.Psth + neuDB.StdDev, '-.c')
    plot(neuDB.Psth - neuDB.StdDev, '-.c')
    line([neuDB.Window(1), neuDB.Window(1)], ylim, 'Color', 'g', 'LineStyle', '--')
    line([neuDB.Window(2), neuDB.Window(2)], ylim, 'Color', 'g', 'LineStyle', '--')
    line(xlim, [neuDB.Background,...
        neuDB.Background] / 50 * 1000, 'Color', 'black')
    a = ylim;
    line(stim, [a(1)-1 a(1)-1], 'color', 'blue', 'LineWidth', 3)
    xlabel('t (ms)')
    ylabel('Rate (Hz)')
    title(sprintf('Latency: %.2f ms\nSNR: %.2f dB', neuDB.Window(1),...
        10 * log10(max(neuDB.Psth) / median(abs(neuDB.Psth)))));
    
end