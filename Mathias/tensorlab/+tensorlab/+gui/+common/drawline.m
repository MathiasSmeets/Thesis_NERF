function drawline(h, pos, width)
    ax = axes(h, 'Visible', 'off', 'Units', 'pixels', 'Position', [pos width 4]);
    axis(ax, [0 1 0 1]);
    line(ax, [0 1], [0.5 0.5], 'LineWidth', 0.5, 'Color', [0.2 0.2 0.2]);
    hold(ax, 'on');
    line(ax, [0 1], [0.4 0.4], 'LineWidth', 0.5, 'Color', [1 1 1]);
    line(ax, [0 1], [0.1 0.1], 'LineWidth', 0.5, 'Color', get(h, 'BackgroundColor'));
end 
