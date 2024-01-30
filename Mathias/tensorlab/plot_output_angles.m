function out_figure = plot_output_angles(u, angles, model_name)
% plot_output_angles  Plots the output angles
%
% plot_output_angles(u, angles)
%
% Author(s):    Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%               Matthieu Vendeville (matthieu.vendeville@esat.kuleuven.be)
%
% Version History:
% - 2018/05/03   MV      Initial version

tmp = cellfun(@(u) normc(u)'*normc(u), u, 'UniformOutput', false);
out_figure = figure('Name', ['Cosine of angles between rank-1 terms - ' model_name]);
% set(gcf,'WindowStyle','modal')
imagesc(angles);
xlabel('Rank-1 term');
ylabel('Rank-1 term');
axy = gca;
axy.YTick = [1:size(tmp{1},1)]; % Backwards compatibility
axy.XTick = [1:size(tmp{1},1)];
title('Cosine of angles between rank-1 terms')
colorbar;
caxis([-1 1])
