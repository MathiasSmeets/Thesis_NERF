function h = plot_angles(U)
% plot_output_angles  Plots the output angles
%
% plot_output_angles(u, angles)
%
% Author(s):    Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/23   NV      Initial version

nrm = cellfun(@(u) sqrt(sum(abs(u).^2,1)), U, 'UniformOutput', false);
U   = cellfun(@(u,n) bsxfun(@rdivide, u, n), U, nrm, 'UniformOutput', false);
W   = cellfun(@(u) real(u'*u), U, 'UniformOutput', false);
congruence = prod(cat(3,W{:}),3);
angles     = abs(acos(congruence)) / pi * 180;
imagesc(angles);
xlabel('Rank-1 term');
ylabel('Rank-1 term');
R = size(U{1},2);
set(gca, 'YTick', 1:R);
set(gca, 'XTick', 1:R);
title('Angles (degrees) between rank-1 terms');
colormap(jet);
colorbar;
caxis([0 180])
