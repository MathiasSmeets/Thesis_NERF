function out_figure = plot_convergence_curves(output, TolFun, TolX, model_name)
% plot_convergence_curves  Plots solution convergence curves
%
% plot_convergence_curves(output, TolFun, TolX)
%
% Author(s):    Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%               Matthieu Vendeville (matthieu.vendeville@esat.kuleuven.be)
%
% Version History:
% - 2018/05/03   MV      Initial version


out_figure = figure('Name',['Convergence algorithm - ' model_name]);
% set(gcf,'WindowStyle','modal');
if ~isfield(output.Refinement, 'relerr')
    semilogy(output.Algorithm.fval); hold on
    %Relative Step Sizes
    semilogy(output.Algorithm.relfval);
    semilogy(output.Algorithm.relstep);
    stoprgb = output.Algorithm.info;
else
    semilogy(output.Refinement.fval); hold on
    %Relative Step Sizes
    semilogy(output.Refinement.relfval);
    semilogy(output.Refinement.relstep);
    stoprgb = output.Refinement.info;
end
colorstmp = get(gca,'colororder');
line(get(gca,'XLim'),[TolFun TolFun],'LineStyle','--','Color',colorstmp(2,:))
line(get(gca,'XLim'),[TolX TolX],'LineStyle','--','Color',colorstmp(3,:))
legend({'fval','relfval','relstep'})
xlabel('Iteration')
xs = get(gca,'XLim');
text(xs(1),TolFun,'TolFun','Color',colorstmp(2,:),'VerticalAlignment', 'bottom')
text(xs(1),TolX,'TolX','Color',colorstmp(3,:),'VerticalAlignment', 'bottom')
switch stoprgb
    case 1
        title('Objective function tolerance reached')
    case 2
        title('Step size tolerance reached')
    case 3
        title('Maximum number of iterations reached')
end
       