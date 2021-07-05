function plot_band(lower, upper, lower_line, upper_line, plot_title, plot_xlabel, plot_ylabel, plot_legend, varargin)

% ticks interpreter
set(gca,'TickLabelInterpreter','latex')

% Plot lines and confidence band for dynamic variance decompositions

maxmin = @(x) max(min(x,1),0);

% Plot lines/bands
hold on;
if ~isempty(lower)
    the_band = plot(1:length(lower), maxmin(lower), '--k');
end
if ~isempty(upper)
    the_band = plot(1:length(upper), maxmin(upper), '--k');
end
if ~isempty(lower_line)
    the_line = plot(1:length(lower_line), maxmin(lower_line), '-k', 'LineWidth', 2);
end
if ~isempty(upper_line)
    the_line = plot(1:length(upper_line), maxmin(upper_line), '-k', 'LineWidth', 2);
end
hold off;

% Titles/labels/legend
title(plot_title, 'interpreter', 'latex');
xlabel(plot_xlabel, 'interpreter', 'latex');
ylabel(plot_ylabel, 'interpreter', 'latex');

% Set optional arguments
if ~isempty(varargin)
    for j=1:length(varargin)/2
        set(gca, varargin{2*j-1}, varargin{2*j});
    end
end

% Legend
if ~isempty(plot_legend)
    the_fontsize = get(gca, 'FontSize');
    legend([the_line the_band], plot_legend, 'FontSize', the_fontsize, 'interpreter', 'latex');
end

end