function save_fig(folder,filename)

    % Save figure in various formats

%     saveas(gcf, fullfile(folder, strcat(filename, '.fig')));
%     saveas(gcf, fullfile(folder, strcat(filename, '.png')));
    saveas(gcf, fullfile(folder, strcat(filename, '.eps')), 'epsc');

end