function save_fig(folder,filename)

    saveas(gcf, fullfile(folder, strcat(filename, '.fig')));
    saveas(gcf, fullfile(folder, strcat(filename, '.png')));
    saveas(gcf, fullfile(folder, strcat(filename, '.eps')), 'epsc');

end