function build_error_boxplots_together(datasets, methods, dataset_labels, ymaxs, dpi)
    
    experiment_types = {'sim_all', 'sim_intra' };    
    sublabels = {'a)', 'b)','c)','d)','e)','f)'};
    
    k = 0;
    for i = 1:length(experiment_types)
        for j = 1:length(datasets)
            data = load_error_results(datasets{j}, experiment_types{i}, methods);
            
            k = k+1; subplot(length(experiment_types),length(datasets),k);
            box_plotting(data, methods, ymaxs(k));
            ylabel('normalized error')
            title([sublabels{k} '        ' dataset_labels{j} '          '])
        end
    end
    
    set(gcf,'PaperUnits', 'points')
    set(gcf,'PaperPosition', [0 0 600 400])
    print('-dtiff', dpi, 'images/error_boxplots_all.tiff');
    close
end

