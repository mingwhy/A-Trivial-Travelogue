function data = load_error_results(dataset, experiment_type, methods)
    data = cell(1, length(methods));
    for i = 1:length(methods)
        filename = sprintf('results/%s_%s_%s.mat', methods{i}, dataset, experiment_type);
        load(filename);
        data{i} = experiment.error_all(experiment.status_all == 1);
    end
end