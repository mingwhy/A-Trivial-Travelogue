%% MAIN FILE:
%
% Use this to run all the tests in the paper.
% Different tests can be run selectively.
% For convenience all results are stored as .mat files.
%
% Author: Daniel Machado, 2013

addpath('tests')
addpath('utils')
addpath('wrappers')
addpath('plotting')
addpath('/Users/ming/Documents/bioinfo_software/GEM_metabolic.models/RAVEN-main')


%% PROBLEM SETUP

% select dataset and organism (see load_dataset.m for details)
%DATASET = 'ishii'; ORGANISM = 'ecoli';
%DATASET = 'holm'; ORGANISM = 'ecoli';
DATASET = 'rintala'; ORGANISM = 'yeast';

% select experiment type
experiment_type = 'sim_all'; % simulate all fluxes 
%experiment_type = 'sim_intra'; % simulate intracellular fluxes 
%experiment_type = 'sim_secr'; % simulate secretion and growth rates

%COBRA_SOLVER = 'gurobi5';
COBRA_SOLVER = 'gurobi';

%initialize cobra
changeCobraSolver(COBRA_SOLVER, 'all');

%initialize cmpi (requires MADE package; see README)
%cmpi.init()

%load model
model = load_model(ORGANISM);
model = creategrRulesField(model); %https://opencobra.github.io/cobratoolbox/stable/modules/base/utilities/index.html
model.grRules{10}

%model521 = readCbModel('models/yeast_5.21_MCISB.xml')
%model521 = creategrRulesField(model521); %https://opencobra.github.io/cobratoolbox/stable/modules/base/utilities/index.html
%model521.grRules{10}

% add `rev` back or not, see https://github.com/SysBioChalmers/RAVEN/issues/184
%Correct rev vector: true if LB < 0 & UB > 0, or it is an exchange reaction:
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
        model.rev(i) = true;
    end
end

%load data set
dataset = load_dataset(DATASET);

options.experiment_type = experiment_type;
options.reestimate_data = true; % fit experimental data to model
%options.reestimate_data = false; 

%% BENCHMARKING

%methods = {'pFBA', 'iMAT', 'E-Flux','Lee-12', 'GIMME','RELATCH', 'MADE', 'GX-FBA'};
%methods = {'iMAT', 'Lee-12'}; % Lee-12, slower
methods = {'pFBA', 'iMAT', 'E-Flux','Lee-12'};
methods = {'iMAT'}; %GIMME, MADE, RELATCH,GX-FBA not work


for i = 1:length(methods)
    %function experiment = benchmark_method(method, model, dataset, options)
	benchmark_method(methods{i}, model, dataset, options);
end


%% SENSITIVITY ANALYSIS

conf = cell(1,6);
conf{1} = {'GIMME', 'OBJ_FRAC', 0, 1, 'lin'};
conf{2} = {'GIMME', 'GIMME_LOWER_QUANTILE', 0, 1, 'lin'}; 
conf{3} = {'iMAT', 'IMAT_LOWER_QUANTILE', 0, 0.75, 'lin'}; 
conf{4} = {'iMAT', 'IMAT_UPPER_QUANTILE', 0.25, 1, 'lin'}; 
conf{5} = {'iMAT', 'IMAT_EPS', -2, 2, 'log'}; 
conf{6} = {'MADE', 'OBJ_FRAC', 0, 1, 'lin'}; 
 
points = 100; 
 
for i = 1:6
    [method, parameter, min_val, max_val, scale] = deal(conf{i}{:}); 
    sensitivity_analysis(model, dataset, method, parameter, min_val, max_val, points, scale, options);
end


%% ROBUSTNESS ANALYSIS

%methods = {'GIMME', 'iMAT', 'E-Flux', 'Lee-12', 'RELATCH', 'MADE', 'GX-FBA'};
%methods = {'pFBA', 'iMAT', 'E-Flux','Lee-12'};
methods = {'iMAT'};
steps = 2; points = 2;

for i = 1:length(methods)
     robustness_analysis(model, dataset, methods{i}, dataset.conditions{2}, dataset.conditions{1}, steps, points, options);
end

%% PLOTTING

%build_figures()
 DPI = '-r300';
    methods = {'pFBA', 'GIMME', 'iMAT', 'MADE', 'E-Flux', 'Lee-12', 'RELATCH', 'GX-FBA'};
    datasets = {'ishii', 'holm', 'rintala'};

    dataset_labels = {'{\it E. coli} (Ishii 2007)', '{\it E. coli} (Holm 2010)', 'Yeast (Rintala 2009)'};
    ymaxs = [3 3 3 3 3 3];
    build_error_boxplots_together(datasets, methods, dataset_labels, ymaxs, DPI)
 
    datasets2 = {'ishii','ishii-protein','rintala-red','rintala-protein'};
    dataset_labels2 = {'Ishii (transcript)', 'Ishii (protein)', 'Rintala (transcript)', 'Rintala (protein)'};
    ymaxs = [4 4 4 4 4 4 4 4];
    build_gene_vs_protein_plots(datasets2, methods(2:end), dataset_labels2, ymaxs, DPI)
  
  
