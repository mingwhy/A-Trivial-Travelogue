%%% MAIN FILE:
%
% Use this to run all the tests in the paper.
% Different tests can be run selectively.
% For convenience all results are stored as .mat files.
%
% Author: Daniel Machado, 2013

addpath
addpath('tests')
addpath('utils')
addpath('wrappers')
addpath('plotting')
%addpath('/Users/ming/Documents/bioinfo_software/GEM_metabolic.models/RAVEN-main')


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


%load data set
dataset = load_dataset(DATASET);

options.experiment_type = experiment_type;
options.reestimate_data = true; % fit experimental data to model
%options.reestimate_data = false; 

%methods = {'pFBA', 'iMAT', 'E-Flux','Lee-12', 'GIMME','RELATCH', 'MADE', 'GX-FBA'};
methods = {'iMAT'}; %GIMME, MADE, RELATCH,GX-FBA not work
method=methods{1}

%benchmark_method(method, model, dataset, options);
% code in benchmark_method.m

    N = length(dataset.conditions);
    
    status_all = zeros(1,N);
    error_all = zeros(1,N);
    runtime_all = zeros(1,N);
    fluxes_exp_all = cell(1,N);
    fluxes_sim_all = cell(1,N);
    
    options.precomputed = [];
        
    %for i=1:N
    i=3;
        condition = dataset.conditions{i};
        ref_condition = dataset.conditions{1};
        %result = evaluate_method(model, dataset, method, condition, ref_condition, options);
        
        % all input parameters of `evaluate_method` already specified
        IMAT_LOWER_QUANTILE = 0.25;
        IMAT_UPPER_QUANTILE = 0.75;
        IMAT_EPS = 1;
    
    cond_idx = strcmp(condition, dataset.conditions);
    gene_exp = dataset.transcriptomics(:,cond_idx);
    
%%%%%%%%%%%%%%%%%%%%%%%
        % use evaluate_method code, in case of 'sim_all'
        options.set_growth_rate = false;
        options.set_secretion_rates = false;
        options.ignore_internal_fluxes = false;

    cond_idx = strcmp(condition, dataset.conditions);
    ref_cond_idx = strcmp(ref_condition, dataset.conditions);

    gene_exp = dataset.transcriptomics(:,cond_idx);
    gene_exp_ref = dataset.transcriptomics(:,ref_cond_idx);
    
    fluxes_exp = dataset.fluxomics(:,cond_idx);
    fluxes_exp_ref = dataset.fluxomics(:,ref_cond_idx);


    pick = ismember(dataset.reactions, model.rxns); %ming 2023-08-11
    pick(pick==0) % check all reactions are mapped
    if length(fluxes_exp)~=sum(pick)
        dataset.reactions = dataset.reactions(pick~=0); %ming 2023-08-11
        fluxes_exp = fluxes_exp(pick~=0); %ming 2023-08-11
        fluxes_exp_ref = fluxes_exp_ref(pick~=0); %ming 2023-08-11
    end


    simulated = dataset.reactions;    
    model_ref = model;
    fluxes_exp_old = fluxes_exp; 
    if options.reestimate_data
        %measured in `fit_fluxes_to_model`
        % fit_fluxes_to_model( model, reaction_ids, measured)    
        % reaction_ids=dataset.reactions
        % measured=fluxes_exp_old
        fluxes_exp = fit_fluxes_to_model(model, dataset.reactions, fluxes_exp_old);
        fluxes_exp_ref = fit_fluxes_to_model(model_ref, dataset.reactions, fluxes_exp_ref);
    end


%%%%%%%%%%%%%%%%%%%%%%%
    args.model = model;
    args.gene_names = dataset.genes;
    args.gene_exp = gene_exp;
    args.gene_exp_sd = [];
    args.iMAT.low_threshold = quantile(gene_exp, IMAT_LOWER_QUANTILE);
    args.iMAT.up_threshold = quantile(gene_exp, IMAT_UPPER_QUANTILE);
    args.iMAT.eps = IMAT_EPS;


% in `call_method.m`: [fluxes, status, runtime] = call_method(method, args, true);
% case 'iMAT'
    tstart = tic();    
    fluxes = call_iMAT(args.model, args.gene_names, args.gene_exp, ...
             args.iMAT.low_threshold, args.iMAT.up_threshold, args.iMAT.eps);
           
    runtime = toc(tstart)


