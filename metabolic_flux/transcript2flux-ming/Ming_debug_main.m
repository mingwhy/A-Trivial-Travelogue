%% MAIN FILE:
%
% Use this to run all the tests in the paper.
% Different tests can be run selectively.
% For convenience all results are stored as .mat files.
%
% Author: Daniel Machado, 2013

clc

addpath('tests')
addpath('utils')
addpath('wrappers')
addpath('plotting')
addpath('Lee12_funcs')
%addpath('/Users/ming/Downloads/Lee-12_bmc_matlab/yeast-GEM-main')
%ddpath('/Users/ming/Downloads/Lee-12_bmc_matlab/yeast-GEM-main/code')
addpath('/Users/ming/Documents/bioinfo_software/GEM_metabolic.models/RAVEN-main')

%COBRA_SOLVER = 'gurobi5';
COBRA_SOLVER = 'gurobi';
%initialize cobra
changeCobraSolver(COBRA_SOLVER, 'all');


%load('/Users/ming/Downloads/Lee-12_bmc_matlab/yeast-GEM-main/model/yeast-GEM.mat')

%% PROBLEM SETUP

% select dataset and organism (see load_dataset.m for details)
%DATASET = 'ishii'; ORGANISM = 'ecoli';
%DATASET = 'holm'; ORGANISM = 'ecoli';
DATASET = 'rintala'; ORGANISM = 'yeast';

% select experiment type
experiment_type = 'sim_all'; % simulate all fluxes 
%experiment_type = 'sim_intra'; % simulate intracellular fluxes 
%experiment_type = 'sim_secr'; % simulate secretion and growth rates


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


%% BENCHMARKING

%methods = {'pFBA', 'GIMME', 'iMAT', 'E-Flux', 'Lee-12', 'RELATCH', 'MADE', 'GX-FBA'};
%methods = {'iMAT', 'Lee-12'};
methods = {'Lee-12'};

%for i = 1:length(methods)
%    benchmark_method(methods{i}, model, dataset, options);
%end

%%%%%%%%%% 
%debug benchmark_method.m   
%function experiment = benchmark_method(method, model, dataset, options)
    
method=methods{1}; 

    N = length(dataset.conditions);
    
    status_all = zeros(1,N);
    error_all = zeros(1,N);
    runtime_all = zeros(1,N);
    fluxes_exp_all = cell(1,N);
    fluxes_sim_all = cell(1,N);
    
    options.precomputed = [];
    
    h = waitbar(0, sprintf('testing %s for %s dataset\n', method, dataset.name));
    set(findall(h,'type','text'),'Interpreter','none');

    %for i=1:N
    i=1;
        condition = dataset.conditions{i};
        ref_condition = dataset.conditions{1};
        %result = evaluate_method(model, dataset, method, condition, ref_condition, options);
%%%%%%%%%%         
%debug evaluate_method.m   
%function result = evaluate_method(model, dataset, method, condition, ref_condition, options)

    OBJ_FRAC = 0.9;
    GIMME_LOWER_QUANTILE = 0.25;
    IMAT_LOWER_QUANTILE = 0.25;
    IMAT_UPPER_QUANTILE = 0.75;
    IMAT_EPS = 1;    
 
        options.set_growth_rate = false;
        options.set_secretion_rates = false;
        options.ignore_internal_fluxes = false;
        
    
    cond_idx = strcmp(condition, dataset.conditions);
    ref_cond_idx = strcmp(ref_condition, dataset.conditions);

    gene_exp = dataset.transcriptomics(:,cond_idx);
    gene_exp_ref = dataset.transcriptomics(:,ref_cond_idx);
    
    fluxes_exp = dataset.fluxomics(:,cond_idx);
    fluxes_exp_ref = dataset.fluxomics(:,ref_cond_idx);

    %pick = ismember(dataset.reactions, model.rxns); %ming 2023-08-11
    %dataset.reactions = dataset.reactions(pick~=0); %ming 2023-08-11
    %fluxes_exp = fluxes_exp(pick~=0); %ming 2023-08-11
    %fluxes_exp_ref = fluxes_exp_ref(pick~=0); %ming 2023-08-11

    simulated = dataset.reactions;
    
    model_ref = model;
    
    biomass = find(model.c);
    biomass_rxn = model.rxns(biomass);

    

    fluxes_exp_old = fluxes_exp; %measured in `fit_fluxes_to_model`
    % fit_fluxes_to_model( model, reaction_ids, measured)    
    % reaction_ids=dataset.reactions
    % measured=fluxes_exp_old

    if options.reestimate_data
        fluxes_exp = fit_fluxes_to_model(model, dataset.reactions, fluxes_exp_old);
        fluxes_exp_ref = fit_fluxes_to_model(model_ref, dataset.reactions, fluxes_exp_ref);
    end
    
    if options.set_secretion_rates
        exchange_rxns = intersect(dataset.reactions, [model.uptk_rxns; model.secr_rxns]);
    else
        exchange_rxns = intersect(dataset.reactions, model.uptk_rxns);
    end

    [~, exchange_idxs_data] = ismember(exchange_rxns, dataset.reactions);
    exchange_rates = fluxes_exp(exchange_idxs_data);
    model  = constrain_model_to_data(model, exchange_rxns, exchange_rates);
    exchange_rates_ref = fluxes_exp_ref(exchange_idxs_data);
    model_ref = constrain_model_to_data(model_ref, exchange_rxns, exchange_rates_ref);
    
    if options.ignore_internal_fluxes
        simulated = intersect(dataset.reactions, [biomass_rxn; model.secr_rxns]);
    else
        simulated = setdiff(simulated, exchange_rxns);
    end

    [~, simulated_idx_data] = ismember(simulated, dataset.reactions);
    [~, simulated_idx_model] = ismember(simulated, model.rxns);
    [~, exchange_idxs_model] = ismember(exchange_rxns, model.rxns);
    
    fluxes_exp_sim = fluxes_exp(simulated_idx_data);

    uptake_rate = abs(fluxes_exp(strcmp(model.subtrate_uptake_rxn, dataset.reactions)));

    args.model = model;
    args.model_ref = model_ref;
    args.gene_names = dataset.genes;
    args.gene_exp = gene_exp;
    args.gene_exp_sd = [];
    args.gene_exp_ref = gene_exp_ref;

    args.external_rxns = exchange_rxns;
    args.external_rates = exchange_rates;
    args.external_rates_ref = exchange_rates_ref;
    args.obj_frac = OBJ_FRAC;

    args.gene_scale_rxn = model.subtrate_transport_rxn;
    args.flux_scale_rxn = model.subtrate_uptake_rxn;
    args.scale_value = uptake_rate;
    
    args.gimme.threshold = quantile(gene_exp, GIMME_LOWER_QUANTILE);
    args.FBA.knockouts = dataset.knockouts{cond_idx};
    args.iMAT.low_threshold = quantile(gene_exp, IMAT_LOWER_QUANTILE);
    args.iMAT.up_threshold = quantile(gene_exp, IMAT_UPPER_QUANTILE);
    args.iMAT.eps = IMAT_EPS;
    
    
% call_method.m
% case 'Lee-12'
    fluxes = call_Lee12(args.model, args.gene_names, args.gene_exp, ...
           args.gene_exp_sd, args.gene_scale_rxn, args.flux_scale_rxn, args.scale_value);
        
    [fluxes, status, runtime] = call_method(method, args, true);
    
    result.status = status;
    result.simulated = simulated;
    result.fluxes = fluxes;
    result.fluxes_exp = fluxes_exp;
    result.runtime = runtime;
    result.fluxes_exp_sim = fluxes_exp_sim;
    result.exchange_rxns = exchange_rxns;


