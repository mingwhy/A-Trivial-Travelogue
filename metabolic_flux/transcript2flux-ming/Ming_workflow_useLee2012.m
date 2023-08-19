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
addpath('/Users/ming/Documents/bioinfo_software/GEM_metabolic.models/RAVEN-main')


%COBRA_SOLVER = 'gurobi5';
COBRA_SOLVER = 'gurobi';
%initialize cobra
changeCobraSolver(COBRA_SOLVER, 'all');

%load('/Users/ming/Documents/bioinfo_software/GEM_metabolic.models/yeast-GEM-main/model/yeast-GEM.mat')

%% PROBLEM SETUP

% select dataset and organism (see load_dataset.m for details)
%DATASET = 'ishii'; ORGANISM = 'ecoli';
%DATASET = 'holm'; ORGANISM = 'ecoli';
DATASET = 'rintala'; ORGANISM = 'yeast';


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


% select experiment type
experiment_type = 'sim_all'; % simulate all fluxes 
%experiment_type = 'sim_intra'; % simulate intracellular fluxes 
%experiment_type = 'sim_secr'; % simulate secretion and growth rates

options.experiment_type = experiment_type;
options.reestimate_data = true; % fit experimental data to model


%% BENCHMARKING

%methods = {'pFBA', 'GIMME', 'iMAT', 'E-Flux', 'Lee-12', 'RELATCH', 'MADE', 'GX-FBA'};
%methods = {'iMAT', 'Lee-12'};
methods = {'Lee-12'};

method=methods{1};  
%benchmark_method(method, model, dataset, options);
% use benchmark_method code   
    dataset.conditions
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
    
    %biomass = find(model.c);
    %biomass_rxn = model.rxns(biomass);
    

    fluxes_exp_old = fluxes_exp; 
    if options.reestimate_data
        %measured in `fit_fluxes_to_model`
        % fit_fluxes_to_model( model, reaction_ids, measured)    
        % reaction_ids=dataset.reactions
        % measured=fluxes_exp_old
        fluxes_exp = fit_fluxes_to_model(model, dataset.reactions, fluxes_exp_old);
        fluxes_exp_ref = fit_fluxes_to_model(model_ref, dataset.reactions, fluxes_exp_ref);
    end
    
    %%%%%%%%%%%%%%%%%%%% can not omit below
    if options.set_secretion_rates % use both uptake and secret reactions or only one
        exchange_rxns = intersect(dataset.reactions, [model.uptk_rxns; model.secr_rxns]);
    else
        % model.uptk_rxns specified load_model.m %only use uptake reactions  
        exchange_rxns = intersect(dataset.reactions, model.uptk_rxns);
    end

    %constrain_model_to_data.m in folder `utils`, only affect identified exchange_rxns in model.lb and model.ub
    % identify the fluxes for scaling used in `call_Lee12.m` line 71-74
    % similar to Lee 2012 `analysis.m` original code
    % gene_to_scale   = 'glucose transport'; 
    % flux_to_scale   = 'D-glucose exchange';
    [~, exchange_idxs_data] = ismember(exchange_rxns, dataset.reactions);
    exchange_rates = fluxes_exp(exchange_idxs_data);
    model  = constrain_model_to_data(model, exchange_rxns, exchange_rates);
    
    exchange_rates_ref = fluxes_exp_ref(exchange_idxs_data);
    model_ref = constrain_model_to_data(model_ref, exchange_rxns, exchange_rates_ref);
    
    if options.ignore_internal_fluxes
        simulated = intersect(dataset.reactions, [biomass_rxn; model.secr_rxns]);
    else
        simulated = setdiff(simulated, exchange_rxns); %difference between the two
    end

    %[~, simulated_idx_data] = ismember(simulated, dataset.reactions);
    %[~, simulated_idx_model] = ismember(simulated, model.rxns);
    %[~, exchange_idxs_model] = ismember(exchange_rxns, model.rxns); 
    %fluxes_exp_sim = fluxes_exp(simulated_idx_data);
    %%%%%%%%%%%%%%%%%%%%

    % model.subtrate_uptake_rxn defined in load_model.m 
    uptake_rate = abs(fluxes_exp(strcmp(model.subtrate_uptake_rxn, dataset.reactions)));

    args.model = model;

    args.gene_names = dataset.genes;
    args.gene_exp = gene_exp;
    args.gene_exp_sd = [];

    args.gene_scale_rxn = model.subtrate_transport_rxn;
    args.flux_scale_rxn = model.subtrate_uptake_rxn; % same as above `uptake_rate`
    args.scale_value = uptake_rate;
        

% in `call_method.m`: [fluxes, status, runtime] = call_method(method, args, true);
% case 'Lee-12'
    tstart = tic();
    fluxes = call_Lee12(args.model, args.gene_names, args.gene_exp, ...
           args.gene_exp_sd, args.gene_scale_rxn, args.flux_scale_rxn, args.scale_value);
    
    runtime = toc(tstart)


