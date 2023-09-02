
%%%%%%%%
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/tests')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/utils')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/wrappers')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/plotting')
%addpath('~/Documents/bioinfo_software/GEM_metabolic.models/RAVEN-main')

%initialize cobra
COBRA_SOLVER = 'gurobi';
changeCobraSolver(COBRA_SOLVER, 'all');

%%%%%%%%
file='log1p_female_body.csv'
dataset.name='fly.body'

%file='log1p_female_head.csv'
%dataset.name='fly.head'

table = importdata(file);
conditions = table.textdata(1, 2:end);
genes = table.textdata(2:end, 1);
transcriptomics = table.data;
transcriptomics(transcriptomics < 1e-12) = 0;

dataset.conditions = conditions;
dataset.genes = genes;
dataset.transcriptomics = transcriptomics;

%%%%%%%%
%load model
load('./Fruitfly-GEM-main/model/Fruitfly-GEM.mat');
model=fruitflyGEM;
model.grRules{10}

% all input parameters of `evaluate_method` already specified
IMAT_LOWER_QUANTILE = 0.25;
IMAT_UPPER_QUANTILE = 0.75;
IMAT_EPS = 1;

N = length(dataset.conditions);
runtime_all = zeros(1,N);
fluxes_exp_all = cell(1,N);

for i = 1:length(conditions)   
%for i = [1,5,9] 
    condition=conditions{i};
    disp(condition)
    cond_idx = strcmp(condition, dataset.conditions);
    gene_exp = dataset.transcriptomics(:,cond_idx);
    
    %disp(sum(gene_exp==0)) %remove zero expressed genes in below
    model = model;
    gene_names = dataset.genes(gene_exp~=0);
    gene_exp = gene_exp(gene_exp~=0);
    gene_exp_sd = [];
    lower_threshold = quantile(gene_exp, IMAT_LOWER_QUANTILE);
    upper_threshold = quantile(gene_exp, IMAT_UPPER_QUANTILE);
    eps_param = IMAT_EPS;

    % in call_iMAT.m, gene_to_reaction_levels.m, replace `:` in gene names to `-`.
    tstart = tic();    
    fluxes = call_iMAT(model, gene_names, gene_exp, ...
             lower_threshold, upper_threshold,eps_param);           
    runtime = toc(tstart);

    fluxes_exp_all{i} = fluxes;
    runtime_all(i) = runtime;        
end

method='iMat';
experiment.name = sprintf('%s_%s', method, dataset.name);
experiment.conditions = dataset.conditions;
experiment.runtime_all = runtime_all;
experiment.fluxes_exp_all = fluxes_exp_all;
save(['results/' experiment.name '.mat'], 'experiment');


