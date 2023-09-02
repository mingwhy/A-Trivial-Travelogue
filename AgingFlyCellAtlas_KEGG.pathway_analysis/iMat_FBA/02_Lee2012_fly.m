
%%%%%%%%
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/tests')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/utils')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/wrappers')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/plotting')
addpath('~/Documents/bioinfo_software/GEM_metabolic.models/transcript2flux-ming/Lee12_funcs')
%addpath('~/Documents/bioinfo_software/GEM_metabolic.models/RAVEN-main')

%initialize cobra
COBRA_SOLVER = 'gurobi';
changeCobraSolver(COBRA_SOLVER, 'all');

%%%%%%%%
%dataset.name='fly.head'
%file='log1p_female_head.csv'
%file2='log1p_female_head_sd.csv'

dataset.name='fly.body'
file='log1p_female_body.csv'
file2='log1p_female_body_sd.csv'

table = importdata(file);
conditions = table.textdata(1, 2:end); %cell type name
genes = table.textdata(2:end, 1);
transcriptomics = table.data;
transcriptomics(transcriptomics < 1e-12) = 0;

table = importdata(file2);
transcriptomicsSD=table.data(:,2:end);
transcriptomics(transcriptomics < 1e-12) = 0;

length(genes)
size(transcriptomics)
size(transcriptomicsSD)
transcriptomicsSD(1:3,1:3)

dataset.conditions = conditions;
dataset.genes = genes;
dataset.transcriptomics = transcriptomics;
dataset.transcriptomicsSD = transcriptomicsSD;

%%%%%%%%
%load model
load('./Fruitfly-GEM-main/model/Fruitfly-GEM.mat');
model=fruitflyGEM;
model.grRules{10}
model.rev(10) % make rue model.rev exist or add them mannully

N = length(dataset.conditions);
runtime_all = zeros(1,N);
fluxes_exp_all = cell(1,N);

for i = 1:length(conditions)   
%for i = [1,5,9] 
    condition=conditions{i};
    disp(condition)
    cond_idx = strcmp(condition, dataset.conditions);
    gene_exp = dataset.transcriptomics(:,cond_idx);
    gene_exp_sd = dataset.transcriptomicsSD(:,cond_idx);
    
    %disp(sum(gene_exp==0)) %remove zero expressed genes in below
    model = model;
    gene_names = dataset.genes(gene_exp~=0);
    gene_exp = gene_exp(gene_exp~=0);
    gene_exp_sd = gene_exp_sd(gene_exp~=0);
    %gene_exp_sd = zeros(length(gene_exp));
   
    % map gene weighting to reaction weighting
    [rxn_exp,rxn_exp_sd] = geneToReaction(model,gene_names,gene_exp,gene_exp_sd);

    % sds 0 -> small
    rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;
   
    tstart = tic();    
    v_gene_exp      = dataToFlux(model,rxn_exp,rxn_exp_sd);   
    runtime = toc(tstart) 
    disp(runtime)

    fluxes_exp_all{i} = v_gene_exp;
    runtime_all(i) = runtime;        
end

method='Lee';
experiment.name = sprintf('%s_%s', method, dataset.name);
experiment.conditions = dataset.conditions;
experiment.runtime_all = runtime_all;
experiment.fluxes_exp_all = fluxes_exp_all;
save(['results/' experiment.name '.mat'], 'experiment');


