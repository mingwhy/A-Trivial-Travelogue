% kieran: 26 apr 12

function [reaction_name,experimental,p_gene_exp,p_standard_fba,p_standard_fba_best,...
    p_gimme,p_shlomi] = ...
    analysis(model_filename,genedata_filename,experimental_fluxes_filename,...
    gene_to_scale,flux_to_scale)

% load model
model       = readCbModel(model_filename);

% load transcript data
genedata	= importdata(genedata_filename);
genenames	= genedata.textdata(:,1);
genenames(1)= [];
gene_exp	= genedata.data(:,1);
gene_exp_sd	= genedata.data(:,2);

% map gene weighting to reaction weighting
[rxn_exp,rxn_exp_sd] = geneToReaction(model,genenames,gene_exp,gene_exp_sd);
% sds 0 -> small
rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;

% scale by uptake reaction
uptake              = find(strcmp(gene_to_scale,model.rxnNames));
rxn_exp_sd          = rxn_exp_sd/rxn_exp(uptake);
rxn_exp             = rxn_exp/rxn_exp(uptake);
model.lb(uptake)	= 1;
model.ub(uptake)	= 1;

% Gene expression constraint FBA
v_gene_exp      = dataToFlux(model,rxn_exp,rxn_exp_sd);

% Standard FBA
solution        = optimizeCbModel(model,[],'one');
v_standard_fba	= solution.x;
fOpt            = solution.f;

% "required metabolic functionalities" (=growth) set to 90% of maximum
model.lb(model.c == 1)	= 0.9*fOpt;
% gimme
v_gimme         = gimme(model,rxn_exp);

% schlomi
v_shlomi    	= shlomi(model,rxn_exp);

% compare
experimental_fluxes = importdata(experimental_fluxes_filename);

reaction_name   = experimental_fluxes.textdata;
experimental	= zeros(size(experimental_fluxes.textdata,1),1);
p_gene_exp      = zeros(size(experimental_fluxes.textdata,1),1);
p_standard_fba	= p_gene_exp;
p_gimme         = p_gene_exp;
p_shlomi     	= p_gene_exp;

flux = strcmp(flux_to_scale,reaction_name);
flux = experimental_fluxes.data(flux,1);

for k = 1:size(experimental_fluxes.textdata,1)
    j = find(strcmp(reaction_name{k},model.rxnNames));
    experimental(k)     = experimental_fluxes.data(k,1);
    p_gene_exp(k)       = flux*abs(v_gene_exp(j));
    p_standard_fba(k)	= flux*abs(v_standard_fba(j));
    p_gimme(k)          = flux*abs(v_gimme(j));
    p_shlomi(k)     	= flux*abs(v_shlomi(j));
end

% remove small entries
p_gene_exp(abs(p_gene_exp)<1e-6)            = 0;
p_standard_fba(abs(p_standard_fba)<1e-6)	= 0;
p_gimme(abs(p_gimme)<1e-6)                  = 0;
p_shlomi(abs(p_shlomi)<1e-6)                = 0;

% find best fit from standard FBA solution
% ... overkill for this problem, but reuses existing method
model.lb(model.c == 1) = fOpt;
data = nan(size(model.rxns));
data(uptake) = 1;
for k = 1:size(experimental_fluxes.textdata,1)
    j = find(strcmp(experimental_fluxes.textdata{k,1},model.rxnNames));
    data(j) = experimental_fluxes.data(k,1)/flux; %#ok<FNDSB>
end
v_standard_fba_best = dataToFlux(model,data,data);
p_standard_fba_best	= zeros(size(experimental_fluxes.textdata,1),1);
for k = 1:size(experimental_fluxes.textdata,1)
    j = find(strcmp(experimental_fluxes.textdata{k,1},model.rxnNames));
    p_standard_fba_best(k) = flux*abs(v_standard_fba_best(j)); %#ok<FNDSB>
end
p_standard_fba_best(abs(p_standard_fba_best)<1e-6)	= 0;

