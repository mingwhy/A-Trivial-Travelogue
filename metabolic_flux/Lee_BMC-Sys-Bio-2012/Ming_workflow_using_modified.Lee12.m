
%% preparation
%# down load files from https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-6-73#Sec13
%# divide `12918_2012_943_MOESM6_ESM.m` (Additional file 6:Analysis) into separate **.m files.
%# rename files following `results.m`, then run matlab code in results.m
%# pretty slow, may take 1.5hr to estimate the flux vector per sample

%%%%%%%%%%%%%%

clc % clear would remove all stored objects in matlab

COBRA_SOLVER = 'gurobi';
%initialize cobra
changeCobraSolver(COBRA_SOLVER, 'all');

model_name      = 'yeast_5.21_MCISB.xml';
gene_to_scale   = 'glucose transport';  %reaction to scale gene expression
flux_to_scale   = 'D-glucose exchange'; %reaction to scale flux distribution

% used in estimation, make sure they exist, `call_Lee12.m`
scale_gene_idx = find(strcmp(gene_to_scale, model.rxnNames)) % or model.rxn
scale_rxn_idx = find(strcmp(flux_to_scale, model.rxnNames))


rsquared = @(f,y)(1 - sum((y-f).^2)/sum((y-mean(y)).^2));

%% 75%

%[reaction_name,experimental,p_gene_exp,p_standard_fba,p_standard_fba_best,p_gimme,p_shlomi] = ...
%    analysis(model_name, 'genedata_75.txt','experimental_fluxes_75.txt',gene_to_scale,flux_to_scale);

model_filename=model_name
genedata_filename='genedata_75.txt'
experimental_fluxes_filename='experimental_fluxes_75.txt'

% load model
model = readCbModel(model_name);
model = creategrRulesField(model); %https://opencobra.github.io/cobratoolbox/stable/modules/base/utilities/index.html
model.grRules{10}

% add `rev` back or not, see https://github.com/SysBioChalmers/RAVEN/issues/184
%Correct rev vector: true if LB < 0 & UB > 0, or it is an exchange reaction:
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
        model.rev(i) = true;
    end
end

% load transcript data
genedata	= importdata(genedata_filename);
gene_names	= genedata.textdata(:,1);
gene_names(1)= [];
gene_exp	= genedata.data(:,1);
gene_exp_sd	= genedata.data(:,2);


% map gene weighting to reaction weighting
[rxn_exp,rxn_exp_sd] = geneToReaction(model,gene_names,gene_exp,gene_exp_sd);

%%%%%%%%%%%%%%%%%%%%%%%
% sds 0 -> small
rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;

% scale by uptake reaction
uptake              = find(strcmp(gene_to_scale,model.rxnNames));
rxn_exp_sd          = rxn_exp_sd/rxn_exp(uptake);
rxn_exp             = rxn_exp/rxn_exp(uptake);
model.lb(uptake)    = 1;
model.ub(uptake)    = 1;

% FIX: scale all bounds accordingly, not just glc uptake 
% in `call_Lee12.m` from from `transcript2flux` https://github.com/cdanielmachado/transcript2flux
%model.lb = model.lb / abs(flux);
%model.ub = model.ub / abs(flux);

%%%%%%%%%%%%%%%%%%%%%%%
% Gene expression constraint FBA (init `COBRA_SOLVER = 'gurobi';` in advance)
% Lee BMC 2012, slow
v_gene_exp      = dataToFlux(model,rxn_exp,rxn_exp_sd);


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% Standard FBA, fast
solution        = optimizeCbModel(model,[],'one');
v_standard_fba  = solution.x;
%sum(v_standard_fba)

fOpt            = solution.f;
% "required metabolic functionalities" (=growth) set to 90% of maximum
model.lb(model.c == 1)  = 0.9*fOpt;
% gimme, fast
v_gimme         = gimme(model,rxn_exp);
%sum(v_gimme)

% schlomi iMat, 2008, slow
v_shlomi        = shlomi(model,rxn_exp);
% sum(v_shlomi)

%%%%%%%%%%%%%%%%%%%%%%%
% compare
experimental_fluxes = importdata(experimental_fluxes_filename);

reaction_name   = experimental_fluxes.textdata;
experimental    = zeros(size(experimental_fluxes.textdata,1),1);
p_gene_exp      = zeros(size(experimental_fluxes.textdata,1),1);
p_standard_fba  = p_gene_exp;
p_gimme         = p_gene_exp;
p_shlomi        = p_gene_exp;

flux = strcmp(flux_to_scale,reaction_name);
flux = experimental_fluxes.data(flux,1);

for k = 1:size(experimental_fluxes.textdata,1)
    j = find(strcmp(reaction_name{k},model.rxnNames));
    experimental(k)     = experimental_fluxes.data(k,1);
    p_gene_exp(k)       = flux*abs(v_gene_exp(j));
    p_standard_fba(k)   = flux*abs(v_standard_fba(j));
    p_gimme(k)          = flux*abs(v_gimme(j));
    p_shlomi(k)         = flux*abs(v_shlomi(j));
end

% remove small entries
p_gene_exp(abs(p_gene_exp)<1e-6)            = 0;
p_standard_fba(abs(p_standard_fba)<1e-6)    = 0;
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
p_standard_fba_best = zeros(size(experimental_fluxes.textdata,1),1);
for k = 1:size(experimental_fluxes.textdata,1)
    j = find(strcmp(experimental_fluxes.textdata{k,1},model.rxnNames));
    p_standard_fba_best(k) = flux*abs(v_standard_fba_best(j)); %#ok<FNDSB>
end
p_standard_fba_best(abs(p_standard_fba_best)<1e-6)  = 0;
