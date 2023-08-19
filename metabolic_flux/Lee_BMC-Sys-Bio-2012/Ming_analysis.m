% kieran: 26 apr 12
% modify original analysis.m code provided in Lee2012 supp files.

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
%[rxn_exp,rxn_exp_sd] = geneToReaction(model,genenames,gene_exp,gene_exp_sd);
%%%%%%%%%%%%%%%%%%%%%%%
m=model;
g=genenames;
t=gene_exp;
t_sd=gene_exp_sd;

% function [r,r_sd] = geneToReaction(m,g,t,t_sd) 
r       = zeros(size(m.rxns));
r_sd    = zeros(size(m.rxns));

for k = 1:length(g)
    g{k} = strrep(g{k},'-','_'); %rename gene names
end

for k = 1:length(m.rxns)
    ga = m.grRules{k};
    %ga = m.rules{k};
    ga = strrep(ga,'-','_');  %consistent with above, rename gene names
    w = regexp(ga,'\<\w*\>','match'); 
    w = setdiff(w,{'and','or'});
    for kk = 1:length(w)
        j = find(strcmp(w{kk},g)); % gene index if m.rules
        %j = w{kk};
        %if j=='x'
        %    continue
        %end
        n = t(j);
        n_sd = t_sd(j);        
        ga = regexprep(ga,['\<',w{kk},'\>'],[num2str(n),'±',num2str(n_sd)]);
    end
    
    %edit addGeneData
    %[n,n_sd] = addGeneData(ga);
    %%%%
    %[n,n_sd] = addGeneData(ga);
    gaa=ga;
    n = nan;
    n_sd = nan;

    ApmB = '[0-9\.]+±[0-9\.]+';

    if ~isempty(gaa)
        while isnan(n)        
            try         
                match_expr      = '([0-9\.])+±([0-9\.]+)';
                g_av            = regexprep(gaa,match_expr,'$1');
                g_sd            = regexprep(gaa,match_expr,'$2');
                n               = eval(g_av);
                n_sd            = eval(g_sd);                                
            catch %#ok<CTCH>            
                % replace brackets
                match_expr      = ['\((',ApmB,')\)'];
                replace_expr    = '$1';
                gaa = regexprep(gaa,match_expr,replace_expr);
                
                % replace and
                match_expr      = ['(',ApmB,') and (',ApmB,')'];
                replace_expr    = '${AandB($1,$2)}';
                gaa = regexprep(gaa,match_expr,replace_expr,'once');
                
                % replace or
                match_expr      = ['(',ApmB,') or (',ApmB,')'];
                replace_expr    = '${AorB($1,$2)}';
                gaa = regexprep(gaa,match_expr,replace_expr,'once');
               
            end
        end
    end
    %%%%
    
    r(k) = n;
    r_sd(k) = n_sd;
end

rxn_exp=r;
rxn_exp_sd=r_sd;
%%%%%%%%%%%%%%%%%%%%%%%

% sds 0 -> small
rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;

% scale by uptake reaction
uptake              = find(strcmp(gene_to_scale,model.rxnNames));
rxn_exp_sd          = rxn_exp_sd/rxn_exp(uptake);
rxn_exp             = rxn_exp/rxn_exp(uptake);
model.lb(uptake)	= 1;
model.ub(uptake)	= 1;

% Gene expression constraint FBA (init `COBRA_SOLVER = 'gurobi';` in advance)
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

