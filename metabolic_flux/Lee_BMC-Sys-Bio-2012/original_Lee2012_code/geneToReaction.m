
function [r,r_sd] = geneToReaction(m,g,t,t_sd)

% kieran: 16 sep 11

%m=model;
%g=genenames;
%t=gene_exp;
%t_sd=gene_exp_sd;

r       = zeros(size(m.rxns));
r_sd    = zeros(size(m.rxns));

for k = 1:length(g)
    g{k} = strrep(g{k},'-','_');
end

for k = 1:length(m.rxns)
    ga = m.grRules{k};
    %ga = m.rules{k};
    ga = strrep(ga,'-','_');
    w = regexp(ga,'\<\w*\>','match'); 
    w = setdiff(w,{'and','or'});
    for kk = 1:length(w)
        j = find(strcmp(w{kk},g));
        n = t(j);
        n_sd = t_sd(j);        
        ga = regexprep(ga,['\<',w{kk},'\>'],[num2str(n),'±',num2str(n_sd)]);
    end

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

