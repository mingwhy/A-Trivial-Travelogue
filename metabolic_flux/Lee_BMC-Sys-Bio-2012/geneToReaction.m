
function [r,r_sd] = geneToReaction(m,g,t,t_sd)
% geneToReaction.m from `transcript2flux` https://github.com/cdanielmachado/transcript2flux
%m=model;
%g=gene_names;
%t=gene_exp;
%t_sd=gene_exp_sd;

    % kieran: 16 sep 11

    r       = zeros(size(m.rxns));
    r_sd    = zeros(size(m.rxns));

    for k = 1:length(g)
        g{k} = strrep(g{k},'-','_');
    end

    for k = 1:length(m.rxns)
        %k=10
        ga = m.grRules{k};
        %ga = m.rules{k}; %ming 2023-08-11
        ga = strrep(ga,'-','_');
        w = regexp(ga,'\<\w*\>','match'); 
        w = setdiff(w,{'and','or'});
        for kk = 1:length(w)
            j = find(strcmp(w{kk},g));
            n = t(j);
            n_sd = t_sd(j);        
            ga = regexprep(ga,['\<',w{kk},'\>'],[num2str(n),'pm',num2str(n_sd)]);
            % ga = regexprep(ga,['\<',w{kk},'\>'],[num2str(n),'±',num2str(n_sd)]); #original Lee function
        end
        [n,n_sd] = addGeneData(ga);
        r(k) = n;
        r_sd(k) = n_sd;
    end

end

