
function v_sol = gimme(model,gene_exp)

% cutoff at lower quartile
cutoff                  = quantile(gene_exp(~isnan(gene_exp)),0.25);

% set up big matrices
model.c = zeros(size(model.c));

for k = 1:length(gene_exp)
    if gene_exp(k) < cutoff
        c = cutoff - gene_exp(k);
        [n1,n2] = size(model.S);
        % v = v+ - v-;
        model.S(n1+1,k)    = 1; model.b(n1+1) = 0;
        model.S(n1+1,n2+1) = -1; model.lb(n2+1) = 0; model.ub(n2+1) = inf; model.c(n2+1) = -c;
        model.S(n1+1,n2+2) = 1; model.lb(n2+2) = 0; model.ub(n2+2) = inf; model.c(n2+2) = -c;
    end
end

solution        = optimizeCbModel(model,[],'one');
v_sol = solution.x(1:length(model.rxns));
