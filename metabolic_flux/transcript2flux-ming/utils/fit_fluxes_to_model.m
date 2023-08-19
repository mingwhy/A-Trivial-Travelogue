function [ fluxes, distance, status, fluxes_all ] = fit_fluxes_to_model( model, reaction_ids, measured)
% Find the nearest flux distribution (euclidean distance) to the given
% measured rates that fit the solution space of the model.
%
% Author: Daniel Machado, 2013

    [~, idxs] = ismember(reaction_ids, model.rxns);

    %measured(idxs==0)=[];  %ming 2023-08-11
    %idxs(idxs==0)=[];  %ming 2023-08-11
    
    [m, n] = size(model.S);
    QPproblem.F = zeros(n, n);
    for i = 1:length(idxs)
        j = idxs(i);
        if j==0   %ming 2023-08-11
            continue
        end
        QPproblem.F(j,j) = 1;
    end
    QPproblem.c = zeros(n, 1);    
    QPproblem.c(idxs) = -measured;
    %QPproblem.c(idxs(idxs~=0)) = -measured(idxs~=0);  %ming 2023-08-11

    QPproblem.A = model.S;
    QPproblem.b = zeros(m,1);
    QPproblem.lb = model.lb(1:n);
    QPproblem.ub = model.ub(1:n);
    QPproblem.csense = repmat('E', 1, m);
    QPproblem.osense = 1;
    
    try
        solution = solveCobraQP(QPproblem);
        status = solution.stat;
    catch e
        disp(e)
        status = 0;
    end

    if status == 1
        fluxes_all = solution.full;
        fluxes = fluxes_all(idxs);
        distance = norm(fluxes - measured);
        %fluxes = fluxes_all(idxs(idxs~=0)); %ming 2023-08-11
        %distance = norm(fluxes(idxs~=0) - measured(idxs~=0)); %ming 2023-08-11
    else
        fluxes = [];
        distance = nan;
    end
end

