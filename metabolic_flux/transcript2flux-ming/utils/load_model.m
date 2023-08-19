function model = load_model(organism)
% Load cobra model for the selected organism.
% Options: 'ecoli', 'yeast'
%
% Author: Daniel Machado, 2013


    ECOLI_MODEL = 'models/iAF1260.xml'; 
    YEAST_MODEL = 'models/iTO977.xml';
    %model_name      = 'yeast_5.21_MCISB.xml'; %in Lee 2012
    %YEAST_MODEL = 'models/yeast_5.21_MCISB.xml'; %# down load files from https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-6-73#Sec13
    
    %organism=ORGANISM
    switch organism
        case 'ecoli'
            model = readCbModel(ECOLI_MODEL);
            model.lb(model.lb < 0) = -1000;
            model.lb(strcmp('EX_co2(e)', model.rxns)) = 0;
            [model.uptk_rxns, model.secr_rxns] = get_exchange_rxns_ecoli(model);
            model.subtrate_transport_rxn = 'GLCptspp';
            model.subtrate_uptake_rxn = 'EX_glc(e)';
            
        case 'yeast'
            model = readCbModel(YEAST_MODEL);
            model.c = strcmp('CBIOMASS', model.rxns);
            model.ub(model.ub > 1) = 1000;
            % uptk_rxns: all reactions names end with 'xtI'
            % secr_rxns: all reactions names end with 'xtO'
            [model.uptk_rxns, model.secr_rxns] = get_exchange_rxns_yeast(model); %see below
            model.subtrate_transport_rxn = 'GAL2_1'; % used in Lee 2012
            model.subtrate_uptake_rxn = 'GLCxtI'; % used in Lee 2012
            model.lb(strcmp('U214_',model.rxns)) = 0; 
                    
            % essential components for anaerobic conditions
            model.ub(strcmp('44DIMZYMSTxtI', model.rxns)) = 1000;
            model.ub(strcmp('C141xtI', model.rxns)) = 1000;
            model.ub(strcmp('C161xtI', model.rxns)) = 1000;
            model.ub(strcmp('C181xtI', model.rxns)) = 1000;
            model.ub(strcmp('ERG572224xtI', model.rxns)) = 1000;
            model.ub(strcmp('LANOSTxtI', model.rxns)) = 1000;
            model.ub(strcmp('ZYMSTxtI', model.rxns)) = 1000;
            
            % fix problems with lipid composition
            model = addExchangeRxn(model,{'m612'}, 0, 1e-6); 
            model = addExchangeRxn(model,{'m709'}, -1000, 0); 
    end
    
end

function [uptk_rxns, secr_rxns] = get_exchange_rxns_ecoli(model)
     uptk_rxns = model.rxns(strncmp('EX', model.rxns, 2) & model.lb < 0 );
     secr_rxns = model.rxns(strncmp('EX', model.rxns, 2) & model.lb == 0 );
end
            
function [uptk_rxns, secr_rxns] = get_exchange_rxns_yeast(model)
    %find(strcmp(model.rxns, '2MBACxtI')) %1287
    uptk_rxns = {};
    secr_rxns = {};
    for i = 1:length(model.rxns)
        rxn = model.rxns{i};
        if length(rxn) > 3 && strcmp(rxn(end-2:end), 'xtI') %reaction name end with 'xtI'
            uptk_rxns = [uptk_rxns; rxn]; %#ok<AGROW>
        elseif length(rxn) > 3 && strcmp(rxn(end-2:end), 'xtO')
            secr_rxns = [secr_rxns; rxn]; %#ok<AGROW>
        end
    end
end
