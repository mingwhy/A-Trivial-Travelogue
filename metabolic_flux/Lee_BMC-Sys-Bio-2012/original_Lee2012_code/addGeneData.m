
function [n,n_sd] = addGeneData(gaa)

% kieran: 22 july 11

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