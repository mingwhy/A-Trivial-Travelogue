
function [n,n_sd] = addGeneData(gaa)
% geneToReaction.m from `transcript2flux` https://github.com/cdanielmachado/transcript2flux
% kieran: 22 july 11

    n = nan;
    n_sd = nan;

    %ApmB = '[0-9\.]+pm[0-9\.]+';
    %FIX Daniel
    ApmB = '([0-9\.\-e])+pm([0-9\.\-e]+)';
    
    tries = 0;

    f_and = @(a,b)AandB(a,b);
    f_or = @(a,b)AorB(a,b);
    
    if ~isempty(gaa)
        while isnan(n)
            
            tries = tries+1;
            if tries > 1000
                fprintf(1, 'warning: stuck at loop evaluating %s\n', gaa);
                break
            end

            try 
                match_expr      = ApmB; %'([0-9\.])+pm([0-9\.]+)';
                g_av            = regexprep(gaa,match_expr,'$1');
                g_sd            = regexprep(gaa,match_expr,'$2');
                n               = eval(g_av);
                n_sd            = eval(g_sd);

            catch %#ok<CTCH>

                % replace brackets 
                match_expr      = ['\(\s*(',ApmB,')\s*\)']; % FIX 
                %match_expr      = ['\((',ApmB,')\)'];

                replace_expr    = '$1';
                gaa = regexprep(gaa,match_expr,replace_expr);

                % replace and
                match_expr      = ['(',ApmB,')\s+and\s+(',ApmB,')']; % FIX 
                %match_expr      = ['(',ApmB,') and (',ApmB,')'];
                replace_expr    = '${f_and($1,$2)}';
                %replace_expr    = '${AandB($1,$2)}';
                gaa = regexprep(gaa,match_expr,replace_expr,'once');

                % replace or
                match_expr      = ['(',ApmB,')\s+or\s+(',ApmB,')']; % FIX 
                %match_expr      = ['(',ApmB,') or (',ApmB,')'];
                replace_expr    = '${f_or($1,$2)}';
                %replace_expr    = '${AorB($1,$2)}';
                gaa = regexprep(gaa,match_expr,replace_expr,'once');

            end
        end
    end

end
