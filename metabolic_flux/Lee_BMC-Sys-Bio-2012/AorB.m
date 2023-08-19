function str = AorB(str1,str2) %#ok<DEFNU>
% geneToReaction.m from `transcript2flux` https://github.com/cdanielmachado/transcript2flux
    %ApmB = '([0-9\.])+pm([0-9\.]+)';

    %FIX Daniel
    ApmB = '([0-9\.\-e])+pm([0-9\.\-e]+)';

    match_expr      = ApmB;
    m1              = eval(regexprep(str1,match_expr,'$1'));
    s1              = eval(regexprep(str1,match_expr,'$2'));
    m2              = eval(regexprep(str2,match_expr,'$1'));
    s2              = eval(regexprep(str2,match_expr,'$2'));

    m = m1 + m2;

    s = sqrt(s1^2 + s2^2);

    str = [num2str(m),'pm',num2str(s)];
end
