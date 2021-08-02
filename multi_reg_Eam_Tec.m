
% Programmer: socexp, WZJ, rewrite 20210723 
% Input: Tec -- Theoretical predictor, eigencycle
%        Eam -- Abed experiment angular momentum 
% Output: multi linear regression results 
%          coef_coefint3col -- the coefficient, 95% confidence intervals 
%                    for the coefficient estimates. 6 raws 3 column
%          resid_residint3col -- the residual, 95% confidence intervals 
%                    for the residuals estimates. 6 rows 3 column
%          R2_F_p_evar -- the R2 statistic, the F statistic 
%                     and its p value, and an estimate 
%                     of the error variance. 1 row 4 column
function [coef_coefint3col,resid_residint3col,R2_F_p_evar] = ...
                       multi_reg_Eam_Tec(Eam,Tec)

    y=Eam;
    x=Tec;
    X = [ones(size(y)) x]; 
    [coef,coefint,resid,residint,R2_F_p_evar] = regress(y,X);
    coef_coefint3col = [coef coefint];
    resid_residint3col = [resid residint];
    
end