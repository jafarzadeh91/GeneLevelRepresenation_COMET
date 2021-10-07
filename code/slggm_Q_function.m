%% In this function we calculate the Q(theta, theta_old)as objective function of EM algorithm. FMI, go to Bishop P. 441
% @ Inputs:
% PI
% Gamma
% X
% L
% Sigma
% logger
% @ Outputs:
% func_value
function func_value = slggm_Q_function(PI, Gamma,W, X, L, Sigma)
    logger_active = 0;
    if(nargin==7)
        logger_active = 1;
    end
    
    if(logger_active==1)
        logger.trace(mfilename,'Started.')
    end
    
    log_PI=log(PI);
    func_value = 0;
    for row = 1:1:size(Gamma,1)
       cols = find(Gamma(row,:));
       row_res = Gamma(row, cols).*(log_PI(cols) + log_mvnpdf_with_numeric_sigma(X(:, row),L(:, cols).*W(row,cols), Sigma(cols)));
       func_value = func_value + sum(row_res);
    end
    
    if(logger_active==1)
        logger.trace(mfilename,'finished.')
    end
end