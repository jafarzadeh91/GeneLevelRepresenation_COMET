%% in this function we calculate P(theta) as it is needed in EM objective function.
% @ Inputs: 
% L
% K
% lambda
% logger
% Note : Sigma has no prior in our model
% @ outputs
% func_value


function func_value = slggm_log_P_theta_function(L, K, lambda, logger)
    logger_active = 0;
    if(nargin==4)
        logger_active=1;
    end
    func_value = log_matrixnormalpdf(L, 0, inv(K));
    if(logger_active==1)
        logger.trace(mfilename,'started.');
    end
    if(logger_active==1)
        logger.trace(mfilename,'finished.');
    end
    if(lambda==0)
       return; 
    end
    sudo_cov = 1/lambda;
    func_value = func_value + (-1*sum(sum(abs(K-0)/sudo_cov))-(size(K,1)^2)*log(2*sudo_cov));

end

