%% In this function we calculate the Gamma matrix which is the P(Z|X,theta) and will be used in M-step of EM algorithm.
% @ inputs:
% PI:prior probability of assignment of a prob to each gene (sudo-counting based)
% X: n*d matrix of subect-probe measurements.
% L: n*q matrx that each column shows the mean of probe measurements for
% each gene
% Sigma: The variance associated with distribution of X.
% binary_Gamma
% logger
% @ Outputs
% Gamma
function Gamma = slggm_E_step(PI, X, W, L, Sigma, binary_Gamma, logger)
    %SLGGM_EXPECTATION_CALCULATION Summary of this function goes here
    %Detailed explanation goes here
    logger_active=0;
    if(nargin==6)
        logger_active = 1;
    end
    
    Gamma = sparse(size(binary_Gamma,1), size(binary_Gamma,2));
    %[row, col] = find(binary_Gamma);
    log_PI=log(PI);
    if(logger_active==1)
        logger.trace(mfilename,'updating each row of gamma matrix...');
    end

    for row = 1:1:size(Gamma,1)
       cols = find(binary_Gamma(row,:));
       Gamma(row, cols) =  log_PI(cols) + log_mvnpdf_with_numeric_sigma(X(:, row),W(row,cols).*L(:, cols), Sigma(cols));
    end
    if(logger_active==1)
        logger.trace(mfilename,'normalizing rows of Gamma matrix...');
    end
     
    max_possible_log=log(realmax)*0.9;

    for i=1:1:size(Gamma, 1)
       non_zeros = find(Gamma(i,:)); 
       Gamma_i_non_zeros = Gamma(i, non_zeros);
       max_val = max(Gamma_i_non_zeros);
       subtracted_value = max_val-max_possible_log; 
       Gamma_i_non_zeros = Gamma_i_non_zeros - subtracted_value;
       %Gamma_i_non_zeros(Gamma_i_non_zeros<min_possible_log) = min_possible_log;
       Gamma_i_non_zeros = exp(Gamma_i_non_zeros);    
       Gamma(i,non_zeros)=Gamma_i_non_zeros;
    end
    
    sum_vec_Gamme = sum(Gamma, 2);
    Gamma = bsxfun(@rdivide,Gamma,sum_vec_Gamme(:));


end

