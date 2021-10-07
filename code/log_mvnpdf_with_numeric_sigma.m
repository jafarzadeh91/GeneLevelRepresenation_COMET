%% this function calculates the logarithm of value for multivariate and single-variable gaussian distribution.
% @ Inputs
% X: a n*1 consists of n samples
% Mu: a n*1 vector showing the average of X random variable
% Sigma: a n*n matrix showing the covariance
% multi-variate  normal distribution with mean and variance
function val = log_mvnpdf_with_numeric_sigma(X, Mu, sigma)
    sigma = sigma'; 
    logdet_value = size(Mu,1)*(log(sigma)+1.8379) ; %% log(2*pi)=1.8379
    %val = -(0.5)*((X-Mu)'/(sigma*eye(size(Mu,1))))*(X-Mu) - (0.5)*(logdet_value);
    val = -(0.5)*sum((X*ones(1,size(Mu,2))-Mu).^2)./sigma - (0.5)*(logdet_value);
    
end

