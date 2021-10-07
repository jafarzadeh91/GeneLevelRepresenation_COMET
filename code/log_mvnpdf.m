%% this function calculates the logarithm of value for multivariate and single-variable gaussian distribution.
% @ Inputs
% X: a n*1 consists of n samples
% Mu: a n*1 vector showing the average of X random variable
% Sigma: a n*n matrix showing the covariance
% multi-variate  normal distribution with mean and variance
function val = log_mvnpdf(X, Mu, Sigma)

    logdet_value = 2 *sum(log(diag(chol(2*pi*Sigma))));
    val = -(0.5)*((X-Mu)'/Sigma)*(X-Mu) - (0.5)*(logdet_value);
end

