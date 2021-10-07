%% this function calculates the logarithm of value for matrix gaussian distribution.
% @ Inputs
% X: a n*p matrix consists of p samples of n-dimensional vectors.
function val = log_matrixnormalpdf(X, Mu, sigma)
    n = size(X,1);
    p = size(X,2);
    val = -(0.5)*trace(inv(sigma)*(X-Mu)'*speye(n)*(X-Mu)) - (n/2)*(p*log(2*pi)+logdet(sigma,'chol')) - (p/2);

end

