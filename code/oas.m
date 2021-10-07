% Oracle Approximating Shrinkage (OAS)
% Input:    X = nxp observation matrix, n = #samples, p = #features
% Output:   C = pxp well-conditioned covariance matrix
%           S = pxp empirical covariance matrix
% Reference: Chen et al., Shrinkage Algorithms for MMSE Covariance Estimation, TSP, 2010
function [C,S] = oas(X)
    [n,p] = size(X);
    S = cov(X,1); % Empirical covariance
    F = trace(S)*eye(p)/p; % Most well-conditioned Estimate
    trS2 = trace(S*S);
    tr2S = trace(S)^2;
    rho = min(((1-2/p)*trS2+tr2S)/((n+1-2/p)*(trS2-tr2S/p)),1); % Relative weighting
    C = (1-rho)*S+rho*F;
end