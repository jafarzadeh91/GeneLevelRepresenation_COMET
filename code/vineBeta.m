%% VineBeta.m
%% this function generate diverse types of PSD matrices with big non-diagonal values. The elements are forming a "probably" normal distribution around 0 and with some variances.
% @ Inputs:
% d : the size of matrix
% betaparam : effect the variance of afforsaid gaussian in some way.(Play with that.)
function S = vineBeta(d, betaparam)
    P = zeros(d);           %// storing partial correlations
    S = eye(d);

    for k = 1:d-1
        for i = k+1:d
            P(k,i) = betarnd(betaparam,betaparam); %// sampling from beta
            P(k,i) = (P(k,i)-0.5)*2;     %// linearly shifting to [-1, 1]
            p = P(k,i);
            for l = (k-1):-1:1 %// converting partial correlation to raw correlation
                p = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
            end
            S(k,i) = p;
            S(i,k) = p;
        end
    end

    %// permuting the variables to make the distribution permutation-invariant
    permutation = randperm(d);
    S = S(permutation, permutation);
end