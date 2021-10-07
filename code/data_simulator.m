%% data_simulator.m
% this function takes the Z matrix(d*q), sigma(q), n, and q 
% to generate synthetic data according to the SLGGM model.
% @ Inputs:
% Z : a binary d*q matrix where d is the number of probes and q is the number of
% genes. Each row of Z includes a 1 value and the remainings are zeros.
% sigma : a vector (of size q) that each element showes the variance for
% the distribution X_d' | L_q(d'), sigma_q(d') for a given probe d'.
% n : the number of subjects
% q : the number of genes
% mean_of_K: the mean of the elements of K matrix, according to our paper it's 0.
% variance_of_K: the variance of the elements of K matrix, which is
% equivalent to lambda in our paper.
% variance_of_K:
% trace : print what is happenning in console or not.
% @ Output:
% X : a n*d matrix showing the expression of probes for each subject
% K : a q*q matrix showing the inverse covariance matrix of different genes
% L : a n*q matrix showing the mean value for each gene in X|L, sigma*I
% distribution
function [X, K, L] = data_simulator(design_matrix, Sigma, n, q,W,K)
    
    %K= sprandsym(q,0.05,0.5);
    %inv_K=vineBeta(q,5);
    inv_K=inv(K);
    
    
    d = size(design_matrix,1);
    
    X = sparse(n,d);
    %X=randn(n,d);
    
    L = mvnrnd(zeros(n, q), inv_K);
    
    %L=X*W;
    
    
        
    for i=1:1:d
        disp(strcat(num2str(i),'/',num2str(d)))
       gene_index = find(design_matrix(i,:)~= 0);
       
       %if(slggm==1)
       %     X(:, i) = sum(mvnrnd(L(:, gene_index).*W(i,gene_index), Sigma(gene_index)),2);
       %elseif(nips==1)
        X(:, i) = sum(mvnrnd(full(L(:, gene_index).*W(i,gene_index)), Sigma(gene_index)),2);    
       %end
    end
    
end

