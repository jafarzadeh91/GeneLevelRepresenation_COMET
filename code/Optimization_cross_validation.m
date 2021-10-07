lambda_arr = [0.01];

for lambda_ind = 1:1:length(lambda_arr)
    for beta_ind = 1:1:length(beta_arr)
        lambda = lambda_arr(lambda_ind);
        Optimization;
    end
end
