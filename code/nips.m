%% slggm.m
% In this function, we implement the SLGGM algorithm.
% @ Inputs
% X: a n*d matrix that show the measurement over all subjects for each
% probe
% Gamma_initial: a d*q matrix
% L_initial: a n*q matrix whihc show the average of probes measurements for
% different genes.
% Sigma_initial: a q dimenstion vector showing the variance for the
% distribtuion of probes measurement for each gene.
% K_initial: a q*q inverse covariance of genes, showing the partial
% correlation between the probes measurement of different genes.
% lambda: for quic or bigquic
% number_of_optimization_iteration: As slggm is an iterative algorithm,
% this parameter showes the number of interations
% number_of_coordinate_ascent_iterations
% optimizer_function: bigquic, quic or oas.
% resources
% logger:
%% @ Outputs
% L_learned_in_diff_iters
% Sigma_learned_in_diff_iters
% K_learned_in_diff_iters
% objective_function_in_diff_iters
function [L_learned_in_diff_iters,  K_learned_in_diff_iters, W_learned_in_diff_iters, objective_function_in_diff_iters] = nips(X, Gamma_initial, W_initial, L_initial,  K_initial, beta,eta, lambda, number_of_optimization_iteration, number_of_coordinate_ascent_iterations, optimizer_function)

        
    L_learned_in_diff_iters = cell(1+number_of_optimization_iteration, 1);
    K_learned_in_diff_iters = cell(1+number_of_optimization_iteration, 1);
    W_learned_in_diff_iters = cell(1+number_of_optimization_iteration, 1);

    objective_function_in_diff_iters = zeros(1+3*number_of_optimization_iteration, 7);

    
    L_learned = L_initial;

    K_learned = K_initial;
    W_learned=W_initial;

    L_learned_in_diff_iters{1} = L_learned;
    K_learned_in_diff_iters{1} = K_learned;
    W_learned_in_diff_iters{1} = W_learned;
    
    number_of_subjects = size(X, 1);
    number_of_probes = size(Gamma_initial,1);

  
    Gamma_reversed = Gamma_initial;
    Gamma_reversed(Gamma_reversed~=0) = 1 ./ Gamma_reversed(Gamma_reversed~=0);
        
     gamma_reversed_normaized = Gamma_reversed./repmat(sum(Gamma_reversed,2),1,size(Gamma_reversed,2));


     objective_function_in_diff_iters(1,1) = sum(sum((X-L_learned*W_learned').^2))/number_of_subjects;
     objective_function_in_diff_iters(1,2) = sum(sum(gamma_reversed_normaized.*abs(W_learned.^2)));
     objective_function_in_diff_iters(1,3) = objective_function_in_diff_iters(1,1) + objective_function_in_diff_iters(1,2);

     objective_function_in_diff_iters(1,4) = trace(L_learned*K_learned*L_learned')-logdet(K_learned);
     objective_function_in_diff_iters(1,5) = lambda*sum(sum(abs(K_learned)));
     objective_function_in_diff_iters(1,6) = (objective_function_in_diff_iters(1,4) + objective_function_in_diff_iters(1,5));



    lambda_max=1;  
    counter = 1;
    

    for iter_number = 1:1:number_of_optimization_iteration
     
        counter=counter+1;

        %logger.trace('slggm.m','K update started.')
        if(strcmp(optimizer_function,'bigquic'))
            K_learned = bigquic(L_learned, strcat('-l ',num2str(lambda),' -n ',num2str(number_of_threads)));
        elseif(strcmp(optimizer_function,'quic'))
            %L_learned_normalized=L_learned/max(max(L_learned));
            %L_learned_normalized_zero_diag=L_learned_normalized'*L_learned_normalized;
            %L_learned_normalized_zero_diag(1:size(K_initial,1)+1:end)=0;
            %lambda_max=max(max(abs(L_learned_normalized_zero_diag)));
            K_learned = QUIC('default', cov(L_learned), lambda, 1e-6, 2, 100);
        elseif(strcmp(optimizer_function,'bcd'))
            K_learned = L1precisionBCD(cov(L_learned), lambda);
        else
           error('optimizer_function was not specified. ') 
        end   
            

             objective_function_in_diff_iters(counter,1) = sum(sum((X-L_learned*W_learned').^2))/number_of_subjects;
             objective_function_in_diff_iters(counter,2) = sum(sum(gamma_reversed_normaized.*abs(W_learned.^2)));
             objective_function_in_diff_iters(counter,3) = objective_function_in_diff_iters(counter,1) + objective_function_in_diff_iters(counter,2);
 
             objective_function_in_diff_iters(counter,4) = trace(L_learned*K_learned*L_learned')-logdet(K_learned);
             objective_function_in_diff_iters(counter,5) = lambda_max*lambda*sum(sum(abs(K_learned)));
             objective_function_in_diff_iters(counter,6) = (objective_function_in_diff_iters(counter,4) + objective_function_in_diff_iters(counter,5));
             
 


              
              %% Learn W
                
             for i=1:1:number_of_probes
                 disp(i); 
                 genes_indices=Gamma_initial(i,:)~=0;
                 Lambda_matrix= diag(gamma_reversed_normaized(i,genes_indices));

                 %[lasso_fit, fit_info]=lasso(L_learned(:,genes_indices)*Lambda_matrix_inverse,X(:,i),'LambdaRatio',0,'cv',3,'NumLambda',10);
                 %W_learned(i,genes_indices)=lasso_fit(:,find(fit_info.MSE==min(fit_info.MSE),1))'*Lambda_matrix_inverse;
                 %L_corr=corrcoef(L_learned(:,genes_indices));
                 %L_corr = L_learned(:,genes_indices)'*L_learned(:,genes_indices);
                 %L_corr = L_corr /max(max(abs(L_corr)));
                
                 f= @(w_vec)(((X(:,i)-L_learned(:,genes_indices)*w_vec')'*(X(:,i)-L_learned(:,genes_indices)*w_vec')/number_of_subjects)+beta*length(genes_indices)*sum((gamma_reversed_normaized(i,genes_indices)*abs(w_vec)')));
                 %W_learned(i, genes_indices) = fmincon(f,W_learned(i,genes_indices),[],[],ones(1,length(find(genes_indices~=0))),1,zeros(length(find(genes_indices)),1));
                 %optimoptions('fmincon','SpecifyConstraintGradient',true);

                 W_learned(i, genes_indices) = fmincon(f,W_learned(i,genes_indices),[],[],[],[],[],[],@L2Norm);

                 %W_learned(i, genes_indices) = W_learned(i, genes_indices) /sum(W_learned(i, genes_indices));
                 %W_learned(i, genes_indices) = (L_corr+beta*Lambda_matrix)^(-1)*L_learned(:,genes_indices)'*X(:,i);
             end
             
            
             %objective_function_in_diff_iters(counter,1) = mean(abs(cellfun(@corr,num2cell(X,1),num2cell(L_learned*W_learned',1))));
             %corr_learned = inv(K_learned);
             %corr_L = corrcoef(L_learned);
             %objective_function_in_diff_iters(counter,2) = mean(abs(cellfun(@corr,num2cell(corr_learned,1),num2cell(corr_L,1))));
                
             counter=counter+1;
             
             objective_function_in_diff_iters(counter,1) = sum(sum((X-L_learned*W_learned').^2))/number_of_subjects;
             objective_function_in_diff_iters(counter,2) = sum(sum(gamma_reversed_normaized.*abs(W_learned.^2)));
             objective_function_in_diff_iters(counter,3) = objective_function_in_diff_iters(counter,1) + objective_function_in_diff_iters(counter,2);
 
             objective_function_in_diff_iters(counter,4) = trace(L_learned*K_learned*L_learned')-logdet(K_learned);
             objective_function_in_diff_iters(counter,5) = lambda_max*lambda*sum(sum(abs(K_learned)));
             objective_function_in_diff_iters(counter,6) = (objective_function_in_diff_iters(counter,4) + objective_function_in_diff_iters(counter,5));
             
 
            

    

            
            %% update L

            %nominator = zeros(size(L_learned,1),size(L_learned,2));
            %denominator=zeros(size(K_learned,1),size(K_learned,2));
            %for i=1:1:size(W_learned,1)
            %    disp(i)
            %   nominator = nominator + X(:,i)*W_learned(i,:);
            %   denominator = denominator + W_learned(i,:)'*W_learned(i,:);
            %end
            
            %denominator = denominator + inv(K_learned);
            
            %L_learned = nominator / denominator;
            %W_corr=corr(W_learned);
            %W_corr=W_corr/max(max(abs(W_corr)));
            %K_learned_normalized = K_learned/max(max(abs(K_learned)));

            L_learned = (X*W_learned)/(W_learned'*W_learned + eta*K_learned);
            %L_learned = (X*W_learned)/(W_learned'*W_learned );
            %L_learned = X*W_learned;

        
         counter = counter +1;
    
             objective_function_in_diff_iters(counter,1) = sum(sum((X-L_learned*W_learned').^2))/number_of_subjects;
             objective_function_in_diff_iters(counter,2) = beta*sum(sum(gamma_reversed_normaized.*abs(W_learned.^2)));
             objective_function_in_diff_iters(counter,3) = objective_function_in_diff_iters(counter,1) + objective_function_in_diff_iters(counter,2);
 
             objective_function_in_diff_iters(counter,4) = trace(L_learned*K_learned*L_learned')-logdet(K_learned);
             objective_function_in_diff_iters(counter,5) = lambda_max*lambda*sum(sum(abs(K_learned)));
             objective_function_in_diff_iters(counter,6) = (eta/number_of_subjects)*(objective_function_in_diff_iters(counter,4) + objective_function_in_diff_iters(counter,5));
             
 
        %% In this step we need to take care of genes with no associated probes showing up with Gamma(:,j) = all zeros , Sigma(j) = Inf and L(:,j) = NaN that completely makes sense btw. 
        L_learned_in_diff_iters{1+iter_number} = L_learned;
        K_learned_in_diff_iters{1+iter_number} = K_learned;
        W_learned_in_diff_iters{1+iter_number} = W_learned;

    end
    
end

