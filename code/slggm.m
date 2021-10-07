%% slggm.m
% In this function, we implement the SLGGM algorithm.
%% @ Inputs
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
%% @ Outputs
% L_learned_in_diff_iters
% Sigma_learned_in_diff_iters
% K_learned_in_diff_iters
% Gamma_learned_in_diff_iters
% objective_function_in_diff_iters
function [L_learned_in_diff_iters, Sigma_learned_in_diff_iters, K_learned_in_diff_iters, Gamma_learned_in_diff_iters,W_learned_in_diff_iters, objective_function_in_diff_iters] = slggm(X, Gamma_initial, W_initial, L_initial, Sigma_initial, K_initial, lambda, number_of_optimization_iteration, number_of_coordinate_ascent_iterations, optimizer_function)

            
    L_learned_in_diff_iters = cell(1+number_of_optimization_iteration, 1);
    Sigma_learned_in_diff_iters = cell(1+number_of_optimization_iteration, 1);
    K_learned_in_diff_iters = cell(1+number_of_optimization_iteration, 1);
    Gamma_learned_in_diff_iters = cell(1+number_of_optimization_iteration, 1);
    W_learned_in_diff_iters = cell(1+number_of_optimization_iteration, 1);

    objective_function_in_diff_iters = zeros(1+number_of_optimization_iteration, 3);


    
    L_learned = L_initial;
    Sigma_learned = Sigma_initial;
    K_learned = K_initial;
    Gamma_learned = Gamma_initial;
    W_learned=W_initial;
    
    pi_learned_unnormalized = sum(Gamma_learned, 1);
    pi_learned = pi_learned_unnormalized/sum(pi_learned_unnormalized);
    
    L_learned_in_diff_iters{1} = L_learned;
    Sigma_learned_in_diff_iters{1} = Sigma_learned;
    K_learned_in_diff_iters{1} = K_learned;
    Gamma_learned_in_diff_iters{1} = Gamma_learned;
    W_learned_in_diff_iters{1} = W_learned;
    
    number_of_genes = size(K_learned, 1);
    number_of_subjects = size(X, 1);


    objective_function_in_diff_iters(1,1) = slggm_Q_function(pi_learned, Gamma_learned,W_learned, X, L_learned, Sigma_learned);
    objective_function_in_diff_iters(1,2) = slggm_log_P_theta_function(L_learned, K_learned, lambda);
    objective_function_in_diff_iters(1,3) = slggm_log_P_theta_function(L_learned, K_learned, 0);
    

    binary_Gamma = Gamma_learned~=0;
   
    for iter_number = 1:1:number_of_optimization_iteration
        
        
       %% E step:
       Gamma_learned = slggm_E_step(pi_learned, cell2mat(X_for_different_chunks_of_subjects), W_learned, L_learned, Sigma_learned, binary_Gamma, logger);
 
       %% add prior to Gamma
       number_of_probes_for_each_gene=sum(Gamma_learned ~= 0, 1);
       if(~isempty(find(number_of_probes_for_each_gene==0, 1)))
          Gamma_learned(Gamma_initial~=0) = Gamma_learned(Gamma_initial~=0) + eps;
       end
       sum_vec_Gamma_learned = sum(Gamma_learned, 2);
       Gamma_learned = bsxfun(@rdivide,Gamma_learned,sum_vec_Gamma_learned(:)) ;
          
       %% update pi
       pi_learned_unnormalized = sum(Gamma_learned, 1);
       pi_learned = pi_learned_unnormalized/sum(pi_learned_unnormalized);
        
       %% M Step:
     
       for cor_asc_iter_ind = 1:1:number_of_coordinate_ascent_iterations

            %% update K            

            if(strcmp(optimizer_function,'bigquic'))
                K_learned = bigquic(L_learned, strcat('-l ',num2str(lambda),' -n ',num2str(number_of_threads)));
            elseif(strcmp(optimizer_function,'quic'))
                K_learned = QUIC('default', cov(L_learned), lambda, 1e-6, 2, 100);
            elseif(strcmp(optimizer_function,'bcd'))
                %K_learned = L1precisionBCD(cov(L_learned), lambda);
            else
               error('optimizer_function was not specified. ') 
            end         
        
            
            %% sigma
            
            for gene_number=1:1:number_of_genes

               mtrx = X(:,:)- L_learned(:,gene_number)*W_learned(:,gene_number)';
               mtrx=sum(mtrx.^2);

               Sigma_learned(gene_number) = sum(Gamma_learned(:,gene_number).*mtrx');
               Sigma_learned(gene_number) = Sigma_learned(gene_number) / (number_of_subjects*sum(Gamma_learned(:,gene_number)));
            end

                
            %% Learn W
            [rows,cols,~]=find(binary_Gamma~=0);
            for i=1:1:length(rows)               
                W_learned(rows(i),cols(i)) = L_learned(:,cols(i))\X(:,rows(i));
            end


            
            %% update L
            for gene_number=1:1:number_of_genes
                L_zero_col = L_learned;
                L_zero_col(:,gene_number)=0;
                L_learned(:, gene_number) = (X*(Gamma_learned(:,gene_number).*W_learned(:,gene_number)) - Sigma_learned(gene_number)*(L_zero_col * K_learned(:,gene_number)))/(sum(Gamma_learned(:,gene_number).*W_learned(:,gene_number).^2)+Sigma_learned(gene_number)*K_learned(gene_number,gene_number));      
            end

       end

        
       
        %% In this step we need to take care of genes with no associated probes showing up with Gamma(:,j) = all zeros , Sigma(j) = Inf and L(:,j) = NaN that completely makes sense btw. 
        L_learned_in_diff_iters{1+iter_number} = L_learned;
        Sigma_learned_in_diff_iters{1+iter_number} = Sigma_learned;
        K_learned_in_diff_iters{1+iter_number} = K_learned;
        Gamma_learned_in_diff_iters{1+iter_number} = Gamma_learned;
        W_learned_in_diff_iters{1+iter_number} = W_learned;
        
        objective_function_in_diff_iters(iter_number+1,1) = slggm_Q_function(pi_learned, Gamma_learned, X,W_learned, L_learned, Sigma_learned);
        objective_function_in_diff_iters(iter_number+1,2) = slggm_log_P_theta_function(L_learned, K_learned, lambda);
        objective_function_in_diff_iters(iter_number+1,3) = slggm_log_P_theta_function(L_learned, K_learned, 0);
    

    end
    
end

