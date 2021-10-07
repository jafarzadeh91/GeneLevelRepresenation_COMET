disp('start optimization.m...') 
libraries_header;
chromosomes_number=22:-1:1;


%optimization_type = 'real'; % or 'synth'
%% dataset parameters

%data_instance_identifier='default_nips_sigma_0.01_beta_2_positive_weights_0';
data_instance_identifier='default_nips';

%% optimization parameters
number_of_optimization_iterations = 20; % plus an initial random intialization
number_of_coordinate_ascent_iterations = 1;


lambda_arr=[0.01];
beta_arr_in_learning = [0 0.01 0.05 0.1 0.5 1 2 4];
eta_arr_in_learning=[eps 0.01 0.05 0.1 0.5 1 2 4];
%wpca_state = [true, false];
data_name='KoborDNAm';
method_name='nips'; %% 'slggm'

%% user selection
%real_data_load 

%for i=1:1:length(chromosomes_number)
        chr_number = chromosomes_number(i);
        disp('optimization started from SLGGM model <3') 
        real_data_load;
        filter_by_chromosome;    
        quality_control_probes_genes
        

        
        if(strcmp(optimization_type,'synth'))
            load(strcat(strcat('data/synth_data/',data_name),'/','chr_',num2str(chr_number),'_',data_instance_identifier,'.mat'))    
        elseif(strcmp(optimization_type,'real'))
            [Gamma_determinestic_current_dataset] =  probe_gene_distance_matrix_to_proximity_distribution(probe_gene_distance_matrix, 30000);

            X_current_dataset = subject_probe_measurement_matrix;

        end


    %% optimization

   
    Gamma_initial=Gamma_determinestic_current_dataset;
    Gamma_reversed = Gamma_initial;
    Gamma_reversed(Gamma_reversed~=0) = 1 ./ Gamma_reversed(Gamma_reversed~=0);
        
    K_initial=eye(number_of_genes);
    
    %% cross-validation
    for lambda_ind=1:1:length(lambda_arr)
        for beta_ind=1:1:length(beta_arr_in_learning)
            for eta_ind=1:1:length(eta_arr_in_learning)
 %               for wpca_state_ind=1:1:length(wpca_state)
 %                   wpca = wpca_state(wpca_state_ind);
                    beta = beta_arr_in_learning(beta_ind);
                    eta = eta_arr_in_learning(eta_ind);
                    lambda = lambda_arr(lambda_ind);
                    %method_instance_identifier=strcat('default_lambda_',num2str(lambda),'_beta_',num2str(beta),'_eta_',num2str(eta),'_wpca_',num2str(wpca));
                    method_instance_identifier=strcat('default_lambda_',num2str(lambda),'_beta_',num2str(beta),'_eta_',num2str(eta));


                %K_initial = QUIC('default', cov(L_initial), lambda, 1e-6, 2, 100);
                
%                 for gene_ind=1:1:length(genes)
%                     str_split = strsplit(genes{gene_ind},{':'});
%                     str_split = strsplit(str_split{2},{'.'});
%                     genes{gene_ind}=str_split{1};
%                 end
                
                    [train_subjects_inds, test_subjects_inds,~]=divideblock(number_of_subjects,0.8,0.2,0);

                    W_initial = Gamma_determinestic_current_dataset;
                    for pr_ind=1:1:number_of_probes
                       W_initial(pr_ind,:)=W_initial(pr_ind,:)./nthroot(sum(W_initial(pr_ind,:).^2),2);
                       
                        
                    end
                    %L_initial = (X_current_dataset(train_subjects_inds,:)*W_initial)/(W_initial'*W_initial);
                    L_initial = L_initialization_value(Gamma_determinestic_current_dataset,X_current_dataset(train_subjects_inds,:),true);
                    L_initial=zscore(L_initial);
                    
                    %W_initial = zeros(number_of_probes,number_of_genes);
                    %for prob_ind=1:1:number_of_probes
                    %    genes_indices=Gamma_initial(prob_ind,:)~=0;
                    %    x = lsqnonneg(L_initial(train_subjects_inds,genes_indices),X_current_dataset(train_subjects_inds,prob_ind));
                    %end
                    
%                  for probe_ind=1:1:number_of_probes
%                      genes_indices=Gamma_initial(probe_ind,:)~=0;
%                      Lambda_matrix= diag(Gamma_initial(probe_ind,genes_indices));
%                      W_initial(probe_ind, genes_indices) = (L_initial(train_subjects_inds,genes_indices)'*L_initial(train_subjects_inds,genes_indices)+beta*Lambda_matrix)^(-1)*L_initial(train_subjects_inds,genes_indices)'*X_current_dataset(train_subjects_inds,i);
%                  end
        

                    %% save optimization results
                    if(strcmp(optimization_type,'synth'))
                        root_addr =  strcat('inference/synth_data/',data_name,'/',data_instance_identifier,'/');
                    elseif(strcmp(optimization_type,'real'))
                        root_addr =  strcat('inference/real_data/',data_name,'/',data_instance_identifier,'/');
                    end
                    mkdir(root_addr)
                

                    [L_learned_in_diff_iters,  K_learned_in_diff_iters, W_learned_in_diff_iters, objective_function_in_diff_iters] = nips(X_current_dataset(train_subjects_inds,:), Gamma_initial,W_initial, L_initial,  K_initial, beta, eta, lambda, number_of_optimization_iterations, number_of_coordinate_ascent_iterations, 'quic');        
                    %load(strcat(root_addr,'chr_',num2str(chr_number),'_',method_instance_identifier,'.mat'));
                    
                    objective_function_in_diff_iters_test = zeros(number_of_optimization_iterations+1,6);
    %               
                    L_learned_in_diff_iters_test = cell(length(L_learned_in_diff_iters),1);
                    Gamma_reversed = Gamma_initial;
                    Gamma_reversed(Gamma_reversed~=0) = 1 ./ Gamma_reversed(Gamma_reversed~=0);
        
                    gamma_reversed_normaized = Gamma_reversed./repmat(sum(Gamma_reversed,2),1,size(Gamma_reversed,2));

                    for counter = 1:1:number_of_optimization_iterations+1

                            L_test = (X_current_dataset(test_subjects_inds,:)*W_learned_in_diff_iters{counter})/(W_learned_in_diff_iters{counter}'*W_learned_in_diff_iters{counter} );
                            %L_test = regress(X_current_dataset(test_subjects_inds,:),W_learned_in_diff_iters{counter});
                            
                            %(*W_learned_in_diff_iters{counter})/('*W_learned_in_diff_iters{counter} );

                            L_learned_in_diff_iters_test{counter} = L_test;
% 
%                             objective_function_in_diff_iters_test(counter,1) = mean(abs(cellfun(@corr,num2cell(X_current_dataset(test_subjects_inds,:),1),num2cell(L_test*W_learned_in_diff_iters{counter}',1))));
%                             corr_learned = inv(K_learned_in_diff_iters{counter});
%                             corr_L = corrcoef(L_test);
%                             objective_function_in_diff_iters_test(counter,2) = mean(abs(cellfun(@corr,num2cell(corr_learned,1),num2cell(corr_L,1))));
% 
                             objective_function_in_diff_iters_test(counter,1) = sum(sum((X_current_dataset(test_subjects_inds,:)-L_test*W_learned_in_diff_iters{counter}').^2))/length(test_subjects_inds);
                             objective_function_in_diff_iters_test(counter,2) = beta*sum(sum(gamma_reversed_normaized.*abs(W_learned_in_diff_iters{counter}.^2)));
                             objective_function_in_diff_iters_test(counter,3) = objective_function_in_diff_iters_test(counter,1) + objective_function_in_diff_iters_test(counter,2);

                             objective_function_in_diff_iters_test(counter,4) = trace(L_test*K_learned_in_diff_iters{counter}*L_test')-logdet(K_learned_in_diff_iters{counter});
                             objective_function_in_diff_iters_test(counter,5) = lambda*sum(sum(abs(K_learned_in_diff_iters{counter})));
                             objective_function_in_diff_iters_test(counter,6) = (eta/length(test_subjects_inds))*(objective_function_in_diff_iters_test(counter,4) + objective_function_in_diff_iters_test(counter,5));



    %                     objective_function_in_diff_iters_test(counter,1) = sum(sum((X_current_dataset(test_subjects_inds,:)-L_test*W_learned_in_diff_iters{counter}').^2));
    %                     objective_function_in_diff_iters_test(counter,2) = beta*sum(sum(Gamma_reversed.*abs(W_learned_in_diff_iters{counter}.^2)));
    %                     objective_function_in_diff_iters_test(counter,3) = objective_function_in_diff_iters_test(counter,1) + objective_function_in_diff_iters_test(counter,2);
    % 
    %                     objective_function_in_diff_iters_test(counter,4) = trace(L_learned_in_diff_iters{counter}*K_learned_in_diff_iters{counter}*L_learned_in_diff_iters{counter}')-logdet(K_learned_in_diff_iters{counter});
    %                     objective_function_in_diff_iters_test(counter,5) = lambda*sum(sum(K_learned_in_diff_iters{counter}));
    %                     objective_function_in_diff_iters_test(counter,6) = objective_function_in_diff_iters_test(counter,4) + objective_function_in_diff_iters_test(counter,5);
                    end


                


                    clear prob_subject expression_and_phenotype  probe_gene
                    save(strcat(root_addr,'chr_',num2str(chr_number),'_',method_instance_identifier,'.mat'));

                %end
            end
        end
    end
%end