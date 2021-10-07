libraries_header;
chromosomes_number=22;
optimization_type = 'synth'; % or 'synth'
%% dataset parameters
data_instance_identifier='default_nips_sigma_0.01_beta_2';

%% optimization parameters
number_of_optimization_iterations = 2; % plus an initial random intialization
number_of_coordinate_ascent_iterations = 1;

lambda_arr=[0.02];
beta_arr_in_learning = [1 0 0.01 0.05 0.1];
eta_arr_in_learning=[0 0.01 0.05 0.1 0.5 1 2 4];

%wpca_state = [true, false];
data_name='rosmap';
method_name='nips'; %% 'slggm'


markers = {'.','+','*','x','square','diamond'};

%% cross-validation
objective_fuct_for_different_params=zeros(length(beta_arr_in_learning),length(eta_arr_in_learning));
for chr_ind=1:1:length(chromosomes_number)
    legends={};
    legend_ind=1;
    chr_number =  chromosomes_number(chr_ind);
    for lambda_ind=1:1:length(lambda_arr)
        for beta_ind=1:1:length(beta_arr_in_learning)
            for eta_ind=1:1:length(eta_arr_in_learning)
                    beta_current = beta_arr_in_learning(beta_ind);
                    eta_current = eta_arr_in_learning(eta_ind);
                    lambda = lambda_arr(lambda_ind);
                    method_instance_identifier=strcat('default_lambda_',num2str(lambda),'_beta_',num2str(beta_current),'_eta_',num2str(eta_current));

%                     load(strcat(strcat('data/synth_data/',data_name),'/','chr_',num2str(chr_number),'_',data_instance_identifier,'.mat'))    
%                     root_addr =  strcat('inference/synth_data/',data_name,'/',data_instance_identifier,'/');
%                     load(strcat(root_addr,'chr_',num2str(chr_number),'_',method_instance_identifier,'.mat'))

                    %% plot objective function
                    corr_L_train_in_iterations=zeros(number_of_optimization_iterations+1,number_of_genes);
                    corr_W_train_in_iterations=zeros(number_of_optimization_iterations+1,number_of_probes);

                    corr_L_test_in_iterations = zeros(number_of_optimization_iterations+1,number_of_genes);
                    
                    AUC_in_different_iterations = zeros(number_of_optimization_iterations+1,1);
                    linearized_gt_K = reshape(K_current_dataset,1,[]);
                    linearized_gt_K(1:357:end)=[];
                    linearized_gt_K (linearized_gt_K~=0)=1;
                    
                    for itr_ind=1:1:number_of_optimization_iterations+1
                        corr_W_train = zeros(number_of_probes,1);
                        corr_L_train = cellfun(@corr,num2cell(L_learned_in_diff_iters{itr_ind}(train_subjects_inds,:),1),num2cell(L_current_dataset(train_subjects_inds,:),1));
                        for prob_ind=1:1:number_of_probes
                            non_zeros_for_this_prb=find(Gamma_binary(prob_ind,:)~=0);
                            corr_W_train(prob_ind)= corr(W_learned_in_diff_iters{itr_ind}(prob_ind,non_zeros_for_this_prb)',W_current_dataset(prob_ind,non_zeros_for_this_prb)');
                        end
                         corr_L_test = cellfun(@corr,num2cell(L_learned_in_diff_iters_test{itr_ind},1),num2cell(L_current_dataset(test_subjects_inds,:),1));
%                         corr_L_train_in_iterations(itr_ind)=mean(abs(corr_L_train));
%                         corr_W_train_in_iterations(itr_ind)=mean(abs(corr_W_train));
%                         corr_L_test_in_iterations(itr_ind)=mean(abs(corr_L_test));
                        corr_L_train_in_iterations(itr_ind,:)=corr_L_train;
                        corr_W_train_in_iterations(itr_ind,:)=corr_W_train;
                        corr_L_test_in_iterations(itr_ind,:)=corr_L_test;
                        
                        linearized_learned_l = K_learned_in_diff_iters{itr_ind};
                        linearized_learned_l = reshape(linearized_learned_l,1,[]);
                        linearized_learned_l(1:357:end) =[];

                        [~,~,~,AUC]=perfcurve(linearized_gt_K,abs(linearized_learned_l),1);    
                        AUC_in_different_iterations(itr_ind)=AUC;
                    end


                    %% plots
                    figure
                    root_addr_res =  strcat('results/synth_data/',data_name,'/',data_instance_identifier,'/');
                    mkdir(root_addr_res)
                    hold on
                    boxplot(corr_L_train_in_iterations','Labels',{1:1:number_of_optimization_iterations+1})
                    plot(mean(corr_L_train_in_iterations'))
                    xlabel('Coordinate Descent Iteration')
                    ylabel('Correlation')
                    saveas(gcf,strcat(root_addr_res,method_instance_identifier,'_L_train_correlation.png'));
                    hold off
                    close
                    
                    figure
                    hold on
                    boxplot(abs(corr_L_test_in_iterations)','Labels',{1:1:number_of_optimization_iterations+1})
                    plot(mean(abs(corr_L_test_in_iterations)'))
                    xlabel('Coordinate Descent Iteration')
                    ylabel('Correlation')
                    saveas(gcf,strcat(root_addr_res,method_instance_identifier,'_L_test_correlation.png'));
                    hold off
                    close
                    
                    figure
                    hold on
                    boxplot(corr_W_train_in_iterations','Labels',{1:1:number_of_optimization_iterations+1})
                    plot(nanmean(corr_W_train_in_iterations'))
                    xlabel('Coordinate Descent Iteration')
                    ylabel('Correlation')
                    saveas(gcf,strcat(root_addr_res,method_instance_identifier,'_W_correlation.png'));
                    hold off
                    close
                    
                    figure
                    hold on
                    x_inds=[1 3:3:size(objective_function_in_diff_iters,1)];
                    plot(objective_function_in_diff_iters(x_inds,1)/number_of_probes)
                    plot(objective_function_in_diff_iters_test(:,1)/number_of_probes)
                    legend('Train Data', 'Test Data')
                    xlabel('Coordinate Descent Iteration')
                    ylabel('Mean Square Error (Obj. Func.)')
                    grid on
                    saveas(gcf,strcat(root_addr_res,method_instance_identifier,'_Objective_Function.png'));

                    hold off
                    close

                    figure
                    hold on
                    plot(AUC_in_different_iterations(2:end))
                    xlabel('Coordinate Descent Iteration')
                    ylabel('Area Under the Curve')
                    grid on
                    saveas(gcf,strcat(root_addr_res,method_instance_identifier,'_K_AUC.png'));

                    hold off
                    close
                    
                    objective_fuct_for_different_params(beta_ind,eta_ind) = min(objective_function_in_diff_iters_test(:,1));
                    %saveas(gcf,strcat('Results/',method_instance_identifier,'_L_corr.png'));
%                     plot(objective_function_in_diff_iters_test)
%                     plot(corr_L_train_in_iterations);
%                     plot(corr_L_test_in_iterations);
%                     plot(corr_W_train_in_iterations); 
                        
            end  
        end
    end
    
    objective_fuct_for_different_params = objective_fuct_for_different_params /number_of_probes;
    h=heatmap(objective_fuct_for_different_params);
    h.XDisplayLabels=beta_arr_in_learning;
    h.YDisplayLabels = eta_arr_in_learning;
    axp = struct(h);       %you will get a warning
    axp.Axes.XAxisLocation = 'top';
    legend(legends,'Interpreter','latex');
    xlabel('\beta')
    ylabel('\eta')
    hold off;
end
   
                
