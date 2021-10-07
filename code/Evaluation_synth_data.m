clear
libraries_header;
chromosomes_number=22;

%% dataset parameters
data_instance_identifier='default_nips';
method_instance_identifier='default_lambda_0.01_beta_0.1_eta_0.1_wpca_0';
data_name='rosmap';

root_result_directory = 'Results/Evaluation_synth_data';
mkdir(root_result_directory);

number_of_optimization_iterations = 10; % plus an initial random intialization

weighted_precision_in_iterations = zeros(number_of_optimization_iterations+1,1);
weighted_recall_in_iterations = zeros(number_of_optimization_iterations+1,1);

weighted_precision_in_iterations_random = zeros(number_of_optimization_iterations+1,1);
weighted_recall_in_iterations_random = zeros(number_of_optimization_iterations+1,1);

for i=1:1:length(chromosomes_number)
    
    disp('Eval for synth data <3')

    chr_number = chromosomes_number(i);
    
    load(strcat(strcat('data/synth_data/',data_name),'/','chr_',num2str(chr_number),'_',data_instance_identifier,'.mat'))    

    root_addr =  strcat('inference/synth_data/',data_name,'/');
    load(strcat(root_addr,'chr_',num2str(chr_number),'_',method_instance_identifier,'.mat'))

    L_current_dataset = L_current_dataset(1:562,:);

    
    corr_L=[];
    corr_W = [];

    
    K_current_dataset_est = QUIC('default', cov(L_current_dataset), 0.01, 1e-6, 2, 100);
   

    
    for j=1:1:number_of_optimization_iterations+1

        temp_L = cellfun(@corr,num2cell(L_learned_in_diff_iters{j},1),num2cell(L_current_dataset,1));

        corr_L = [corr_L temp_L];

        temp_W = cellfun(@corr,num2cell(full(W_learned_in_diff_iters{j}),1),num2cell(full(W_current_dataset),1));

        corr_W = [corr_W temp_W];


        Target_linearized = reshape(K_current_dataset_est,[],1);
        Prediction_linearized = reshape(K_learned_in_diff_iters{j},[],1);
        
        Target_linearized(1:number_of_genes+1:end)=[];
        Prediction_linearized(1:number_of_genes+1:end)=[];


        % remove diagonal elements
        
        
         indx = randperm(length(Target_linearized));
         
         Target_linearized = abs(Target_linearized);
         Prediction_linearized = abs(Prediction_linearized);
         Target_linearized(Target_linearized~=0)=1;
         Prediction_linearized(Prediction_linearized~=0)=1;
        
         TP = sum((Target_linearized~=0 & Prediction_linearized~=0).*abs(Prediction_linearized));
         FP = sum((Target_linearized==0 & Prediction_linearized~=0).*abs(Prediction_linearized));         
         FN = sum((Target_linearized~=0 & Prediction_linearized==0).*abs(Target_linearized));
         TP_for_recall = sum((Target_linearized~=0 & Prediction_linearized~=0).*abs(Target_linearized));       
         weighted_precision_in_iterations(j)=TP/(TP+FP);
         weighted_recall_in_iterations(j)=TP_for_recall/(TP_for_recall+FN);

         Target_linearized_rnd = Target_linearized(indx);
         TP_rnd = sum((Target_linearized_rnd~=0 & Prediction_linearized~=0).*abs(Prediction_linearized));
         FP_rnd = sum((Target_linearized_rnd==0 & Prediction_linearized~=0).*abs(Prediction_linearized));
         FN_rnd = sum((Target_linearized_rnd~=0 & Prediction_linearized==0).*abs(Target_linearized_rnd));
         TP_for_recall_rnd = sum((Target_linearized_rnd~=0 & Prediction_linearized~=0).*abs(Target_linearized_rnd));        
         weighted_precision_in_iterations_random(j) = TP_rnd/(TP_rnd+FP_rnd);
         weighted_recall_in_iterations_random(j) = TP_for_recall_rnd/(TP_for_recall_rnd+FN_rnd);
        
        
    end

    
    number_of_probes_for_each_gene=sum(Gamma_binary,2);

    boxplot(abs(corr_L),ceil([1:1:length(corr_L)]./number_of_genes))
    xlabel('Coordinate Ascend Iteration');
    ylabel('Correlation(Learned-L,Ground Truth)');
    hold on
    plot(nanmean(abs(reshape(corr_L,[],number_of_optimization_iterations+1))),'g--*')
    hold off
    print(strcat(root_result_directory,'/','L-Correlation-BoxPlot-Synth-Chr',num2str(chr_number)),'-dpdf','-fillpage')

    
    boxplot(abs(corr_W),ceil([1:1:length(corr_W)]./number_of_genes))
    xlabel('Coordinate Ascend Iteration');
    ylabel('Correlation(Learned-W,Ground Truth)');
    hold on        
    plot(nanmean(abs(reshape(corr_W,[],number_of_optimization_iterations+1))),'g--*')
    legend('Mean')
    hold off
    print(strcat(root_result_directory,'/','W-Correlation-BoxPlot-Synth-Chr',num2str(chr_number)),'-dpdf','-fillpage')

    
end