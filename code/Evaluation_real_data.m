chr_numbers=22;

chunck_size = 30000;


%% dataset parameters
method_instance_identifier='default_lambda_0.01_beta_0.1_eta_0.01';


data_name='rosmap';


root_result_directory = 'Results/Evaluation_real_data';
mkdir(root_result_directory);

for chr_ind = 1:1:length(chr_numbers)
    
    chr_number = chr_numbers(chr_ind);
    real_data_rosmap_header
    filter_by_chromosome
    quality_control_probes_genes

    [~,common_users_meth,common_users_exp]=intersect(prob_subject.m.id,expression_and_phenotype.data.projid);

    root_addr =  strcat('inference/real_data/',data_name,'/');
    load(strcat(root_addr,'chr_',num2str(chr_number),'_',method_instance_identifier,'.mat'));

    corr_slggm_in_iterations = zeros(length(L_learned_in_diff_iters),length(number_of_genes));
    L_current_dataset = expression_and_phenotype.data.data(common_users_exp,genes_with_probe_index);
    
    corr_matrix = zeros(number_of_optimization_iterations+1,number_of_genes);
    for iter = 1:1:number_of_optimization_iterations+1
       current_L_zscored = zscore(L_learned_in_diff_iters{iter});
       current_L_zscored = current_L_zscored(common_users_exp,:);
       corr_vec=cellfun(@corr,num2cell(current_L_zscored,1),num2cell(L_current_dataset,1));
       corr_matrix(iter,:)=corr_vec;
        
    end
    corr_matrix = abs(corr_matrix);
end


