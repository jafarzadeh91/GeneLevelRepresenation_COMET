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
    quality_control_probes_genes

    [~,common_users_meth,common_users_exp]=intersect(prob_subject.m.id,expression_and_phenotype.data.projid);

    root_addr =  strcat('inference/real_data/',data_name,'/');
    load(strcat(root_addr,'chr_',num2str(chr_number),'_',method_instance_identifier,'.mat'));
    %common_users_meth=1:1:length(expression_and_phenotype.data.projid);
    
    
    method_instance_identifier='default';
    
    corr_baselines = zeros(4,length(number_of_genes));
    corr_baselines = zeros(5,length(number_of_genes));

    corr_slggm_in_iterations = zeros(length(L_learned_in_diff_iters),length(number_of_genes));
    %% construct Z and Gamma matrices; Z = d*q shows which gene each prob belongs to, Gamma is row-wise normalized version of Z
    [Gamma_determinestic] =  probe_gene_distance_matrix_to_proximity_distribution(probe_gene_distance_matrix, 30000);
    non_zero_Gammas = find(Gamma_determinestic~=0);
    for gene_ind=1:1:length(genes)
       disp(gene_ind);
       prob_subject_for_this_gene = prob_subject.m.data(Gamma_determinestic(:,gene_ind)~=0,:)';
       if(size(prob_subject_for_this_gene,2)<2)
           continue;
       end
       prob_subject_for_this_gene_mean = mean(prob_subject_for_this_gene,2);
       %prob_subject_for_this_gene_mean = prob_subject_for_this_gene_mean(common_users_meth);
       %prob_subject_for_this_gene_centerized = prob_subject_for_this_gene(common_users_meth) - repmat(prob_subject_for_this_gene_mean,1,size(prob_subject_for_this_gene(common_users_meth,:),2));
       
       prob_subject_for_this_gene_mean = prob_subject_for_this_gene_mean(common_users_meth);
       prob_subject_for_this_gene_centerized = prob_subject_for_this_gene(common_users_meth,:) - repmat(prob_subject_for_this_gene_mean,1,size(prob_subject_for_this_gene,2));
       
       
       [coeff,score,latent] = pca(prob_subject_for_this_gene_centerized'); 

        first_pc = coeff(:,1);
        [~,ind] = max(Gamma_determinestic(:,gene_ind));
        first_prob = prob_subject_for_this_gene(common_users_meth,1);
        aggregated_probes = prob_subject_for_this_gene(common_users_meth,:)*Gamma_determinestic(Gamma_determinestic(:,gene_ind)~=0,gene_ind);
        
        W_for_weighted_PCA = zeros(size(prob_subject_for_this_gene,2),size(prob_subject_for_this_gene,2));
        W_for_weighted_PCA(1:size(prob_subject_for_this_gene,2)+1:end) = Gamma_determinestic(Gamma_determinestic(:,gene_ind)~=0,gene_ind);
        W_for_weighted_PCA(1:size(prob_subject_for_this_gene,2)+1:end) = W_for_weighted_PCA(1:size(prob_subject_for_this_gene,2)+1:end) / sum(W_for_weighted_PCA(1:size(prob_subject_for_this_gene,2)+1:end));
        
        prob_subject_for_this_gene_weighted_mean=mean(prob_subject_for_this_gene * W_for_weighted_PCA,2);
        prob_subject_for_this_gene_weighted_mean = prob_subject_for_this_gene_weighted_mean(common_users_meth);
        prob_subject_for_this_gene_centerized_weighted = prob_subject_for_this_gene(common_users_meth,:) - repmat(prob_subject_for_this_gene_weighted_mean,1,size(prob_subject_for_this_gene,2));
        Weighted_covariance = prob_subject_for_this_gene_centerized_weighted*W_for_weighted_PCA*prob_subject_for_this_gene_centerized_weighted';
        [U,S,V] = svd(Weighted_covariance);
        first_pc_weighted = U(:,1);
        %prob_subject_for_this_gene_centerized_weighted = prob_subject_for_this_gene_centerized * W_for_weighted_PCA;
        %[coeff_w,score_w,latent_w] = pca(prob_subject_for_this_gene_centerized_weighted'); 
        %first_pc_w = coeff(:,1);

        
        
        %% correlation calculation
        corr_baselines(1,gene_ind) = corr(first_pc,expression_and_phenotype.data.data(common_users_exp,gene_ind));
        corr_baselines(2,gene_ind) = corr(first_pc_weighted,expression_and_phenotype.data.data(common_users_exp,gene_ind));

        corr_baselines(3,gene_ind) = corr(first_prob,expression_and_phenotype.data.data(common_users_exp,gene_ind));
        corr_baselines(4,gene_ind) = corr(aggregated_probes,expression_and_phenotype.data.data(common_users_exp,gene_ind));
        corr_baselines(5,gene_ind) = corr(prob_subject_for_this_gene_mean,expression_and_phenotype.data.data(common_users_exp,gene_ind));
        
          for iter=1:1:length(L_learned_in_diff_iters)
                %corr_slggm_in_iterations(iter,gene_ind) = corr(L_learned_in_diff_iters{iter}(common_users_meth,gene_ind),expression_and_phenotype.data.data(common_users_exp,genes_with_probe_index(gene_ind)));         
                corr_slggm_in_iterations(iter,gene_ind) = corr(L_learned_in_diff_iters{iter}(common_users_meth,gene_ind),expression_and_phenotype.data.data(common_users_exp,genes_with_probe_index(gene_ind)));         

          end  
    end
    
    %% plot
    plot(log(objective_function_in_diff_iters(:,1)),'b--*')
    xlabel('Coordinate Ascend Iteration')
    ylabel('Log Objective Function Value')
    print(strcat(root_result_directory,'/','Log-Objective-Function-Value-Chr',num2str(chr_number)),'-dpdf','-fillpage')
    
    %% plot
    
    x = -1:0.05:1;
    [counts_bl]=hist(corr_baselines',x);
    hold on
    %bar(x,counts_bl(:,1)','FaceAlpha',0.8);
    %bar(x,counts_bl(:,2)','FaceAlpha',0.8);
    %bar(x,counts_bl(:,3)','FaceAlpha',0.8);
    
    counts_slggm = hist(corr_slggm_in_iterations(11,:),x);
    
    %bar(x,counts_slggm','FaceAlpha',0.8);
    
    count_mtrx = [counts_bl';counts_slggm]';
    
    count_mtrx_normalized  = count_mtrx;
    for i=1:1:size(count_mtrx,2)
        count_mtrx_normalized(:,i) = count_mtrx_normalized(:,i) / number_of_genes;      
    end
    
    figure;
    h = bar(categorical(x),count_mtrx_normalized);
    l = cell(1,6);
    l{1}='First PC'; l{2}='First PC W.';l{3}='Closest Prb.'; l{4}='Aggregated Prb.';l{5}='Average Prb.'; l{6}='Our Method';    
    legend(h,l);
    
    count_mtrx_normalized_max = max(count_mtrx_normalized,[],2);
    count_mtrx_max = max(count_mtrx,[],2);

    for i=1:1:length(x)
        if(count_mtrx_max(i)~=0)
            text(i,count_mtrx_normalized_max(i),strcat(num2str(count_mtrx_max(i))),'Rotation',90,'FontWeight','bold')
        end    
    end
end