 clear;
libraries_header;
chromosomes_number=[22];

%% load real data to simulate data based on that!
data_name='rosmap'; %% or 'KoborDNAm'
method_name = 'nips'; % or 'slggm'

for i=1:1:length(chromosomes_number)
    chr_number = chromosomes_number(i);
    real_data_load
    filter_by_chromosome;    
    quality_control_probes_genes

    number_of_genes_current_dataset = length(genes);
    
    
    K=generate_sparse_psd_block_matrix(number_of_genes_current_dataset,0.5,0.5,floor(number_of_genes_current_dataset/30));

    %% dataset parameters
%     sigma_scaling_factor_arr =[0.1 0.01 0.1 0.01] ;
%     beta_arr=[1 1 5 5];

    sigma_scaling_factor_arr =0.01 ;
    beta_arr=2;

    
    for data_ind=1:1:length(sigma_scaling_factor_arr)
            sigma_scaling_factor=sigma_scaling_factor_arr(data_ind);
            beta_factor = beta_arr(data_ind);
            Data_generator
    end
end

