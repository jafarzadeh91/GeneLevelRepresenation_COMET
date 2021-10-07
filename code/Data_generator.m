positive_weights=0;

data_instance_identifier=strcat('default_nips_sigma_',num2str(sigma_scaling_factor),'_beta_',num2str(beta_factor),'_positive_weights_',num2str(positive_weights));

number_of_genes_current_dataset = length(genes);

disp(strcat('probes_are_simulating_for_datasetomosome:',num2str(chr_number)));

%% construct or simulate other parts of the model        
[Gamma_determinestic_current_dataset] = probe_gene_distance_matrix_to_proximity_distribution(probe_gene_distance_matrix, 30000);          
Z_current_dataset = random_Z_sampler_from_Gamma(Gamma_determinestic_current_dataset);
Gamma_binary=Gamma_determinestic_current_dataset~=0;
non_zeros=find(Gamma_binary~=0);

W_current_dataset = sparse(size(Gamma_determinestic_current_dataset,1),size(Gamma_determinestic_current_dataset,2));

if(strcmp(method_name,'nips'))
                 Sigma_current_dataset = ones(1, number_of_genes_current_dataset)*sigma_scaling_factor;

    %for j=1:1:length(non_zeros)
         %W_current_dataset(non_zeros(j))=laprnd(1,1,0,nthroot(2,2) / Gamma_determinestic_current_dataset(non_zeros(j)));
         %W_current_dataset(non_zeros(j)) = ((rand(length(non_zeros(j)),1) > 0.5)*2 - 1).*Gamma_determinestic_current_dataset(non_zeros(j));

         %W_current_dataset(non_zeros(j)) = ((rand(length(non_zeros(j)),1) > 0.5)*2 - 1).*Gamma_determinestic_current_dataset(non_zeros(j));
         %W_current_dataset(non_zeros) = ((rand(length(non_zeros(j)),1) > 0.5)*2 - 1)*randn(length(non_zeros),1)./beta_factor+ Gamma_determinestic_current_dataset(non_zeros);
         %_current_dataset(non_zeros) = randn(length(non_zeros),1)./(beta_factor*Gamma_determinestic_current_dataset(non_zeros));
         W_current_dataset(non_zeros) = ((rand(length(non_zeros),1) > 0.5)*2 - 1).*(Gamma_determinestic_current_dataset(non_zeros)/(beta_factor));

         if(positive_weights==1)
                W_current_dataset(non_zeros)=abs(W_current_dataset(non_zeros));
         end
    %end

    %for pr=1:1:number_of_probes
    %   W_current_dataset(pr,:)=W_current_dataset(pr,:)/sum(abs(W_current_dataset(pr,:))); 
    %end
    design_matrix = Gamma_binary;


elseif(strcmp(method_name,'slggm'))
    for j=1:1:length(non_zeros)
         W_current_dataset(non_zeros(j))=randn();
         Sigma_current_dataset = rand(1, number_of_genes_current_dataset)*sigma_scaling_factor;

    end
    design_matrix = Z_current_dataset;
end


[X_current_dataset, K_current_dataset, L_current_dataset]=data_simulator(design_matrix, Sigma_current_dataset, number_of_subjects, number_of_genes,W_current_dataset,K);

%% data save
root_addr=strcat('data/synth_data/',data_name);
mkdir(root_addr);
clear i j;    
clear expression_and_phenotype prob_subject probe_gene subject_probe_measurement_matrix;
save(strcat(root_addr,'/','chr_',num2str(chr_number),'_',data_instance_identifier,'.mat'),'-v7.3');

