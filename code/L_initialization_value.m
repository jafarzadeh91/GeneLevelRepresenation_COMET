function L_initial = L_initialization_value(Gamma_determinestic_current_dataset,X_current_dataset,weighted)
%L_INITIALIZATION_VALUE Summary of this function goes here
%   Detailed explanation goes here
    number_of_genes = size(Gamma_determinestic_current_dataset,2);
    number_of_subjects= size(X_current_dataset,1);
    L_initial = zeros(number_of_subjects,number_of_genes);
    for gene_ind=1:1:number_of_genes
       probes_in_the_vicinity = find(Gamma_determinestic_current_dataset(:,gene_ind)~=0);
       distances_in_the_vicinity = Gamma_determinestic_current_dataset(probes_in_the_vicinity,gene_ind);
       distances_in_the_vicinity=abs(distances_in_the_vicinity)/sum(abs(distances_in_the_vicinity));
       X_local = X_current_dataset(:,probes_in_the_vicinity);
       X_local=full(transpose(X_local));
       
       if(weighted==true)
            X_local = full(diag(distances_in_the_vicinity)*X_local);
            X_local = X_local - repmat(mean(X_local),size(X_local,1),1);
       end
       if(size(X_local,1)==1)
            L_initial(:,gene_ind) = X_local';
       else
            wcoeff = pca(X_local,'VariableWeights','variance');
            coefforth = inv(diag(std(X_local)))* wcoeff;
            L_initial(:,gene_ind) = coefforth(:,1);
       end
    end 
end
