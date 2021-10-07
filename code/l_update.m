function L_learned = l_update(X, Gamma_learned, L_learned, K_learned, Sigma_learned, number_of_genes)            
    for gene_number=1:1:number_of_genes
        L_zero_col = L_learned;
        L_zero_col(:,gene_number)=0;
        L_learned(:, gene_number) = (X*Gamma_learned(:,gene_number) - Sigma_learned(gene_number)*(L_zero_col * K_learned(:,gene_number)))/(sum(Gamma_learned(:,gene_number))+Sigma_learned(gene_number)*K_learned(gene_number,gene_number));      
    end
end
            