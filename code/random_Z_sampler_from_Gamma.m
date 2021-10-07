function Z = random_Z_sampler_from_Gamma(Gamma)
    Z = sparse(size(Gamma,1), size(Gamma,2));
    for row_ind = 1:1:size(Gamma,1)
        row = Gamma(row_ind,:);
        row_non_zeros=find(row~=0);
        row_non_zeros_values=full(row(row_non_zeros));
        gene_index = discretesample(row_non_zeros_values,1);
        Z(row_ind, row_non_zeros(gene_index)) = 1;
    end
end

