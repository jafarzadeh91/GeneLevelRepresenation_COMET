function K = generate_sparse_psd_block_matrix(number_of_elements,value_cut_off,min_eigen_value, number_of_blocks)
    K = sparse(number_of_elements, number_of_elements);
    blocks_indices=[1 randperm(number_of_elements,number_of_blocks) number_of_elements+1];
    blocks_indices = sort(blocks_indices);
    
    for i=1:1:length(blocks_indices)-1
        start_ind=blocks_indices(i);
        end_ind=blocks_indices(i+1);
        block_matrix=sprandsym(end_ind-start_ind, end_ind-start_ind);
        block_matrix(abs(block_matrix)<value_cut_off)=0;
        diagonal_indices= 1==eye(size(block_matrix));
        block_matrix(diagonal_indices)=rand(length(diagonal_indices),1);
        expected_min_eig_value=min_eigen_value;
        real_min_eig_value=min(eig(block_matrix));    
        block_matrix(diagonal_indices)=block_matrix(diagonal_indices)+(expected_min_eig_value-real_min_eig_value);
        K(start_ind:end_ind-1, start_ind:end_ind-1)=block_matrix;
    end
end

