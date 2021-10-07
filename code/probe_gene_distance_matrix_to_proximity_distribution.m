%% this function converts the prb_gene_dis_mtrx(d*q) to binary matrix Z (d*q) and its associated normalized Gamma (d*q) matrix.
% @ Inputs
% prb_gene_dis_mtrx: (d*q) matrix shows the distance of each gene from all probes in
% its neighborhood (e.g. 1MB window aroound the gene).
% chunk_size: chunck size is an arbitrary value 1<chunk_size<d to help the
% memory performance of method. Set it to the biggest number that the memory of your machine can handle :)
% @ Outputs
% Z: (d*q) binary matrix with just one 1 value in each row showing the
% nearest gene to each prob
% Gamma: a normalized version of Z with sum of one in each row that show
% the distance-based association of each prob to a set of genes. for genes
% that the prob is not in its neighborhood, the value is 0.
function [Gamma] = probe_gene_distance_matrix_to_proximity_distribution(prb_gene_dis_mtrx, chunck_size)
    number_of_probes = size(prb_gene_dis_mtrx,1);
    
    Gamma = sparse(size(prb_gene_dis_mtrx,1),size(prb_gene_dis_mtrx,2));
    number_of_chunks = ceil(number_of_probes/chunck_size);

    for chunk_index=1:1:number_of_chunks
       
       interval = (chunk_index-1)*chunck_size+1:min(chunk_index*chunck_size, number_of_probes);
       rows = full(prb_gene_dis_mtrx(interval,:));
       rows_abs = abs(rows);
       rows_abs(rows_abs~=0)=1./rows_abs(rows_abs~=0);

       sum_vec_rows_abs = sum(rows_abs,2);
       Gamma(interval,:) = bsxfun(@rdivide,rows_abs,sum_vec_rows_abs(:));           
    
    end
end

