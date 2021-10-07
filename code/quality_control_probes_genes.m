%% We have some probes with zero gene OR genes with zero probes in their neighborhood. We will remove them first
genes_with_probe_index = find(any(probe_gene_distance_matrix,1));
probes_with_gene_index = find(any(probe_gene_distance_matrix(:,genes_with_probe_index),2));
probe_gene_distance_matrix = probe_gene_distance_matrix(probes_with_gene_index, genes_with_probe_index);
subject_probe_measurement_matrix = subject_probe_measurement_matrix(:, probes_with_gene_index);
probes = probes(probes_with_gene_index,:);
chrs = chrs(probes_with_gene_index);
genes = genes(genes_with_probe_index);

%% define X matrix, number of genes, probes and subjects
number_of_genes = size(probe_gene_distance_matrix,2); % q=13484 genes
number_of_probes = size(probe_gene_distance_matrix,1); % d=420103 probes
number_of_subjects = size(subject_probe_measurement_matrix,1); % n=702 subjects

