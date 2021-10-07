%% load Data
addpath(genpath('../../rosmapAD')) % dataset folder
prob_subject=load('methylationSNMnorm.mat');
probe_gene=load('methySNMtoGene100KB.mat');
expression_and_phenotype = load('expressionAndPhenotype');

% important: we modify the name to consistent (and well, more meaningfull) names before working with them
probe_gene_distance_matrix = probe_gene.xToGene;
probe_gene_distance_matrix(abs(probe_gene_distance_matrix)>99999)=0;
subject_probe_measurement_matrix = prob_subject.m.data';
probes = prob_subject.m.rowlabels;
genes = expression_and_phenotype.data.collabels;

%  rnd_perm_genes_index=randperm(length(genes));
%  genes=genes(rnd_perm_genes_index);
%  probe_gene_distance_matrix = probe_gene_distance_matrix(:,rnd_perm_genes_index);

chrs=prob_subject.m.chr;

