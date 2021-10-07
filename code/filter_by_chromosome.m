probes_on_this_chromosome=find(chrs==chr_number);
probe_gene_distance_matrix=probe_gene_distance_matrix(probes_on_this_chromosome,:);
subject_probe_measurement_matrix=subject_probe_measurement_matrix(:,probes_on_this_chromosome);
probes=probes(probes_on_this_chromosome);
chrs=chrs(probes_on_this_chromosome);

