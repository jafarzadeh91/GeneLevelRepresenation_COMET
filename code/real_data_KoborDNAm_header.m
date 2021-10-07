%% load Data
data_name='KoborDNAm';

prob_subject=load('MSK_probes_subjects_new_age_study_regressed_out.mat');
prob_subject = prob_subject.prob_subject_age_study_regressed_out;

%prob_subject=load('MSK_probes_subjects_new.mat');

prob_gene=load('probe_gene_distance_matrix_KoborDNAm_100_kb.mat');
subject_phenotype = load('MSK_subjects_phenotypes_new');

%load('../../UCSC/HumanMethylation45015017482v12hg19.mat');

%% important: we modify the name to consistent (and well, more meaningfull) names before working with them
probe_gene_distance_matrix = prob_gene.probe_gene_distance_matrix;
subject_probe_measurement_matrix = prob_subject.combat_probes';

%% gender specific
%gender = 'F'; % or 'F'
% if(lower(gender)=='f' || lower(gender)=='m')
%     data_name = strcat(data_name,lower(gender));
%     subjects_with_specified_gender = find(subject_phenotype.MSKmeta.Sex == upper(gender));
%     subject_probe_measurement_matrix = subject_probe_measurement_matrix(subjects_with_specified_gender, :);
%     
% end
chrs = grp2idx( prob_gene.rows_metadata.CHR)-3;
%chrs = ones(length(prob_gene.rows_metadata.CHR),1)*chr_number;

if(chr_number ~=0)
   prbs_on_chr=find(prob_gene.rows_metadata.CHR == categorical(chr_number));
   probe_gene_distance_matrix = probe_gene_distance_matrix(prbs_on_chr,:);
   subject_probe_measurement_matrix = subject_probe_measurement_matrix(:,prbs_on_chr);
   probes = prob_subject.rows(prbs_on_chr);
   chrs = chrs(prbs_on_chr);
   %data_name = strcat(data_name,'_chr_',num2str(chr_number));
else 
    probes = prob_subject.rows;
end


genes = prob_gene.cols;

% 
% if(signed_effect==1)
%     chrs=repelem(chrs,2);
%     genes=repelem(genes,2);
%     K_prior=repelem(K_prior,2,2);
% 
%     probe_gene_distance_matrix=repelem(probe_gene_distance_matrix,1,2);
%     positive_elements=1:2:size(probe_gene_distance_matrix,2);
%     negative_elements=2:2:size(probe_gene_distance_matrix,2);
%     
%     temp=probe_gene_distance_matrix(:,positive_elements);
%     temp(temp<0)=0;
%     probe_gene_distance_matrix(:,positive_elements)=temp;
%     
%     temp=probe_gene_distance_matrix(:,negative_elements);
%     temp(temp>0)=0;
%     probe_gene_distance_matrix(:,negative_elements)=temp;
%     
%     clear temp positive_elements negative_elements
% end


clear prob_subject probe_gene;