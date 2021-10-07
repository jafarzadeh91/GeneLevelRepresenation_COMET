clear;
libraries_header;
chromosomes_number=[22];

%% dataset parameters
data_instance_identifier='default_nips';
sigma_scaling_factor = 0.00000000000000001;
%% load real data to simulate data based on that!
data_name='rosmap'; %% or 'KoborDNAm'
method_name = 'nips'; % or 'slggm'

for i=1:1:length(chromosomes_number)
        chr_number = chromosomes_number(i);

        %% data save
        root_addr=strcat('data/synth_data/',data_name);
        file_addr  = strcat(root_addr,'/','chr_',num2str(chr_number),'_',data_instance_identifier);
        load(strcat(file_addr,'.mat'));
        
        mkdir(file_addr);
        imagesc(abs(K_current_dataset));
        xlabel('Genes #')
        ylabel('Genes #')
        print(strcat(file_addr,'/K_current_dataset'),'-dpdf','-bestfit')
        imagesc(abs(K_current_dataset)>0)
        xlabel('Genes #')
        ylabel('Genes #')
        print(strcat(file_addr,'/K_current_dataset_binary'),'-dpdf','-bestfit')
end
