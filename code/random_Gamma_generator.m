%% random_Gamma_generator.m
% In this file, the method receive the Z assignment matrix and construct an
% artificial random Gamma function by sampling from a flat dirichlet
% distribuion.
% @ Inputs
% Z : the matrix showing the assinment of each probe to the genes which
% have the prob in their's neighborhood
% @ Ouputs
% random_Gamma : 
function random_Gamma = random_Gamma_generator(Gamma,gamma_variance, trace)
%RANDOM_GAMMA_GENERATOR Summary of this function goes here
%   Detailed explanation goes here
    random_Gamma = Gamma;
    
    for i=1:1:length(random_Gamma)
        elements_ind = find(Gamma(i, :));
        number_of_elements = length(elements_ind);
        while(true)
            rnd_vec = randn(1,number_of_elements)*gamma_variance;
            new_vec = Gamma(i,elements_ind) + rnd_vec;
            non_zeros_inds=find(new_vec~=0);
            new_vec(new_vec<=0) = 0.01;
            new_vec(non_zeros_inds) = new_vec(non_zeros_inds) - (sum(new_vec)-1)/length(non_zeros_inds);
            
            if(~isempty(find(new_vec<=0, 1))) % repeat
                if(trace == 1)
%                   disp(strcat('an attemp failed to generate high quality random Gamma... Please wait.,,')) 
                end
                continue;
            else
                if(trace == 1)
                   disp(strcat('random Gamma simulation: ', num2str(i),'/',num2str(length(random_Gamma)))) 
                end
                random_Gamma(i,elements_ind)=new_vec;
                break;
            end
        end
    end
end

