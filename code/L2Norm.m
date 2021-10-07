function [c,ceq,c_gred,ceq_gred] = L2Norm(w_vec)
%L1NORM Summary of this function goes here
%   Detailed explanation goes here
ceq=sum(abs(w_vec).^2)-1;

c=[];
c_gred=[];
ceq_gred = 2*sum(w_vec);
end

