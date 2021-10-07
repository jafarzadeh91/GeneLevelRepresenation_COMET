function [c,ceq] = L1Norm(w_vec)
%L1NORM Summary of this function goes here
%   Detailed explanation goes here
ceq=1;
if(sum(abs(w_vec))==1)
   ceq=0; 
end
c=[];
end

