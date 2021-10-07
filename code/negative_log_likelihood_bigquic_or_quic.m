%% negative_log_likelihood_bigquic.m
% this function calculates the objective function of BIGQUIC for arbitary
% inputs. 
% @ Inputs:
% L 
% K 
% lambda 
% @ Outputs:
% val 
function val = negative_log_likelihood_bigquic_or_quic(L, K, lambda)

    val=-1*log(det(K)) + trace(L'*L*K) + lambda * sum(sum(abs(K)));
end

