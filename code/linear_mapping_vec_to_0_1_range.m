function out= linear_mapping_vec_to_0_1_range(inp)
    min_val = min(inp);
    out = inp - min_val;
    out = out/max(out);            
end

