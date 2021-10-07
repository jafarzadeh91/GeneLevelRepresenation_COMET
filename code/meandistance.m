function res = meandistance(vec1,vec2,power)
diff=abs(vec1-vec2);
diff=diff.^power;
res = sum(diff)/length(vec1);   
end

