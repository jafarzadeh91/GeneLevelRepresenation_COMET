
%L(:,2)=L(:,1)+randn()*0.2;
%L(:,3)=L(:,1)*randn();
%K=[5 3 2;3 6 3;2 3 4];
K=[1 0 0;0 1 0; 0 0 1];
L = mvnrnd(zeros(20, 3),K );

    
W=10000*randn(5,3);
X=L*W';
X_zscored=zscore(X);
[coeff,score,latent,tsquared,explained,mu] =pca(X_zscored');
L_learned=randn(20,3);
W_learned=randn(5,3);
for i=1:1:10
    disp(sum(sum(abs(X-L_learned*W').^2)))
    
 W_learned= (L_learned'*L_learned)^(-1)*L_learned'*X;
 W_learned=W_learned';
    disp(sum(sum(abs(X-L_learned*W').^2)))
 L_learned = (X*W_learned)/(W_learned'*W_learned);
disp(sum(sum(abs(X-L_learned*W_learned').^2)))
end

corr_matrix_ours=zeros(3,3);
corr_matrix_ours_w=zeros(3,3);
corr_matrix_pca=zeros(3,3);
for i=1:1:3
    for j=1:1:3
        corr_matrix_ours(i,j)=corr(L_learned(:,i),L(:,j));
                corr_matrix_ours_w(i,j)=corr(W_learned(:,i),W(:,j));

        corr_matrix_pca(i,j)=corr(coeff(:,i),L(:,j));

    end
end


corrs_ours=(corr_matrix_ours(sub2ind(size(corr_matrix_ours),[1 2 3],munkres(1-abs(corr_matrix_ours)))))
corrs_ours_w=(corr_matrix_ours_w(sub2ind(size(corr_matrix_ours),[1 2 3],munkres(1-abs(corr_matrix_ours)))))

corrs_pca=(corr_matrix_pca(sub2ind(size(corr_matrix_pca),[1 2 3],munkres(1-abs(corr_matrix_pca)))))
