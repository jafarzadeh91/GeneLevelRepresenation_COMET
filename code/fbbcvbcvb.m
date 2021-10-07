L_correlation = {
[0.6661    0.7501    0.7734    0.7845    0.7906    0.7943    0.7966    0.7981    0.7991    0.7997    0.8000]
[ 0.6661    0.7128    0.7242    0.7270    0.7276    0.7272    0.7267    0.7261    0.7253    0.7245    0.7236]
[0.4970    0.5491    0.5777    0.6037    0.6272    0.6473    0.6641    0.6776    0.6890    0.6990    0.7079]
 [0.4970    0.5626    0.6005    0.6203    0.6320    0.6395    0.6451    0.6495    0.6530    0.6558    0.6580]
};
W_correlation = {
[ 0.8506    0.8833    0.8977    0.9052    0.9093    0.9117    0.9130    0.9135    0.9137    0.9134  ]
[ 0.8506    0.8810    0.8919    0.8956    0.8953    0.8931    0.8906    0.8880    0.8852    0.8823]
[0.6634    0.7054    0.7369    0.7590    0.7750    0.7885    0.8005    0.8109    0.8192    0.8255]
[0.6634    0.7149    0.7399    0.7536    0.7626    0.7688    0.7734    0.7769    0.7795    0.7814]
};

figure(1)
X = 1:1:11;
plot(X,L_correlation{1},X,L_correlation{2},X,L_correlation{3},X,L_correlation{4});
ylim([0 1])
title('Avg. Corr. between Estimated L and Ground Truth L')
xlabel('Coordinate Descent Iteration')
ylabel('Correlation')
grid on
legend({'wPCA init.', 'wPCA init. and $K=\O$','PCA init.','PCA init. with $K=\O$'},'Interpreter','latex')

figure(2)
X = 1:1:10;
plot(X,W_correlation{1},X,W_correlation{2},X,W_correlation{3},X,W_correlation{4});
ylim([0 1])
title('Avg. Corr. between Estimated W and Ground Truth W')
xlabel('Coordinate Descent Iteration')
ylabel('Correlation')
grid on
legend({'wPCA init.', 'wPCA init. and $K=\O$','PCA init.','PCA init. with $K=\O$'},'Interpreter','latex')


K_precision_wpca = [0.463070881342590;0.512526529819307;0.525910499353563;0.530551485541587;0.531451802743032;0.530375529891842;0.528360484307643;0.525979610017247;0.523471043169961;0.521168397038507];
K_recall_wpca = [0.727422398155205;0.724576335270676;0.707776810148764;0.698642243373783;0.695196249399037;0.693322093393925;0.691088004796872;0.691893056750258;0.691721197898699;0.691482934633112];
K_fscore_wpca = 2*(K_precision_wpca.*K_recall_wpca)./(K_precision_wpca+K_recall_wpca);


K_precision_pca = [0.347735150033456;0.372521062835037;0.391624645324604;0.409642558695802;0.423446753113011;0.432462361887890;0.437274678099976;0.440533645309824;0.443227173120830;0.445686111676398];
K_recall_pca = [0.321763850616205;0.345678227765666;0.354878392038868;0.371600914499505;0.398744013388066;0.425840145622244;0.451081901531284;0.475684177482468;0.499086526053735;0.523767604111011];
K_fscore_pca = 2*(K_precision_pca.*K_recall_pca)./(K_precision_pca+K_recall_pca);


figure(3)
X = 1:1:10;
plot(X,K_precision_wpca,X,K_precision_pca);
ylim([0 1])
title('Weighted Precision of K')
xlabel('Coordinate Descent Iteration')
ylabel('Precision')
grid on
legend('wPCA init.', 'PCA init.')


figure(4)
X = 1:1:10;
plot(X,K_recall_wpca,X,K_recall_pca);
ylim([0 1])
title('Weighted Recall of K')
xlabel('Coordinate Descent Iteration')
ylabel('Recall')
grid on
legend('wPCA init.', 'PCA init.')

