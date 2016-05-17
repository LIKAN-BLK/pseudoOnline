function [W] = train(X1,X0)
%
% X0: [Nsamples0 * Nfeats]
% X1: [Nsamples1 * Nfeats]
%

N0 = size(X0, 1);
N1 = size(X1, 1);
X = [X0;X1]; 
Y = [ones(1,N0) 2*ones(1,N1)]';
for i = 1:CV.NumTestSets
    trIdx = CV.training(i);
    tstIdx = CV.test(i);
    Ytr = Y(trIdx, :);
    Ytst = Y(tstIdx, :);
            
    Xtr = X(trIdx, :);
    Xtst = X(tstIdx, :);  
    %chosing R2 relevant feats
%     ind_of_correl_feats = find_meaningful_feats(Xtr,Ytr,50);
%     Xtr = Xtr(:,ind_of_correl_feats);
%     Xtst  = Xtst(:,ind_of_correl_feats);
    
    N0tr = sum(Ytr == 1);
    N1tr = sum(Ytr == 2);
    N0tst = sum(Ytst == 1);
    N1tst = sum(Ytst == 2);

    % train
    obj = train_shrinkage(Xtr, Ytr);
    W(:,:,i) = obj.W;
    %W(:,:,i) = rand(size(obj.W))*0.02-0.01;
    
    % calc acc on train sample
    Q = Xtr*W(:,:,i);    
    Q0 = Q(Ytr == 1);
    Q1 = Q(Ytr == 2);
%     
%     ths = Q + eps;
%     ths = sort(ths);
%     for k = 1:length(ths)                
%         sens_tr(k) = length(find(Q1 <= ths(k))) / N1tr;
%         spec_tr(k) = length(find(Q0 > ths(k))) / N0tr;    
%     end;
%     idx = find(spec_tr >= 0.99, 1, 'last');
%     th_opt(i) = ths(idx); 
%     sens_tr(i) = sens_tr(idx);
%     spec_tr(i) = spec_tr(idx);    
    [aucXtr,aucYtr, ~, auc_tr(i)] = perfcurve([ones(N1tr,1); zeros(N0tr,1)], [Q1; Q0], 0);
    [aucXtr, index]=unique(aucXtr);
    meanAucYtr = meanAucYtr + interp1(aucXtr,aucYtr(index),meanAucX);
    
    % test
    Q = Xtst*W(:,:,i);
    Q0 = Q(Ytst == 1);
    Q1 = Q(Ytst == 2);
%     sens_tst(i) = length(find(Q1 <= th_opt(i))) / N1tst;
%     spec_tst(i) = length(find(Q0 > th_opt(i))) / N0tst;    
%     acc_tst(i) = (sens_tst(i) * N1tst + spec_tst(i) * N0tst) / (N1tst + N0tst);
    [aucXtst,aucYtst, ~, auc_tst(i)] = perfcurve([ones(N1tst,1); zeros(N0tst,1)], [Q1; Q0], 0);
    [aucXtst, index]=unique(aucXtst);
    meanAucYtst = meanAucYtst + interp1(aucXtst,aucYtst(index),meanAucX);
end

% spec.tr = [mean(spec_tr) std(spec_tr)];
% sens.tr = [mean(sens_tr) std(sens_tr)];
spec.tr = 0;
sens.tr = 0;
acc.tr = 0;
% spec.tst = [mean(spec_tst) std(spec_tst)];
% sens.tst = [mean(sens_tst) std(sens_tst)];
% acc.tst = [mean(acc_tst),std(acc_tst)];
spec.tst = 0;
sens.tst = 0;
acc.tst = 0;
auc.x = meanAucX;
auc.tr.square = [mean(auc_tr) std(auc_tr)];
auc.tr.y = meanAucYtr/nfolds;
auc.tst.square = [mean(auc_tst) std(auc_tst)];
auc.tst.y = meanAucYtst/nfolds;
auc.all = 0;


params.W = mean(W, 3);
%params.th = mean(th_opt);

% fid = fopen('../res/acc3.txt', 'w');
% fprintf(fid, 'Train:\n');
% fprintf(fid, ' Sensitivity: %f +- %f\n', sens.tr);
% fprintf(fid, ' Specificity: %f +- %f\n', spec.tr);
% fprintf(fid, '\n');
% fprintf(fid, 'Test:\n');
% fprintf(fid, ' Sensitivity: %f +- %f\n', sens.tst);
% fprintf(fid, ' Specificity: %f +- %f\n', spec.tst);
% fclose(fid);

