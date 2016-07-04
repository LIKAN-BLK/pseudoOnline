function [prin_comp,classifier,opt_thr ] = learn( X1, X0)
%get_classifier_mat Returns classifier matrix
    
    N0 = size(X0, 1);
    N1 = size(X1, 1);
    
    num_of_pca = 150;
    [prin_comp,X] = princomp([X0;X1],num_of_pca);
    X = X(:,1:num_of_pca);
    prin_comp = prin_comp(:,1:num_of_pca);
    obj = train_shrinkage(X,[ones(N0,1);2*ones(N1,1)]);
    classifier = obj.W;
    Q0 = (X0*prin_comp)*classifier;
    Q1 = (X1*prin_comp) * classifier;
%     [Xroc,Yroc,T,square,opt_roc_point] = perfcurve([ones(N1,1); zeros(N0,1)], [Q1; Q0], 0);
    
%     ind = ([Xroc,Yroc] == repmat(opt_roc_point,size(Xroc,1),1));
%     opt_thr = T(find(ind(:,1) & ind(:,2)));
    opt_thr=[];

end

