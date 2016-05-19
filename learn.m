function [ W,mu,sigma,opt_thr ] = learn( eegT, eegNT,w_size)
%get_classifier_mat Returns classifier matrix

    tmp = reshape(cat(3,eegNT,eegT),size(eegNT,1),size(eegNT,2)*(size(eegNT,3)+size(eegT,3)));
    mu = mean(tmp,2);
    sigma = std(tmp,0,2);
    clear tmp;
    
    eegT = permute(eegT,[2,1,3]);  %now data dims is time x channel x trial
    eegNT = permute(eegNT,[2,1,3]);
    
    
    
    fs = 500;
    [X1] = get_feats(eegT, w_size);
    [X0] = get_feats(eegNT, w_size);
    
    N0 = size(X0, 1);
    N1 = size(X1, 1);
    obj = train_shrinkage([X0;X1],[ones(N0,1);2*ones(N1,1)]);
    W = obj.W;
    Q0 = X0*W;
    Q1 = X1 * W;
    [Xroc,Yroc,T,square,opt_roc_point] = perfcurve([ones(N1,1); zeros(N0,1)], [Q1; Q0], 0);
    
    ind = ([Xroc,Yroc] == repmat(opt_roc_point,size(Xroc,1),1));
    opt_thr = T(find(ind(:,1) & ind(:,2)));

end

