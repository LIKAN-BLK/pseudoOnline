function [ auc ] = customAUC( intervals,thresholds)
%customAUC calc auc due to one classifier out to interval
for thr = thresholds
    Q=[];
    Labels=[];
    for interval = intervals
        for epoch = epochs
            if(epoch.Q1<thr)
                Q=[Q,epoch.Q1]
                Labels=[Labels,epoch.rp]
            end
        end
    end
    TPR(thr)=calcTPR(Q,Labels);
    FPR(thr)=calcTPR(Q,Labels);
    


end

