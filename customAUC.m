function [ auc ] = customAUC( intervals,intervals_with_rp)
%customAUC calc auc due to one classifier out to interval
% Q=zeros(size(intervals));
% Label=zeros(size(intervals));
thresholds=-10:0.5:10;
% for interval=intervals
%     for epoch = deal(interval{:})
%         thresholds=[thresholds,epoch.Q1];
%     end
% end
% thresholds = unique(thresholds);
% thresholds = thresholds+eps;

% TPR=repmat(struct('value',1.0,'threshold',1.0),1,length(thresholds));
% FPR=repmat(struct('value',1.0,'threshold',1.0),1,length(thresholds));
TPR=[];
FPR=[];

for i = 1:length(thresholds)
    Q=[];
    Labels=[];
    for interval = intervals
        for epoch = deal(interval{:})
            if(~isnan(epoch.Q))&&(epoch.Q<thresholds(i))
                Q=[Q,epoch.Q];
                Labels=[Labels,epoch.rp];
                break;
            end
        end
    end
    
    
    tmp_TPR.value=length(Q(Labels == 1)<thresholds(i))/sum(Labels==1);
    tmp_TPR.threshold = thresholds(i);
    tmp_FPR.value=length(Q(Labels == 0)<thresholds(i))/sum(Labels==0);
    tmp_FPR.threshold = thresholds(i);
    TPR = [TPR,tmp_TPR];
    FPR = [FPR,tmp_FPR];
end
auc = trapz(FPR.value,TPR.value);
