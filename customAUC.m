% function [ auc ] = customAUC( intervals,intervals_with_rp)
% %customAUC calc auc due to one classifier out to interval
% % Q=zeros(size(intervals));
% % Label=zeros(size(intervals));
% thresholds=-10:0.5:10;
% % for interval=intervals
% %     for epoch = deal(interval{:})
% %         thresholds=[thresholds,epoch.Q1];
% %     end
% % end
% % thresholds = unique(thresholds);
% % thresholds = thresholds+eps;
% 
% % TPR=repmat(struct('value',1.0,'threshold',1.0),1,length(thresholds));
% % FPR=repmat(struct('value',1.0,'threshold',1.0),1,length(thresholds));
% TPR=[];
% FPR=[];
% 
% for i = 1:length(thresholds)
%     Q=[];
%     Labels=[];
%     for interval = intervals
%         for epoch = deal(interval{:})
%             if(~isnan(epoch.Q))&&(epoch.Q<thresholds(i))
%                 Q=[Q,epoch.Q];
%                 Labels=[Labels,epoch.rp];
%                 break;
%             end
%         end
%     end
%     
%     
%     tmp_TPR.value=length(Q(Labels == 1)<thresholds(i))/sum(Labels==1);
%     tmp_TPR.threshold = thresholds(i);
%     tmp_FPR.value=length(Q(Labels == 0)<thresholds(i))/sum(Labels==0);
%     tmp_FPR.threshold = thresholds(i);
%     TPR = [TPR,tmp_TPR];
%     FPR = [FPR,tmp_FPR];
% end
% auc = trapz(FPR.value,TPR.value);
function [ auc ] = customAUC( intervals,intervals_with_rp)
%customAUC calc auc due to one classifier out to interval
thresholds=-10:0.5:10;
TPR = [];
FPR = [];
for  thres= thresholds
    TP=0;
    FP=0;
    FN=0;
    TN=0;
    for i = 1:length(intervals)
        interval = intervals{i};
        classifer_triggered = false;
        for epoch = interval
            if(epoch.Q < thres)
                switch (epoch.rp)
                    case 1
                        TP = TP + 1;
                    case 0
                        FP = FP + 1;
                    case -1
                        FN = FN + 1;
                end
                classifer_triggered = true;
                break;
            end
        end
        
        if(~classifer_triggered)
            if(intervals_with_rp(i) == 0)
                TN = TN + 1;
            else
                FN = FN + 1;
            end
        end    
    end
    
    
    tmp_TPR.value=TP/(TP+FN);
    tmp_TPR.threshold = thres;
    tmp_FPR.value=FP/(FP+TN);
    tmp_FPR.threshold = thres;
    TPR = [TPR,tmp_TPR];
    FPR = [FPR,tmp_FPR];
end
plot(thresholds,[FPR.value],'b');
hold on, plot(thresholds,[TPR.value],'r'),hold off;
auc = trapz([FPR.value],[TPR.value]);
