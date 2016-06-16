function [ auc ] = customAUC( intervals,intervals_with_rp)
%customAUC calc auc due to one classifier out to interval
thresholds=-10:0.5:10;
TPR = [];
FPR = [];
max_F1=0;
max_ACC=0;
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
    curr_F1 = 2*TP/(2*TP + FN +FP);
    curr_ACC = (TP+TN)/(TP+TN+FP+FN);
    if(curr_F1 > max_F1)
        F1_threshold = thres;
        max_F1 = curr_F1;
    end
    if(curr_ACC > max_ACC)
        ACC_threshold = thres;
        max_ACC = curr_ACC;
    end
    tmp_TPR.value=TP/(TP+FN);
    tmp_TPR.threshold = thres;
    tmp_FPR.value=FP/(FP+TN);
    tmp_FPR.threshold = thres;
    TPR = [TPR,tmp_TPR];
    FPR = [FPR,tmp_FPR];
end
plot(thresholds,[FPR.value],'b',thresholds,[TPR.value],'r'),legend('FPR','TPR');
auc = trapz([FPR.value],[TPR.value]);



% function [ auc ] = customAUC( intervals,intervals_with_rp)
% %customAUC calc auc due to one classifier out to interval
% thresholds=-10:0.5:10;
% TPR = [];
% FPR = [];
% for  thres= thresholds
%     TP=0;
%     FP=0;
%     FN=0;
%     TN=0;
%     for i = 1:length(intervals)
%         interval = intervals{i};
%         classifer_triggered = false;
%         for j = 1:3:size(interval,2)-(3-1)
%             epochs_window = interval(j:j+3-1);
%             if (sum([epochs_window.Q] < thres) > 2)
%                 rp_rate = mode([epochs_window.rp]); %To convert from logical to double
%                 switch (rp_rate)
%                     case 1
%                         TP = TP + 1;
%                     case 0
%                         FP = FP + 1;
%                         
%                     case -1
%                         FN = FN + 1;
%                 end
%                 classifer_triggered = true;
%                 break;
%             end
%         end
%         
%         if(~classifer_triggered)
%             if(intervals_with_rp(i) == 0)
%                 TN = TN + 1;
%             else
%                 FN = FN + 1;
%             end
%         end    
%     end
%     
%     
%     tmp_TPR.value=TP/(TP+FN);
%     tmp_TPR.threshold = thres;
%     tmp_FPR.value=FP/(FP+TN);
%     tmp_FPR.threshold = thres;
%     TPR = [TPR,tmp_TPR];
%     FPR = [FPR,tmp_FPR];
% end
% plot(thresholds,[FPR.value],'b',thresholds,[TPR.value],'r'),legend('FPR','TPR');
% auc = trapz([FPR.value],[TPR.value]);
