function [ACC_threshold,F1_threshold,hist_F1_threshold] = statistics(save_path,parameters_string, intervals,intervals_with_rp )
%Save all statistics (2 types hists and pseudoAUC to files
%
mkdir([save_path '/' parameters_string '/']);
thresholds=-10:0.5:10;
TPR = [];
FPR = [];
max_F1 = 0;
max_ACC = 0;
max_hist_F1 = 0;
for  thres= thresholds
    TP=0;
    FP=0;
    FN=0;
    TN=0;
    tmp_hist=[];
    for i = 1:length(intervals)
        interval = intervals{i};
        classifer_triggered = false;
        for epoch = interval
            if(epoch.Q < thres)
                tmp_hist = [tmp_hist,epoch.dt_before_mov];
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
                tmp_hist = [tmp_hist,5];
                FN = FN + 1;
            end
        end    
    end
    if ~isempty(tmp_hist)
        histogram(tmp_hist,max(5,round(abs(max(tmp_hist) - min(tmp_hist))*10)));
        title(sprintf('Threshold = %f\n',thres));
        saveas(gcf,[save_path '/' parameters_string '/' num2str(thres) '.png']);
        close;
        curr_hist_F1 = length(tmp_hist(tmp_hist > -2 & tmp_hist <-0.5))/length(tmp_hist);
        if(curr_hist_F1 > max_hist_F1)
            hist_F1_threshold = thres;
            max_hist_F1 = curr_hist_F1;
        end
    end
    
    
    curr_F1 = TP/(TP + FN +FP);
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
title('TPR_FPR');
saveas(gcf,[save_path 'TPR_FPRs/' parameters_string 'TPR_FPR.png']);
close;
auc = trapz([FPR.value],[TPR.value]);
end

