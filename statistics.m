function [] = statistics(save_path,parameters_string, intervals,intervals_with_rp )
%Save all statistics (2 types hists and pseudoAUC to files
%
mkdir([save_path '/' parameters_string '/']);
thresholds=-10:0.5:10;
TPR = [];
FPR = [];
for  thres= thresholds
    TP=0;
    FP=0;
    FN=0;
    TN=0;
    hist=[];
    for i = 1:length(intervals)
        interval = intervals{i};
        classifer_triggered = false;
        for epoch = interval
            if(epoch.Q < thres)
                hist = [hist,epoch.dt_before_mov];
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
    histogram(hist,max(10,abs(max(hist) - min(hist))*10));
    title(sprintf('Threshold = %f\n',thres));
    saveas(gcf,[save_path '/' parameters_string '/' num2str(thres) '.png']);
    close;
    tmp_TPR.value=TP/(TP+FN);
    tmp_TPR.threshold = thres;
    tmp_FPR.value=FP/(FP+TN);
    tmp_FPR.threshold = thres;
    TPR = [TPR,tmp_TPR];
    FPR = [FPR,tmp_FPR];
end
plot(thresholds,[FPR.value],'b',thresholds,[TPR.value],'r'),legend('FPR','TPR');
title('TPR_FPR.png');
saveas(gcf,[save_path parameters_string 'TPR_FPR.png']);
close;
auc = trapz([FPR.value],[TPR.value]);
end

