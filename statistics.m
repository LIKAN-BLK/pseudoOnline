function [ACC_threshold,F1_threshold,hist_F1_threshold] = statistics(save_path,parameters_string, intervals,intervals_with_rp )
%Save all statistics (2 types hists and pseudoAUC to files
%
mkdir([save_path '/' parameters_string '/']);
thresholds=-10:0.5:2;
TPR = [];
FPR = [];
F1 = [];
ACC = [];
hist_F1=[];

TP=[];
FP=[];
FN=[];
TN=[];
for  thres= thresholds
    curr_TP=0;
    curr_FP=0;
    curr_FN=0;
    curr_TN=0;
    tmp_hist=[];
    for i = 1:length(intervals)
        interval = intervals{i};
        classifer_triggered = false;
        for epoch = interval
            if(epoch.Q < thres)
                tmp_hist = [tmp_hist,epoch.dt_before_mov];
                switch (epoch.rp)
                    case 1
                        curr_TP = curr_TP + 1;
                    case 0
                        curr_FP = curr_FP + 1;
                    case -1
                        curr_FN = curr_FN + 1;
                end
                classifer_triggered = true;
                break;
            end
        end
        
        if(~classifer_triggered)
            if(intervals_with_rp(i) == 0)
                curr_TN = curr_TN + 1;
            else
                tmp_hist = [tmp_hist,5];
                curr_FN = curr_FN + 1;
            end
        end    
    end
    
    TP=[TP,curr_TP];
    FP=[FP,curr_FP];
    FN=[FN,curr_FN];
    TN=[TN,curr_TN];
    
    
    curr_F1=1*curr_TP/(1*curr_TP + curr_FN +curr_FP);
    curr_ACC = (curr_TP+curr_TN)/(curr_TP+curr_TN+curr_FP+curr_FN);
    F1 = [F1,curr_F1];
    ACC = [ACC,curr_ACC];
    
    curr_hist_F1 = length(tmp_hist(tmp_hist >= -2 & tmp_hist <=-0.5))/length(tmp_hist);
    hist_F1=[hist_F1,curr_hist_F1];
    
    histogram(tmp_hist,max(5,round(abs(max(tmp_hist) - min(tmp_hist))*10)));
    title(sprintf('Threshold = %f,F1=%f,ACC=%f,pseudoF1=%f\n',thres,curr_F1,curr_ACC,curr_hist_F1));
    saveas(gcf,[save_path '/' parameters_string '/' num2str(thres) '.png']);
    close;
    

end


[~,ACC_index]=max(ACC);
ACC_threshold = thresholds(ACC_index);
[~,F1_index]=max(F1);
F1_threshold = thresholds(F1_index);
[~,hist_F1_index]=max(hist_F1);
hist_F1_threshold = thresholds(hist_F1_index);
end

