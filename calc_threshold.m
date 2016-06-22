function [hist_F1_threshold] = calc_threshold(intervals,intervals_with_rp)
%CALC_THRESHOLD Summary of this function goes here
%   Detailed explanation goes here
thresholds=-10:0.5:10;
max_hist_F1 = 0;
for  thres= thresholds
    tmp_hist=[];
    for i = 1:length(intervals)
        interval = intervals{i};
        classifer_triggered = false;
        for epoch = interval
            if(epoch.Q < thres)
                tmp_hist = [tmp_hist,epoch.dt_before_mov];
                classifer_triggered = true;
                break;
            end
        end
        
        if(~classifer_triggered)
            if(intervals_with_rp(i) == 1)
                tmp_hist = [tmp_hist,5]; % FN considered as VERY late classifier triggering
            end
        end    
    end
    if ~isempty(tmp_hist)
        curr_hist_F1 = length(tmp_hist(tmp_hist > -2 & tmp_hist <-0.5))/length(tmp_hist);
        if(curr_hist_F1 > max_hist_F1)
            hist_F1_threshold = thres;
            max_hist_F1 = curr_hist_F1;
        end
    end
end
end

