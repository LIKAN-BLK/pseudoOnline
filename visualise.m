% 
% output1 = [];
% for i=1:size(classifierOutput1,2)
%     if (size(classifierOutput1{i},2)==322)
%         output1 = [output1;classifierOutput1{i}];
%     end
% end
% 
% output2 = [];
% for i=1:size(classifierOutput2,2)
%     if (size(classifierOutput2{i},2)==322)
%         output2 = [output2;classifierOutput2{i}];
%     end
% end
% time =(1:size(output1,2))/200 - size(output1,2)/200;
% plot(time,mean(output1,1))
% plot(time,mean(output1,1))
% 
% cumm1=zeros(size(output1));
% for i=1:5:size(output1,2)-5
%     cumm1(:,i:i+5-1)=repmat(sum(output1(:,i:i+5-1),2),1,5);
% end
% surf(time,1:size(cumm1,1),cumm1);
% 
% 
% surf(time,1:size(output1,1),output1);

%  tmp_start = find(~isnan(clfOut1));
%    if(~isempty(tmp_start))
%     tmp_end = tmp_start(end);
%     tmp_start = tmp_start(1);
%     if(all(~isnan(clfOut1(tmp_start:tmp_end))))
%         classifierOutput1{i} = clfOut1(tmp_start:tmp_end);
%     end
%    end
%    
%    if(~isempty(tmp_start))
%     tmp_start = find(~isnan(clfOut2));
%     tmp_end = tmp_start(end);
%     tmp_start = tmp_start(1);
%     if(all(~isnan(clfOut2(tmp_start:tmp_end))))
%         classifierOutput2{i} = clfOut2(tmp_start:tmp_end);
%     end
%    end
function [] = visualise(intervals,intervals_with_rp,type)
    %Histogram for times when classifier triggers
    if(strcmp(type,'histogram'))
        hists = [];
        thresholds = -10:0.5:10;
        for  thres= thresholds
            hist=[];
            for i = 1:length(intervals)
                interval = intervals{i};
                classifer_triggered = false;
                for epoch = interval
                   if(epoch.Q<thres)
                        hist = [hist,epoch.dt_before_mov];
                        classifer_triggered = true;
                        break;
                   end
                end
                if(~classifer_triggered) && intervals_with_rp(i)
                    hist = [hist,15];                           %if classifier not triggered, we consider, that it triggered VERY late (10s)
                end
            end
            hists_elem.thres = thres;
            hists_elem.hist = hist;
            histogram(hist);
            title(sprintf('Threshold = %f',thres));
            hists = [hists,hist];
        end
    end
    
    if(strcmp(type,'hist_clf_output'))
        target_hist=[];
        nontarget_hist=[];
        for i = 1:length(intervals)
            interval = intervals{i};
            for epoch = interval
               if(~isnan(epoch.Q))
                   switch(epoch.rp)
                       case 1
                           target_hist = [target_hist,epoch.Q];
                       case 0
                           nontarget_hist=[nontarget_hist,epoch.Q];
                   end
               end
            end            
        end
        histogram(target_hist),hold on,histogram(nontarget_hist);
        legend('Target','NonTarget')
%         histogram(target_hist);
%         title(sprintf('Target'));
%         figure();
%         histogram(nontarget_hist);
%         title(sprintf('NonTarget'));
    end   
end




