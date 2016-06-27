function [] = pseudoOnlineExp3(data_path,w_size_time,target_start,target_end,ch_to_use,rel_thres)


%%LEARNING CLASSIFIER
fs = 200; %Hz
baseline_time=200; %200ms
[prin_comp,classifier] = learn_classifier(data_path,fs,baseline_time,w_size_time,target_start,target_end,ch_to_use,rel_thres);


end

function [prin_comp,classifier] = learn_classifier(data_path,fs,baseline_time,w_size_time,target_start,target_end,ch_to_use,rel_thres)

[eegTRP,eegNT] = cut_epochs_4learn(data_path,fs,w_size_time,...
    baseline_time,target_start,target_end,ch_to_use,rel_thres);
[prin_comp,classifier, opt_thr] = learn(eegTRP,eegNT,w_size_time,fs);
end


function [] = apply_classifier(data_path,w_size_time,rp_start,rp_end,target_start,target_end,ch_to_use,rel_thres)
%%APPLYING CLASSIFIER

parameters_string = sprintf('s_%de_%dw_%d',target_start,target_end,w_size_time);
save_path = [data_path, 'results/'];
mkdir(save_path);
mkdir([save_path 'clf_out_hist_pics/']);
mkdir([save_path 'clf_out_hist_data/']);
mkdir([save_path 'TPR_FPRs/']);


[data,EOG,mask] = loaddata(data_path,ch_to_use);

[~,grad] = gradient(data);


rp_len = 2000 * fs /1000;

starts = find(mask == 12);
ends = find(mask == 13);

w_step = 10 * fs/1000; %10ms

intervals_rp_mask = nan(1,length(ends));
hist_clf_output_t=[];
hist_clf_output_nt=[];

for i=1:size(ends,2)
   interval = data(starts(i):ends(i),:);
   movement = find(mask(starts(i):ends(i)) == 10);
   if ((movement - rp_len +1) < 0)
       continue
   end
   
   rp_mask = zeros(size(mask(starts(i):ends(i))));
   if ~isempty(movement)
    rp_start = max([1,movement + rp_start*fs + 1]);
    rp_indexes = [rp_start:movement+rp_end*fs]; %Our interval of interests ends 0.5s before movement.
    rp_mask(rp_indexes) = 1;
    rp_mask(movement)=10;
   end
   grad_interval = grad(starts(i):ends(i),:);
   EOG_interval = EOG(starts(i):ends(i));
   
       
       
   [intervals{i},intervals_rp_mask(i),tmp_hist_clf_output_t,tmp_hist_clf_output_nt] = ...
       process_interval(interval,rp_mask,grad_interval,EOG_interval,w_size_time,baseline_time,fs,w_step,prin_comp,classifier,rel_thres);
    
    hist_clf_output_t = [hist_clf_output_t,tmp_hist_clf_output_t];
    hist_clf_output_nt = [hist_clf_output_nt,tmp_hist_clf_output_nt];
       
end
[~,~,~,auc] = perfcurve([zeros(size(hist_clf_output_t)),ones(size(hist_clf_output_nt))],[hist_clf_output_t,hist_clf_output_nt],1);

pseudoFisher = (mean(hist_clf_output_t)-mean(hist_clf_output_nt))^2/(cov(hist_clf_output_t)+cov(hist_clf_output_nt));

histogram(hist_clf_output_t),hold on,histogram(hist_clf_output_nt);
legend('Target','NonTarget')
title(sprintf('auc=%f,pFisher=%f',auc,pseudoFisher));
saveas(gcf,[save_path 'clf_out_hist_pics/' parameters_string '_hist.png']);
close;
save([save_path 'clf_out_hist_data/'  parameters_string '_clf_out.mat'], ...
    'hist_clf_output_t','hist_clf_output_nt','target_start','target_end','w_size_time');



tmp_intervals = intervals(~cellfun(@isempty, intervals));
tmp_intervals_rp_mask = intervals_rp_mask(~isnan(intervals_rp_mask));
[ACC_threshold,F1_threshold,hist_F1_threshold] = statistics(save_path,parameters_string, tmp_intervals,tmp_intervals_rp_mask);
sprintf('ACC = %f\n',hist_F1_threshold)
sprintf('F1 = %f\n',F1_threshold)
sprintf('hist_f1 = %f\n',hist_F1_threshold)
fileID = fopen([data_path, 'results/' 'AUCs.txt'],'a');
fprintf(fileID,'%s, AUC = %f, pFisher = %f, Opt_thres=%f \n\r',parameters_string,auc,pseudoFisher,hist_F1_threshold);
fclose(fileID);
% visualise(tmp_intervals,tmp_intervals_rp_mask,'hist_clf_output')
% [ auc ] = customAUC( intervals,intervals_rp_mask);    
end


function [epochs,contain_event,hist_target,hist_non_target,counter] ...
    = process_interval(data,rp_mask,grad,eog,w_size_time,baseline_time,fs,w_step,prin_comp,classifier,rel_thres)
    hist_target = [];
    hist_non_target = [];
    baseline_size = baseline_time*fs/1000;
    w_size = w_size_time*fs/1000;
%     classifierOutput = nan(size(rp_mask));
    epochs=[];
    rp = find(rp_mask == 1);
    if isempty(rp)
        rp_start = size(data,1)+1;
        rp_end = rp_start; % for epochs without RP, for correct calc of dt_before_mov
    else
        rp_start = rp(1);
        rp_end = rp(end);
        movement = find(rp_mask==10);
    end
    
    contain_event = (sum(rp_mask)>0);
    if (sum(rp_mask)>=0)        %sum(rp_mask)>=0 - for all intervals, sum(rp_mask)>0 - for intervals with RP
        for i = 1:w_step:length(data) - w_size - baseline_size + 1
            base_line = data(i:i+baseline_size-1,:);
            epoch_start = i+baseline_size;
            epoch_end = i+baseline_size+w_size-1;
            epoch_data = data(epoch_start:epoch_end,:)-repmat(mean(base_line,1),w_size,1);
            if (is_relevant(epoch_data,grad(epoch_start:epoch_end,:),eog(epoch_start:epoch_end),rel_thres))
                [X] = get_feats(epoch_data,fs, 0, w_size_time);  %arguments is (data,fs,learn_start,learn_end) learn_start,learn_end - start and end of the interval for learning in ms  fs = 200    
                epoch.Q = (X * prin_comp)* classifier;
                if ~isempty(movement)
                    epoch.dt_before_mov = (epoch_end - movement)/fs; %because our rp region of interest ends 0.5s before movement
                else
                    epoch.dt_before_mov=-6; %long before movement
                end
                if (epoch_end >= rp_start) && (epoch_end<=rp_end)   %If Epoch BEFORE RP, we mark it by 0 label, if epoch AFTER rp, we mark it by -1 label
                    epoch.rp = 1;
                    hist_target=[hist_target,epoch.Q];
                    counter.one = counter.one + 1;
                else
                    debug_flag=false;
                    if rp_start > epoch_end
                        debug_flag=true;
                        epoch.rp = 0;
                        hist_non_target=[hist_non_target,epoch.Q];
                        
                    end
                    if epoch_end > rp_end
                        if debug_flag
                            disp('ERROR, two cond')
                        end
                        debug_flag=true;
                        epoch.rp = -1;
                    end
                    if ~debug_flag
                       disp('ERROR,no con') 
                    end
                end          
                epochs = [epochs,epoch];
            end
        end
    end
  
end



function [is_relevant] = is_relevant(baseline_corrected,grad,eog,rel_thres)
    baselineBlow = sum(max(abs(baseline_corrected),[],1) > rel_thres.base_line.thres) ...
        > rel_thres.base_line.num_channels;
    gradBlow = sum(mean(abs(grad),1) > rel_thres.grad.thres) > rel_thres.grad.num_channels;
    eog_blow = max(abs(eog)) < rel_thres.EOG_thres;
    is_relevant = ~(baselineBlow | gradBlow) | eog_blow ;
end


