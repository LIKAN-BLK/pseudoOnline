function [] = pseudoOnlineExp3(data_path,save_path,w_size_time,rp_start,rp_end)

fs = 200; %Hz
baseline_time=200; %200ms

w_size = w_size_time * fs/1000;
[eegTRP,eegNT] = cut_epochs_4learn([],fs,w_size,baseline_time,rp_start,rp_end);
[prin_comp,classifier, opt_thr] = learn(eegTRP(:,:,1:200),eegNT(:,:,1:500),w_size,w_size_time,fs);


[data,EOG,mask] = loaddata(data_path);

[~,grad] = gradient(data);


rp_len = 2000 * fs /1000;

starts = find(mask == 12);
ends = find(mask == 13);

% threshold1= 0:-1:-5;
% threshold2= 0:-1:-5;
% hist1 = zeros(length(threshold1),2*round(ends(1) - starts(1))); %because rp could be at the beginning of interval and could be at the end of interval
% hist2 = zeros(length(threshold2),2*round(ends(1) - starts(1)));
% hist_mv_ind = round(size(hist1,2)/2);

w_step = 10 * fs/1000; %10ms

intervals_rp_mask = nan(1,length(ends));
hist_target=[];
hist_non_target=[];
hist_non_target2=[];
counter.one=0;
counter.zero=0;
counter.minus=0;
counter.all=0;
for i=1:size(ends,2)
   interval = data(starts(i):ends(i),:);
   movement = find(mask(starts(i):ends(i)) == 10);
   if ((movement - rp_len +1) < 0)
       continue
   end
   
   rp_mask = zeros(size(mask(starts(i):ends(i))));
   if ~isempty(movement)
    rp_start = max([1,movement - rp_len + 1]);
    rp_indexes = [rp_start:movement];
    rp_mask(rp_indexes) = 1;
   end
   grad_interval = grad(starts(i):ends(i),:);
   EOG_interval = EOG(starts(i):ends(i));
   
       
       
   [intervals{i},intervals_rp_mask(i),tmp_hist_target,tmp_hist_non_target,tmp_hist_non_target2,tmp_counter] = ...
       process_interval(interval,rp_mask,grad_interval,EOG_interval,w_size_time,baseline_time,fs,w_step,prin_comp,classifier);
    
    hist_target = [hist_target,tmp_hist_target];
    hist_non_target = [hist_non_target,tmp_hist_non_target];
    hist_non_target2 = [hist_non_target2,tmp_hist_non_target2];
    
    counter.one=counter.one+tmp_counter.one;
    counter.zero=counter.zero+tmp_counter.zero;
    counter.minus=counter.minus+tmp_counter.minus;
    counter.all=counter.all+tmp_counter.all;
       
end
histogram(hist_target),hold on,histogram(hist_non_target);
legend('Target','NonTarget')
% saveas(gcf,[save_path 'clf_output_hist.png']);
close;

tmp_intervals = intervals(~cellfun(@isempty, intervals));
tmp_intervals_rp_mask = intervals_rp_mask(~isnan(intervals_rp_mask));
visualise(tmp_intervals,tmp_intervals_rp_mask,'hist_clf_output')
% [ auc ] = customAUC( intervals,intervals_rp_mask);
end

function [epochs,contain_event,hist_target,hist_non_target,hist_non_target2,counter] ...
    = process_interval(data,rp_mask,grad,eog,w_size_time,baseline_time,fs,w_step,prin_comp,classifier)
    counter.one=0;
    counter.zero=0;
    counter.minus=0;
    counter.all=0;
    hist_target = [];
    hist_non_target = [];
    hist_non_target2 = [];
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
    end
    
    
    
    contain_event = (sum(rp_mask)>0);
    if (sum(rp_mask)>=0)        %sum(rp_mask)>=0 - for all intervals, sum(rp_mask)>0 - for intervals with RP
        for i = 1:w_step:length(data) - w_size - baseline_size + 1
            base_line = data(i:i+baseline_size-1,:);
            epoch_start = i+baseline_size;
            epoch_end = i+baseline_size+w_size-1;
            epoch_data = data(epoch_start:epoch_end,:)-repmat(mean(base_line,1),w_size,1);
            if (is_relevant(epoch_data,grad(epoch_start:epoch_end,:),eog(epoch_start:epoch_end)))
                [X] = get_feats(epoch_data,fs, 0, w_size_time);  %arguments is (data,fs,learn_start,learn_end) learn_start,learn_end - start and end of the interval for learning in ms  fs = 200    
                epoch.Q = (X * prin_comp)* classifier;
                if all(rp_mask(epoch_start:epoch_end))   %If Epoch BEFORE RP, we mark it by 0 label, if epoch AFTER rp, we mark it by -1 label
                    epoch.rp = 1;
                    hist_target=[hist_target,epoch.Q];
                    counter.one=counter.one+1;
                else
                    debug_flag=false;
                    if rp_start > epoch_start
                        epoch.rp = 0;
                        hist_non_target=[hist_non_target,epoch.Q];
                        counter.zero=counter.zero+1;
                        debug_flag=true;
                    end
                    if epoch_end > rp_end
                        epoch.rp = -1;
                        hist_non_target2=[hist_non_target2,epoch.Q];
                        counter.minus=counter.minus+1;
                        debug_flag=true;
                    end
                    if ~debug_flag
                        disp('error')
                    end
                end
                epoch.dt_before_mov = (epoch_end - rp_end)/fs;          
                epochs = [epochs,epoch];
                counter.all=counter.all+1;
            end
        end
    end
  
end



function [is_relevant] = is_relevant(baseline_corrected,grad,eog)
    baselineBlow = sum(max(abs(baseline_corrected),[],1) > 70) > 3;
    gradBlow = sum(mean(abs(grad),1) > 2) > 2;
    eog_blow = max(abs(eog)) < 200;
    is_relevant = ~(baselineBlow | gradBlow) | eog_blow ;
end
