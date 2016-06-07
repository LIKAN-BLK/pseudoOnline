function [] = pseudoOnlineExp3()
data_path = '../exp3/data/';
fs = 200; %Hz
baseline_time=200; %200ms
w_size_time = 400;%400ms
w_size = w_size_time * fs/1000;
rp_start = -1500;
rp_end=-800;
[ eegTRP,eegNT] = loaddata(data_path,fs,w_size,baseline_time,rp_start,rp_end);

[prin_comp,classifier, opt_thr] = learn(eegTRP(:,:,1:200),eegNT(:,:,1:500),w_size,w_size_time,fs);


data = load([data_path 'rawdata.mat']);
data = (data.eegdata)';
mask = load([data_path 'relmaska.mat']);
mask = mask.relmaska;

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
for i=1:size(ends,2)
   interval = data(starts(i):ends(i),:);
   movement = find(mask(starts(i):ends(i)) == 10);
   if ((movement - rp_len +1) < 0)
       continue
   end
   rp_start = max(1,movement - rp_len + 1);
   
   rp_indexes = [rp_start:movement];
   rp_mask = zeros(size(mask(starts(i):ends(i))));
   rp_mask(rp_indexes) = 1;

       
       
   [clfOut,intervals{i},intervals_rp_mask(i)] = ...
       process_interval(interval,rp_mask,grad(starts(i):ends(i),:),w_size_time,fs,w_step,prin_comp,classifier);
    
%    clf_mv_ind = find(rp_mask == 1);
%    
%    if(~isempty(clf_mv_ind))
%        clf_mv_ind=clf_mv_ind(end);
%        for thres_ind=1:length(threshold1)
%             activation_index1 = find(clfOut1 < threshold1(thres_ind),1);
%             distance_mv1 = clf_mv_ind - activation_index1;
% 
%             activation_index2 = find(clfOut2 < threshold2(thres_ind),1);
%             distance_mv2 = clf_mv_ind - activation_index2;
%             hist1(thres_ind, hist_mv_ind - distance_mv1) = hist1(thres_ind, hist_mv_ind - distance_mv1) + 1;
%             hist2(thres_ind, hist_mv_ind - distance_mv2) = hist2(thres_ind, hist_mv_ind - distance_mv2) + 1;
%        end
%    end
end
% time = (1:size(hist1,2))/fs;
% time = time - time(length(time))/2;
% plot(time,hist1(1,:));
intervals = intervals(~cellfun(@isempty, intervals));
intervals_rp_mask = intervals_rp_mask(~isnan(intervals_rp_mask));
visualise(intervals,intervals_rp_mask,'histogram')
[ auc ] = customAUC( intervals,intervals_rp_mask);
end

function [classifierOutput,epochs,contain_event] = process_interval(data,rp_mask,grad,w_size_time,fs,w_step,prin_comp,classifier)
    w_size = w_size_time*fs/1000;
    classifierOutput = nan(size(rp_mask));
    epochs=[];
    rp = find(rp_mask == 1);
    if isempty(rp)
        rp_start = length(data);
        rp_end = 1; % for epochs without RP, for correct calc of dt_before_mov
    else
        rp_start = rp(1);
        rp_end = rp(end);
    end
    
    
    
    contain_event = (sum(rp_mask)>0);
    if (sum(rp_mask)>=0)        %sum(rp_mask)>=0 - for all intervals, sum(rp_mask)>0 - for intervals with RP
        baseline_size = w_size;
        
        for i = 1:w_step:length(data) - w_size - baseline_size + 1
            base_line = data(i:i+baseline_size-1,:);
            epoch_start = i+baseline_size;
            epoch_end = i+baseline_size+w_size-1;
            epoch.data = data(epoch_start:epoch_end,:)-repmat(mean(base_line,1),w_size,1);
            if (is_relevant(epoch.data,grad(epoch_start:epoch_end,:)))
                [X] = get_feats(epoch.data,200, 0, w_size_time);  %arguments is (data,fs,learn_start,learn_end) learn_start,learn_end - start and end of the interval for learning in ms  fs = 200    
                if all(rp_mask(epoch_start:epoch_end))   %If Epoch BEFORE RP, we mark it by 0 label, if epoch AFTER rp, we mark it by -1 label
                    epoch.rp = 1.0;
                else
                    if rp_start < (epoch_start-w_size) 
                        epoch.rp = 0.0;
                    else
                        epoch.rp = -1.0;
                    end
                end
                epoch.dt_before_mov = (epoch_end - rp_end)/fs;          
               
                epoch.Q = (X * prin_comp)* classifier;
                epochs = [epochs,epoch];
                
            end
        end
    end
  
end



function [is_relevant] = is_relevant(baseline_corrected,grad)
    baselineBlow = sum(max(abs(baseline_corrected),[],1) > 70) > 3;
    gradBlow = sum(mean(abs(grad),1) > 2) > 2;
    is_relevant = ~(baselineBlow | gradBlow);
    is_relevant = true;
end
