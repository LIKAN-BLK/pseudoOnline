function [] = pseudoOnlineExp3()
data_path = '../exp3/data/';
fs = 200; %Hz
w_size = 400 * fs/1000; %400ms
[ eegTRP1,eegTRP2, eegNT] = loaddata(data_path,fs,w_size);

[prin_comp1,classifier1, opt_thr1] = learn(cat(3,eegTRP2(:,:,1:200),eegTRP1(:,:,1:200)),eegNT(:,:,1:500),w_size);
[prin_comp2,classifier2,opt_thr2] = learn(eegTRP2(:,:,1:200),eegTRP1(:,:,1:400),w_size);

data = load([data_path 'rawdata.mat']);
data = (data.eegdata)';
mask = load([data_path 'relmaska.mat']);
mask = mask.relmaska;

[~,grad] = gradient(data);

rp1_len = 1500*fs/1000;
rp2_len = 500*fs/1000;
rp_len = 2000 * fs /1000;

starts = find(mask == 12);
ends = find(mask == 13);

threshold1= 0:-1:-5;
threshold2= 0:-1:-5;
hist1 = zeros(length(threshold1),2*round(ends(1) - starts(1))); %because rp could be at the beginning of interval and could be at the end of interval
hist2 = zeros(length(threshold2),2*round(ends(1) - starts(1)));
hist_mv_ind = round(size(hist1,2)/2);

w_step = 10 * fs/1000; %10ms

classifierOutput1 = {};
classifierOutput2 = {};
% res = zeros(size(mask));
for i=1:size(ends,2)
   interval = data(starts(i):ends(i),:);
   movement = find(mask(starts(i):ends(i)) == 10);
   if ((movement - rp_len +1) < 0)
       continue
   end
   rp1_start = max(1,movement - rp_len + 1);
   rp2_start = max(1,movement - rp2_len + 1);
   rp1_indexes = [rp1_start:rp1_start+rp1_len-1];
   rp2_indexes = [rp2_start:movement];
   rp_mask = zeros(size(mask(starts(i):ends(i))));
   rp_mask(rp1_indexes) = 1;
   rp_mask(rp2_indexes) = 2;
       
       
   [clfOut1,clfOut2,intervals{i}] = process_interval(interval,rp_mask,grad(starts(i):ends(i),:),w_size,w_step,prin_comp1,prin_comp2,classifier1,classifier2);
    
%    clf_mv_ind = find(rp_mask == 2);
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
Q=[];
Labels = [];
for i = 1:length(intervals)
    for epoch=intervals{i}
        if ~isnan(epoch.Q1)
            Q = [Q,epoch.Q1];
            Labels = [Labels,epoch.rp];
        end
    end
end
[~,~,~,auc,~,~,~] = perfcurve(Labels,Q,1);
end

function [classifierOutput1,classifierOutput2,epochs] = process_interval(data,rp_mask,grad,w_size,w_step,prin_comp1,prin_comp2,classifier1,classifier2)
    classifierOutput1 = nan(size(rp_mask));
    classifierOutput2 = nan(size(rp_mask));
    epochs = [];
    if(sum(rp_mask)>0)
        baseline_size = w_size;
        
        for i = 1:w_step:length(data) - w_size - baseline_size + 1
            base_line = data(i:i+baseline_size-1,:);
            epoch_start = i+baseline_size;
            epoch_end = i+baseline_size+w_size-1;
            epoch.data = data(epoch_start:epoch_end,:)-repmat(mean(base_line,1),w_size,1);

            [X] = get_feats(epoch.data,200, 0, 400);  %arguments is (data,fs,learn_start,learn_end) learn_start,learn_end - start and end of the interval for learning in ms      
            epoch.rp = all(rp_mask(epoch_start:epoch_end));           
            
            if (is_relevant(epoch.data,grad(epoch_start:epoch_end,:)))
                classifierOutput1(epoch_end:epoch_end+w_step-1) = (X *prin_comp1)* classifier1;
                classifierOutput2(epoch_end:epoch_end + w_step-1) = (X * prin_comp2) * classifier2;
            else
                classifierOutput1(epoch_end:epoch_end+w_step-1) = nan;
                classifierOutput2(epoch_end:epoch_end+w_step-1) = nan;
            end      
            
            
            if(is_relevant(epoch.data,grad(epoch_start:epoch_end,:)))
                epoch.Q1 = (X * prin_comp1)* classifier1;
                epoch.Q2 = (X * prin_comp2) * classifier2;
            else
                epoch.Q1 = nan;
                epoch.Q2 = nan;
            end
            epoch.is_relevant = is_relevant(epoch.data,grad(epoch_start:epoch_end,:));
%             if(epoch.rp)
            epochs = [epochs,epoch];
%             end
        end
    end
   
end



function [is_relevant] = is_relevant(baseline_corrected,grad)
    baselineBlow = sum(max(abs(baseline_corrected),[],1) > 70) > 3;
    gradBlow = sum(mean(abs(grad),1) > 2) > 2;
    is_relevant = ~(baselineBlow | gradBlow);  
end
