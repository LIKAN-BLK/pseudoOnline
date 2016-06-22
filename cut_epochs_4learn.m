function [ eegTRP1, eegNT] = cut_epochs_4learn(datapath,fs,w_size_time,baseline_time,start_rp_time,end_rp_time,ch_to_use,rel_thres)
%loaddata Cut epochs frow semi-raw data
window_size = w_size_time * fs/1000;
[data,eog,mask] = loaddata(datapath,ch_to_use);
[~,grad] = gradient(data);
% [irrelevant_mask] = search_irrelevant_data(mask,data);

tstarts = find(mask == 12);
tends = find(mask == 13);
tstarts = tstarts(1:size(tends,2));

rp1_epochs = [];
for i =1:size(tends,2)
    [tmp_rp1_epochs] = ...
        proc_target_interval(data(tstarts(i):tends(i),:),grad(tstarts(i):tends(i),:),mask(tstarts(i):tends(i)),eog(tstarts(i):tends(i)),window_size,fs, ...
        baseline_time,start_rp_time,end_rp_time,rel_thres);
  
    if ~isempty(tmp_rp1_epochs)
        rp1_epochs = cat(3,rp1_epochs,tmp_rp1_epochs);
    end
end

ntstarts = find(mask == 14);
ntends = find(mask == 15);
ntstarts = ntstarts(1:size(ntends,2));

nt_epochs = [];
for i =1:size(ntends,2)
    [tmp_nt_epochs] = ...
        proc_non_target_interval(data(ntstarts(i):ntends(i),:),grad(ntstarts(i):ntends(i),:),...
            eog(ntstarts(i):ntends(i)),window_size,fs, ...
            baseline_time,rel_thres);
  
    if ~isempty(tmp_nt_epochs)
        nt_epochs = cat(3,nt_epochs,tmp_nt_epochs);
    end
end

eegTRP1 = rp1_epochs;
eegNT = nt_epochs(:,:,randperm(size(nt_epochs,3)));

% [params, spec, sens, acc, auc] = test_loading_alg(cat(3,rp1_epochs(:,:,1:250),rp2_epochs),cat(3,nt1,nt2),fs);

end

function [rp1_epochs] = proc_target_interval(interval_data,grad,epoch_mask,eog,w_size,fs,baseline_time,start_rp_time,end_rp_time,rel_thres)
%function for extracting target epochs from intervals with movement     
    bline_width = baseline_time*fs/1000;
    rp1_epochs = [];
    rp_length = 2;%RP length - 2seconds
    if((find(epoch_mask == 10,1) > (rp_length*fs))) % 12+13+10, so interval contain movement,
        %And this movement 3s after beginning of the interval
        %We have only one rp2 epoch in whole interval
        movement = find(epoch_mask==10);      
        start_rp1_data = movement + start_rp_time * fs /1000;
        end_rp1_data = movement + end_rp_time * fs /1000;
        rp1_epochs = cat(3,rp1_epochs,make_epochs(interval_data(start_rp1_data-bline_width+1:end_rp1_data,:), ... %additional window for baseline
            grad(start_rp1_data-bline_width+1:end_rp1_data,:),eog(start_rp1_data-bline_width+1:end_rp1_data),w_size,bline_width,rel_thres));
    end     
end

function [nt_epochs] = proc_non_target_interval(interval_data,grad,eog,w_size,fs,baseline_time,rel_thres)
    
    bline_width = baseline_time*fs/1000;
    start_rel_data = 30*fs; %We will use data 30s from 14 label
    end_rel_data = size(interval_data,1) - (2000 *fs/1000);  %We will use data 2s before movement
    nt_epochs =  make_epochs(interval_data(start_rel_data:end_rel_data,:), ...
        grad(start_rel_data:end_rel_data,:),eog(start_rel_data:end_rel_data),w_size,bline_width,rel_thres);
end


function [epochs] = make_epochs(data,grad,eog,w_size,bline_width,rel_thres) 
    epochs= [];
    for i = size(data,1):-w_size:(bline_width+w_size) %last window for baseline
        tmp_epoch = data(i-w_size+1:i,:);
        baseline = mean(data(i-(bline_width+w_size) + 1:i-w_size,:),1);
        bcorrected_epoch = tmp_epoch - repmat(baseline,size(tmp_epoch,1),1);
        if ~isempty(tmp_epoch)
            if is_relevant(bcorrected_epoch,eog(i-(bline_width+w_size) + 1:i),grad(i-(bline_width+w_size) + 1:i,:),rel_thres)
                epochs = cat(3,epochs,bcorrected_epoch); 
            end
        end
    end
end

% function [is_relevant] = is_relevant(grad,baseline_corrected,eog,rel_thres)
%     baseline_blow = sum(max(abs(baseline_corrected),[],1) > 70) > 3;
%     grad_blow = sum(mean(abs(grad),1) > 2) > 2;
%     eog_blow = max(abs(eog)) < 200;
%     is_relevant = ~(baseline_blow | grad_blow) | eog_blow; 
% end

function [is_relevant] = is_relevant(baseline_corrected,grad,eog,rel_thres)
    baseline_blow = sum(max(abs(baseline_corrected),[],1) > rel_thres.base_line.thres) ...
        > rel_thres.base_line.num_channels;
    grad_blow = sum(mean(abs(grad),1) > rel_thres.grad.thres) > rel_thres.grad.num_channels;
    eog_blow = max(abs(eog)) < rel_thres.EOG_thres;
    is_relevant = ~(baseline_blow | grad_blow) | eog_blow ;
end
