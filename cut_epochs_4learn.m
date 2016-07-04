function [ eegTRP1, eegNT,intervals] = cut_epochs_4learn(datapath,fs,w_size_time,baseline_time,start_rp_time,end_rp_time,ch_to_use,rel_thres)
%loaddata Cut epochs frow semi-raw data
window_size = w_size_time * fs/1000;
[data,eog,mask] = loaddata(datapath,ch_to_use);
[~,grad] = gradient(data);
% [irrelevant_mask] = search_irrelevant_data(mask,data);

tstarts = find(mask == 12);
tends = find(mask == 13);
tstarts = tstarts(1:size(tends,2));

intervals = {};
target_epochs = [];
for i =1:size(tends,2)
    [tmp_intervals,tmp_target_epochs] = ...
        proc_target_interval(data(tstarts(i):tends(i),:),grad(tstarts(i):tends(i),:),mask(tstarts(i):tends(i)),eog(tstarts(i):tends(i)),window_size,fs, ...
        baseline_time,start_rp_time,end_rp_time,rel_thres);
  
    if ~isempty(tmp_intervals)
        intervals = [intervals,{tmp_intervals}];
        target_epochs = [target_epochs,tmp_target_epochs];
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
        nt_epochs = cat(1,nt_epochs,tmp_nt_epochs);
    end
end

eegTRP1 = cat(1,target_epochs.feats);
eegNT = nt_epochs(randperm(size(nt_epochs,1)),:);

% [params, spec, sens, acc, auc] = test_loading_alg(cat(3,rp1_epochs(:,:,1:250),rp2_epochs),cat(3,nt1,nt2),fs);

end

function [intervals,target_epochs] = proc_target_interval(interval_data,grad,epoch_mask,eog,w_size,fs,baseline_time,start_rp_time,end_rp_time,rel_thres)
%function for extracting target epochs from intervals with movement     
    bline_width = baseline_time*fs/1000;
    intervals = [];
    target_epochs=[];
    rp_length = 2;%RP length - 2seconds
    w_step = 50 * fs/1000; %50ms
    if((find(epoch_mask == 10,1) > (rp_length*fs)))
        
        movement = find(epoch_mask==10);      
        start_rp1_data = movement + start_rp_time * fs /1000;
        end_rp1_data = movement + end_rp_time * fs /1000;
        correct_dt = @(x) setfield(x,'dt_before_movement',x.dt_before_movement/fs + end_rp_time/1000);
        set_rp_true = @(x) setfield(correct_dt(x),'rp',1);
        target_epochs = make_epochs(interval_data(start_rp1_data-bline_width+1:end_rp1_data,:), ... %additional window for baseline
            grad(start_rp1_data-bline_width+1:end_rp1_data,:),...
            eog(start_rp1_data-bline_width+1:end_rp1_data),w_size,w_step,bline_width,rel_thres,set_rp_true);
        set_rp_false = @(x) setfield(correct_dt(x),'rp',0);
        tmp_non_rp_epochs = make_epochs(interval_data(1:start_rp1_data-bline_width,:), ... %additional window for baseline
            grad(1:start_rp1_data-bline_width,:),eog(1:start_rp1_data-bline_width),w_size,w_step,bline_width,rel_thres,set_rp_false);
        
        intervals = [tmp_non_rp_epochs,target_epochs];
    end
end

function [nt_epochs] = proc_non_target_interval(interval_data,grad,eog,w_size,fs,baseline_time,rel_thres)
    w_step = 100 * fs/1000;
    bline_width = baseline_time*fs/1000;
    start_rel_data = 30*fs; %We will use data 30s from 14 label
    end_rel_data = size(interval_data,1) - (2000 *fs/1000);  %We will use data 2s before movement
    empty_func = @(x) x;
    nt_epochs =  make_epochs(interval_data(start_rel_data:end_rel_data,:), ...
        grad(start_rel_data:end_rel_data,:),eog(start_rel_data:end_rel_data),w_size,w_step,bline_width,rel_thres,empty_func);
    nt_epochs = cat(1,nt_epochs.feats);
end


function [epochs] = make_epochs(data,grad,eog,w_size,w_step,bline_width,rel_thres,set_field_func) 
    epochs= [];
    
    for i = size(data,1):-w_step:(bline_width+w_size) %last window for baseline
        tmp_epoch = data(i-w_size+1:i,:);
        baseline = mean(data(i-(bline_width+w_size) + 1:i-w_size,:),1);
        bcorrected_epoch = tmp_epoch - repmat(baseline,size(tmp_epoch,1),1);
        res.dt_before_movement = i-size(data,1);                       %range between epoch end and data end. Needed for target intervals to calc time before movement
        tmp_eog = eog(i-w_size+1:i) - mean(eog(i-(bline_width+w_size) + 1:i-w_size),1);
        
        if is_relevant(bcorrected_epoch,grad(i-(bline_width+w_size) + 1:i,:),tmp_eog,rel_thres)
            res.feats = get_feats(bcorrected_epoch);
            res = set_field_func(res);
            epochs = [epochs,res]; 
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
    baseline_blow = sum(max(abs(baseline_corrected),[],1) > rel_thres.base_line.thres) > ...
         rel_thres.base_line.num_channels;
    grad_blow = sum(mean(abs(grad),1) > rel_thres.grad.thres) > rel_thres.grad.num_channels;
    eog_blow = max(abs(eog)) < rel_thres.EOG_thres;
    is_relevant = (~(baseline_blow | grad_blow)) & eog_blow ;
end
