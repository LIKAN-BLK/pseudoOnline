function [] = pseudoOnlineExp3()
data_path = '../exp3/data/';
fs = 200; %Hz
w_size = 400 * fs/1000; %400ms
[ eegTRP1,eegTRP2, eegNT] = loaddata(data_path,fs,w_size);

[W1, ~, ~,opt_thr] = learn(eegTRP1,cat(3,eegTRP2,eegNT),w_size);
[W2, ~, ~,opt_thr] = learn(eegTRP2,cat(3,eegTRP1,eegNT),w_size);

data = load([data_path 'rawdata.mat']);
data = (data.eegdata)';
mask = load([data_path 'relmaska.mat']);
mask = mask.relmaska;

[~,grad] = gradient(data);

rp1_len = 1500*fs/1000;
rp2_len = 500*fs/1000;

starts = find(mask == 12);
ends = find(mask == 13);

w_step = 10 * fs/1000; %10ms


% res = zeros(size(mask));
for i=1:size(ends,2)
   interval = data(starts(i),ends(i));
   movement = find(mask(starts(i),ends(i)));
   if ((movement - rp1_len - rp2_len +1) < 0)
    continue
   end
   rp1_start = max(1,movement - rp1_len - rp2_len + 1);
   rp2_start = max(1,movement - rp2_len + 1);
   rp1_indexes = [rp1_start:rp1_start+rp1_len-1];
   rp2_indexes = [rp2_start:movement];
   rp_mask = zeros(mask(starts(i),ends(i)));
   rp_mask(rp1_indexes) = 1;
   rp_mask(rp2_indexes) = 2;
       
       
   epochs(i,:) = process_interval(interval,rp_mask,grad(starts(i),ends(i)),w_size,w_step,W1,W2);
   
   
end

end

function [epochs] = process_interval(data,rp_mask,grad,w_size,w_step,W1,W2)
    baseline_size = w_size;
    epochs = [];
    for i = 1:w_step:length(data) - w_size - baseline_size + 1
        base_line = data(i:i+baseline_size-1,:);
        epoch_start = i+baseline_size;
        epoch_end = i+baseline_size+w_size-1;
        epoch.data = data(epoch_start:epoch_end,:)-repmat(mean(base_line,1),w_size,1);
        
        if(is_relevant(epoch.data,grad(epoch_start:epoch_end)))
            [X] = get_feats(epoch.data,w_size);        
            epoch.rp1 = sum(rp_mask(epoch_start:epoch_end,:)== 1) / w_size;
            epoch.rp2 = sum(rp_mask(epoch_start:epoch_end,:)== 2) / w_size;
            epoch.Q1 = X * W1;
            epoch.Q2 = X * W2;
            epochs = [epochs,epoch];
        end       
    end    
end

function [is_relevant] = is_relevant(baseline_corrected,grad)
    baselineBlow = sum(max(abs(baseline_corrected),[],1) > 70) > 3;
    gradBlow = sum(mean(abs(grad),1) > 2) > 2;
    is_relevant = ~(baselineBlow | gradBlow); 
end
