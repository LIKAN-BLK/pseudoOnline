% function [ X ] = get_feats(eeg, w_size)
% %Get ampl features
%     desired_feature_num = 15; 
%     
%     filter_width = round(w_size/desired_feature_num); 
%     filter_step = round(0.5 * filter_width);
%     
%     % features
%     
%     ind_beg = [1:filter_step:w_size-filter_width];
%     ind_end = ind_beg + filter_width;
%     
%     
%     N = size(eeg, 3);
%     X = [];
%     for i = 1:N
%         x = [];
%         for t = 1:length(ind_beg)    
%             x = [x; mean(eeg(ind_beg(t):ind_end(t), :, i), 1)];
%         end
%         X(i, :) = x(:);
%     end
% end
function [ X ] = get_feats(eeg, fs, learn_start, learn_end,base_start,base_end,need_preprocess)
%Get ampl features
    if (need_preprocess)
        eeg = preprocess(eeg, fs,base_start,base_end);
    end
    desired_feature_num = 15; 
    full_window = learn_end - learn_start;
    % feat_num = number_of_steps;
    % filter_step = (window_size - filter_width)/num_of_steps
    % filter_step = 0.5*filter_width
    
    filter_width = 2*full_window/(desired_feature_num + 2);
    filter_step = 0.5 * filter_width; 
%     filter_width = roundn(full_window/desired_feature_num,-2); 
%     filter_step = 0.5 * filter_width;
    
    % features
    
    times_beg = [learn_start:filter_step:learn_end-filter_width];
    times_end = times_beg + filter_width;
    ts_beg = round((times_beg) .* fs);
    if (ts_beg(1) == 0)
        ts_beg(1) = 1;
    end
    ts_end = round((times_end) .* fs);
    
    N = size(eeg, 3);
    X = [];
    for i = 1:N
        x = [];
        for t = 1:length(ts_beg)    
            x = [x; mean(eeg(ts_beg(t):ts_end(t), :, i), 1)];
        end
        X(i, :) = x(:);
    end


end
function [eeg] = preprocess(eeg, fs,base_start,base_end)
% baseline correction
base_end_time = base_end;
base_start_time = base_start;
base_end = (base_end_time) * fs;
base_start = (base_start_time) * fs+1;
eeg = eeg - repmat(mean(eeg(base_start:base_end,:,:),1), [size(eeg,1), 1, 1]);
end


