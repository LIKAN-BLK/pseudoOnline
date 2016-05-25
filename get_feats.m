function [ X ] = get_feats(eeg, fs, learn_start, learn_end)
%Get ampl features, learn_start, learn_end in ms
    learn_start = learn_start/1000;
    learn_end = learn_end/1000;

    desired_feature_num = 19; 
    full_window = learn_end - learn_start;

       
    filter_width = full_window/(0.5*desired_feature_num + 1);
    filter_step = 0.5 * filter_width;

    
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