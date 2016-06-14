function [] = run()
data_path='../exp4/';
%DEFAULTS
% save_path=[data_path 'results/'];
% rp_start = -1500;
% rp_end = -800;
% w_size_time = 400;%400ms
% pseudoOnlineExp3(data_path,save_path,w_size_time,rp_start,rp_end)

w_size_time = 400;%400ms
for rp_end = -2000+w_size_time:100:-500
    for rp_start = -2000:100:rp_end - w_size_time 
        sprintf('rp_end=%f,rp_start=%f\n\r',rp_end,rp_start)
        save_path=[data_path 's_' num2str(rp_start) 'e_' num2str(rp_end) 'w_' num2str(w_size_time) '/'];
        mkdir(save_path);
        pseudoOnlineExp3(data_path,save_path,w_size_time,rp_start,rp_end)
    end
end
end

