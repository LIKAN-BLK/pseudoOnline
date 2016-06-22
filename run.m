function [] = run()

%Clean results directory
data_path='../exp5/';
d = dir([data_path,'results/']); %delete all content from results/ dir
ifile = ~[d(:).isdir];
delete([data_path,'results/' d(ifile).name]);


isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i=1:size(nameFolds,1)
    dir2rm = fullfile([data_path,'results/'],nameFolds{i});
    rmdir(dir2rm, 's');
end



%DEFAULTS
rp_start = -2000;
rp_end = -500;
w_size_time = 1000;%400ms
ch_to_use.EEGch = [1:56];
ch_to_use.EOGch = [57,58];
rel_thres.base_line.thres = 70;
rel_thres.base_line.num_channels = 3;
rel_thres.grad.thres=2;
rel_thres.grad.num_channels=2;
rel_thres.EOG_thres=200;

pseudoOnlineExp3(data_path,w_size_time,rp_start,rp_end,ch_to_use,rel_thres)

% w_size_time = 400;%400ms
% for rp_end = -2000+w_size_time:100:-500
%     for rp_start = -2000:100:rp_end - w_size_time 
%         disp(sprintf('rp_start = %d,rp_end=%d\n\r',rp_start,rp_end))
%         pseudoOnlineExp3(data_path,w_size_time,rp_start,rp_end)
%     end
% end
% end

