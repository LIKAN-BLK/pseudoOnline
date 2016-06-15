function [] = run()
data_path='../exp4/';
d = dir([data_path,'results/']); %delete all content from results/ dir
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i=1:size(nameFolds,1)
    dir2rm = fullfile([data_path,'results/'],nameFolds{i});
    rmdir(dir2rm, 's');
end
%DEFAULTS
rp_start = -1700;
rp_end = -500;
w_size_time = 400;%400ms
pseudoOnlineExp3(data_path,w_size_time,rp_start,rp_end)

% w_size_time = 400;%400ms
% for rp_end = -2000+w_size_time:100:-500
%     for rp_start = -2000:100:rp_end - w_size_time 
%         disp(sprintf('rp_start = %d,rp_end=%d\n\r',rp_start,rp_end))
%         pseudoOnlineExp3(data_path,w_size_time,rp_start,rp_end)
%     end
% end
% end

