function [data,EOG,mask] = loaddata(data_path)
%load_data Get data, eog and mask from files.
if isempty(data_path)
    data_path = '../exp4/';
end
mask = load([data_path 'relmaska.mat']);
mask = mask.relmaska;
if exist([data_path 'rawdataAllChnl.mat'],'file') ==2
    data = load([data_path 'rawdataAllChnl.mat']);
    EOG = data.eegdata(37:38,:);
    EOG = EOG';
    EOG = EOG(:,1)-EOG(:,2);
else
    EOG = zeros(size(mask));
end

data = load([data_path 'rawdata.mat']);
data = (data.eegdata)';


end

