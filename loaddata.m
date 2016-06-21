function [EEG,EOG,mask] = loaddata(data_path,EEGch,EOGch)
%load_data Get data, eog and mask from files.
if isempty(data_path)
    data_path = '../exp5/';
end
mask = load([data_path 'maska.mat']);
mask = mask.maska;


data = load([data_path 'rawdata.mat']);
EEG = (data.eegdata(EEGch,:))';
EOG = (data.eegdata(EOGch(1),:) - data.eegdata(EOGch(2),:))';
end

