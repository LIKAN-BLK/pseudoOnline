function [EEG,EOG,mask] = loaddata(data_path,ch_to_use)
%load_data Get data, eog and mask from files.
mask = load([data_path 'maska.mat']);
mask = mask.maska;


data = load([data_path 'rawdata.mat']);
EEG = (data.eegdata(ch_to_use.EEGch,:))';
EOG = (data.eegdata(ch_to_use.EOGch(1),:) - data.eegdata(ch_to_use.EOGch(2),:))';
end

