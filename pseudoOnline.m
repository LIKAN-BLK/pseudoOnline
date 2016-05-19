function [] = pseudoOnline()
data_path = '../exp3/data/';
fs = 200; %Hz
[W, mu, sigma,opt_thr] = learn(data_path,learn_start,learn_end,base_start,base_end);

data = load([data_path 'rawdata.mat']);
data = (data.eegdata)';
mask = load([data_path 'relmaska.mat']);
mask = mask.relmaska;


new_mask = zeros(size(mask));
for j = find(mask) 
    new_mask(j-1000:j)=1;
end
plot(new_mask)


w_step = 0.01 * fs;
w_length = (learn_end - base_start)*fs;
% res = zeros(size(mask));
for i = 1:w_step:(size(data,1) - w_length)
        [Xtst] = get_feats(data(i:(i+w_length-1),:), fs, learn_start-base_start, learn_end-base_start,0,base_end - base_start);
        res(round(i+w_length-1):round(i+w_length-1)+w_step-1) = Xtst * W; %res filled with current classifier output until we get output from next window    
end
plot((1:size(res,2))*2,res),hold on, plot((1:size(mask,2))*2,mask);
end

function res = is_relevant(data,mu,sigma)
   res = (all(max(data,[],1)' <= mu+3*sigma,1)) && (all(min(data,[],1)' >= mu-3*sigma,1));        
end