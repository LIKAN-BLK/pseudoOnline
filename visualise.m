% 
% output1 = [];
% for i=1:size(classifierOutput1,2)
%     if (size(classifierOutput1{i},2)==322)
%         output1 = [output1;classifierOutput1{i}];
%     end
% end
% 
% output2 = [];
% for i=1:size(classifierOutput2,2)
%     if (size(classifierOutput2{i},2)==322)
%         output2 = [output2;classifierOutput2{i}];
%     end
% end
% time =(1:size(output1,2))/200 - size(output1,2)/200;
% plot(time,mean(output1,1))
% plot(time,mean(output1,1))
% 
% cumm1=zeros(size(output1));
% for i=1:5:size(output1,2)-5
%     cumm1(:,i:i+5-1)=repmat(sum(output1(:,i:i+5-1),2),1,5);
% end
% surf(time,1:size(cumm1,1),cumm1);
% 
% 
% surf(time,1:size(output1,1),output1);

%  tmp_start = find(~isnan(clfOut1));
%    if(~isempty(tmp_start))
%     tmp_end = tmp_start(end);
%     tmp_start = tmp_start(1);
%     if(all(~isnan(clfOut1(tmp_start:tmp_end))))
%         classifierOutput1{i} = clfOut1(tmp_start:tmp_end);
%     end
%    end
%    
%    if(~isempty(tmp_start))
%     tmp_start = find(~isnan(clfOut2));
%     tmp_end = tmp_start(end);
%     tmp_start = tmp_start(1);
%     if(all(~isnan(clfOut2(tmp_start:tmp_end))))
%         classifierOutput2{i} = clfOut2(tmp_start:tmp_end);
%     end
%    end




