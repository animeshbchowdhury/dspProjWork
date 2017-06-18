function [g1,d1,d2] = checkValidity_MultiFrquency(X1,Y1,f1,NFFT1)
d1 = 0;
d2 = 0;

k = 1;
g2 = 0;
% j = 1;
% flag1 = 0;
peak1 = sort(findpeaks(2*X1(1:NFFT1/2+1)),'descend');
for i = 1:length(peak1)
chk_val1(i,:) = peak1(1,i)/mean(peak1(i+1:length(peak1)));
    if(chk_val1(i,:) > 2)
        index1(i,:) = find(X1(1:NFFT1/2+1) == peak1(:,i)/2);
        d1(1,k) = f1(index1(i,:));
        k = k+1;
    end
end
mean1 = mean(peak1(1:length(peak1)));
% mean = mean(2*X1(1:NFFT1/2+1));

% for i = 1:length(peak1)
%     ratio1(i,:) = peak1(i)/mean1;
%     if(ratio1(i,:) > 2)
%         rat1(j,:) = ratio1(i,:);
%         j = j+1;
%         flag1 = 1;
%     end
% end

 k = 1;
% flag2 = 0;
peak2 = sort(findpeaks(2*Y1(1:NFFT1/2+1)),'descend');
for i = 1:length(peak2)-1
    chk_val(i,:) = peak2(1,i)/mean(peak2(i+1:length(peak2)));
    if(chk_val(i,:) > 3)
        g2(:,k) = chk_val(i,:);
        index2(i,:) = find(Y1(1:NFFT1/2+1) == peak2(:,i)/2);
        d2(1,k) = f1(index2(i,:));
        k = k+1;
    end
end
mean2 = mean(peak2(1:length(peak2)));
% mean2(1,1:NFFT1/2+1) = 2*Y1(1,1:NFFT1/2+1);
% mean(mean2)
ratio2 = peak2(1,1)/mean2;
% if(length(d1) == 0)
%     d1 = 0;
% end
% if(length(d2) == 0)
%     d2 = 0;
% end

if(g2 >= 0)
    g1 = mean(g2);
else
    g1 = mean(g2);
end
    

% for i = 1:length(peak2)
%     ratio2(i,:) = peak2(i)/mean2;
%     if(ratio2 > 2)
%         rat2(k,:) = ratio2(i,:);
%         k = k+1;
%         flag2 = 1;
%     end
% end

% if((flag1 == 0) && (flag2 == 0))
%     g1 = 0;
% 
% elseif((flag1 == 1) && (flag2 == 1))
%     for i = 1:max(length(rat1),length(rat2))
%         if (max(length(rat1),length(rat2)) == length(rat1))
%             g1(i,:) = rat1(i,:);
%         else
%             g1(i,:) = rat2(i,:);
%         end
% % g1(i,:) = rat1(i,:)/min(length(rat1),length(rat2)) + rat2(i,:)/min(length(rat1),length(rat2));
%     end
%     
% elseif((flag1 == 1) && (flag2 == 0))
%     rat2(k,:) = 0;
%     for i = 1:max(length(rat1),length(rat2))
%         if (max(length(rat1),length(rat2)) == length(rat1))
%             g1(i,:) = rat1(i,:);
%         else
%             g1(i,:) = rat2(i,:);
%         end
% % g1(i,:) = rat1(i,:)/min(length(rat1),length(rat2)) + rat2(i,:)/min(length(rat1),length(rat2));
%     end
% 
% 
% else
%     rat1(k,:) = 0;
%     for i = 1:max(length(rat1),length(rat2))
%         if (max(length(rat1),length(rat2)) == length(rat1))
%             g1(i,:) = rat1(i,:);
%         else
%             g1(i,:) = rat2(i,:);
%         end
% % g1(i,:) = rat1(i,:)/min(length(rat1),length(rat2)) + rat2(i,:)/min(length(rat1),length(rat2));
%     end
% end
% 



% for i = 1:max(length(rat1),length(rat2))
%     if (max(length(rat1),length(rat2)) == length(rat1))
%         g1(i,:) = rat1(i,:);
%     else
%         g1(i,:) = rat2(i,:);
%     end
% % g1(i,:) = rat1(i,:)/min(length(rat1),length(rat2)) + rat2(i,:)/min(length(rat1),length(rat2));
% end
