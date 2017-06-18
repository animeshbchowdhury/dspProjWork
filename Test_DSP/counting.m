function counting(frequency1,frequency2)
q = vertcat(frequency1,frequency2);
maxVal = unique(q);
out = [maxVal,histc(q(:),maxVal)];
% [N,edges] = histcounts(out,75)
% max(out(:,2));
% c = 0;
% y = size(q);
% for i = 1:y(:,1)
%     for j = 1:y(:,2)
%     if(q(i,j) == 100)
%         c = c+1;
%     end
%     end
% end

%mode(q)
minbin = min(out(:,1))+2;
maxbin = max(out(:,1));

frequencyRep = 0;
k = 1;
binNo = 2;
while(maxbin > minbin)
for i = 1:length(out)
    if(out(i,1)>= minbin && out(i,1)<= minbin + binNo)
        frequencyRep(k,1) = frequencyRep(k,1) + out(i,2);
        frequencyRep(k,2) = (2*minbin+binNo)/2;
    end
end
k = k+1;
minbin = minbin+binNo;
frequencyRep(k,1) = 0;
end


cell1 = frequencyRep(1:end,1);
Newcell1 = cell1(cell1 ~= 0);

cell2 = frequencyRep(1:end,2);
Newcell2 = cell2(cell2 ~= 0);

NewfrequencyRep = [Newcell1 Newcell2];
% f=fit(frequencyRep(1:end,2),frequencyRep(1:end,1),'smoothingspline');
plot(NewfrequencyRep(1:end,2),NewfrequencyRep(1:end,1))


% clf;
% figure(1)
% plot(frequencyRep(1:end,2),frequencyRep(1:end,1));
% hold on

% [~,locs] = findpeaks(frequencyRep(1:end,1));
% plot(frequencyRep(locs(),2),frequencyRep(locs,1),'rs','MarkerFaceColor','b')
% hold on

% for p = 1:length(locs)
% % index1(p,:) = find(frequencyRep(1:end/2,2) == locs);
%     d1(p,1) = frequencyRep(locs(p,1),2);
% 
% end

f=fit(NewfrequencyRep(1:end,2),NewfrequencyRep(1:end,1),'smoothingspline');
figure(1)
plot(NewfrequencyRep(1:end,2),NewfrequencyRep(1:end,1));
figure(2)
plot(f,NewfrequencyRep(1:end,2),NewfrequencyRep(1:end,1));
% hold on
% plot(frequencyRep(locs(),2),frequencyRep(locs,1),'rs','MarkerFaceColor','b')
% hold on
% text(frequencyRep(locs,2),frequencyRep(locs,1),num2str(d1(:,1)),'HorizontalAlignment','right');

xlabel('Frequency (in Hz)')
ylabel('Amplitude')
title('SHAKE METER: Autonomous Vibration Detection System')


% frequencyRep = 0;
% 
% for i = 1:length(out)
%     if(out(i,1)>= 148 && out(i,1)<= 152)
%         frequencyRep = frequencyRep + out(i,2);
%     end
% end
% disp('range')
% disp(frequencyRep)