% 
% %combine the spat info of all trajs within each ep for each cell
% %create single num index
clear tsi fu tmpnum spatmodinfo_alltraj dta tsi c ia ic spatmodinfo_alltraj
% dta = spatinfo{1}.output{1};
% tsi = spatinfo{1}.output{1}(:,1:4);
% for fu = 1:length(tsi);
% tmpnum(fu,1) = [tsi(fu,1)*1000+tsi(fu,2)*100+tsi(fu,3)*10+tsi(fu,4)];

%combine the spat info of all trajs within each ep for each cell

dta = spatinfo{1}.output{1};
tmpnum = nan(length(dta(:,1)),1);
for hu = 1:length(dta(:,1));
    x = dta(hu,1:4);
    y = sprintf('%d',x);
    z = str2num(y);
    tmpnum(hu) = z;
end
[c ia ic] = unique(tmpnum, 'stable');
spatmodinfo_alltraj = [];
for gu = 1:length(ia);
    if (ia(gu) < length(ia)); %if it's not the last one
        spatmodinfo_alltraj(gu,5) = nansum(dta(ia(gu):(ia(gu+1)-1),6));
        spatmodinfo_alltraj(gu,1:4) = dta(ia(gu),1:4); %take the index of first spt info number
    else
        spatmodinfo_alltraj(gu,5) = nansum(dta(ia(gu):end,6));
        spatmodinfo_alltraj(gu,1:4) = dta(ia(gu),1:4);
    end
end



%add back the indices
% for bu = 1:length(c);
%     %rebuild combined index
%      spatmodinfo_alltraj(bu,4) = rem(c(bu),10);
%      spatmodinfo_alltraj(bu,3) = (rem(c(bu),100)-rem(c(bu),10))/10;
%      spatmodinfo_alltraj(bu,2) = ((rem(c(bu),1000))-rem(c(bu),100))/100;
%      spatmodinfo_alltraj(bu,1) = ((rem(c(bu),10000))-rem(c(bu),1000))/1000;
%     spatmodinfo_alltraj(bu,5) = tmpsum(bu,1);
% end
  
    

