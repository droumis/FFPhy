

load times

disp('Loaded Times and Generating Ranges Field');
for n=2:6       
    currnames=names{n};
    [T,R]=strtok(currnames,'-');
    ranges(n,:)=timetrans({T(end-7:end) R(2:end)}, 10000, 2);
end


for i=1:6, names{i}, end

save times names ranges
disp('Saved Times File with the above "names" fields');



