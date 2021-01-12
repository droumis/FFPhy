function o = getFieldIndex(allfields, fields2get)
f = split(allfields, ' ');
o = [];
for g = 1:length(fields2get)
    o = [o find(ismember(f, fields2get{g}))];
end
end