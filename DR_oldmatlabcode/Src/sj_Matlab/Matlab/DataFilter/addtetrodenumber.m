function addtetrodelocation(animdirect,fileprefix)
% This function adds tetrode number to the cellinfo structure


load([animdirect, fileprefix,'cellinfo']);
o = cellfetch(cellinfo,'numspikes');

for i = 1:size(o.index, 1);
    in = o.index(i,:);
    cellinfo{in(1)}{in(2)}{in(3)}{in(4)}.tetnum = in(3);
end

save([animdirect, fileprefix,'cellinfo'], 'cellinfo')

