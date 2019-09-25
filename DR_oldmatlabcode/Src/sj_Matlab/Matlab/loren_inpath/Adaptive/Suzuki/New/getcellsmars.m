%cd d:\Data\;
cd c:\matlab_stuff\cindy\cindymars\data;

d = dir;
ncells = length(dir) - 2;
dataset = cell(ncells,1);
for i = 1:ncells,
    dataset{i}.name = d(i+2).name;
    dataset{i}.cobj = lorencortex(dataset{i}.name, [1,2]);
end
