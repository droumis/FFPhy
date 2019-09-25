cd Data;
d = dir;
dataset = cell(33,1);
for i = 1:33,
    dataset{i}.name = d(i+2).name;
    dataset{i}.c{1} = cortex(dataset{i}.name, [1,2]);
end


