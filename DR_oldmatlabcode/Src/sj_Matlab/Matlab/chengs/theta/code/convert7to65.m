function convert7to65

% find all files
f= dir('behav*.mat');
nfiles=size(f,1);
[files{1:nfiles}]= deal(f.name);

for i=1:nfiles
    load(f(i).name);
    save(f(i).name, 'behavdata', '-v6');
end
