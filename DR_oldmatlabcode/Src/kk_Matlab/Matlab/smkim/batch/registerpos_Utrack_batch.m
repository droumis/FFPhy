
load('Utrack_template');

l = dir('rawpos*.mat');
for f = 1:length(l)
    infilename = l(f).name;
    load(infilename);
    outfilename = [ 'registeredpos_' infilename(8:end) ];
    disp(outfilename)
    registeredpos = registerpos_Utrack(rawpos,Utrack_template);
    save(outfilename,'registeredpos');
    clear('rawpos','registeredpos');
end
        
