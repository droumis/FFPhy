function loadFile(fname)
%function loadFile(fname)
% Loads file fname only if not already loaded

global fmaux select

found= 0;
nfiles= length(fmaux.loaded);
for k=1:nfiles
    if strcmp(fname, fmaux.loaded(k).var)
        found= k;
        break;
    end
end

if ~found
    disp(sprintf('loading %s ...', fname))
    load(fname);
    fmaux.loaded(nfiles+1).var= fname;
    fmaux.loaded(nfiles+1).day= 0;
end

