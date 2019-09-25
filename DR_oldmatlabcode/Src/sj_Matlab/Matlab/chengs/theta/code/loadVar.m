function [success,loaded]= loadVar(dir,varname,day,prefix,rename,filename)
%function [success,loaded]= loadVar(dir,varname,day,prefix,rename,filename)
%
% day= 0 means there is one file for whole experiment
%

silent= 0;

global fmaux

loaded= false;
success= 1;

if nargin < 3; day= 0; end
if nargin < 4
    prefix= [];
end
if nargin < 5
    rename= false;
end

invarname= [prefix varname];
if rename
    outvarname= varname;
else
    outvarname= invarname;
end

if nargin<6
    if day~= 0
        fname= sprintf('%s/%s%.2d.mat',dir,invarname,day);
    else
        fname= sprintf('%s/%s.mat',dir,invarname);
    end
else
    fname= sprintf('%s/%s.mat',dir,filename);
end

if strcmp(dir, '.')
    fname= [pwd fname(2:end)];
elseif fname(1)~='/'; 
    fname= which(fname); 
end

found= 0;
loadfile= 1;
nfiles= 0;
if isfield(fmaux, 'loaded')
    nfiles= length(fmaux.loaded);
    for k=1:nfiles
        if strcmp(outvarname, fmaux.loaded(k).var)
            found= k;
            if strcmp(fname, fmaux.loaded(k).fname) & day==fmaux.loaded(k).day
                loadfile= 0;
            end
            break;
        end
    end
end

if loadfile
    eval(['global ' outvarname]);

    if exist(fname, 'file')
        if ~silent; disp(sprintf('loading %s ...', fname)); end
        load(fname);
        loaded= true;
    else
        disp(sprintf('File %s does not exist (,yet).', fname))
        success= 0;
    end
    if found
        fmaux.loaded(found).fname= fname;
        fmaux.loaded(found).day= day;
    elseif success
        fmaux.loaded(nfiles+1).fname= fname;
        fmaux.loaded(nfiles+1).var= outvarname;
        fmaux.loaded(nfiles+1).day= day;
    end

    if rename & success
        eval([outvarname '= ' invarname ';']);
    end
end

