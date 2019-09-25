function err= forallrats(fctname, selectid, dirname, prefix)
%function err= forallrats(fctname, selectid, dirname)
%function err= forallrats(fctname, selectid, dirname, prefix)
%
% [prefix]= {'kyl', 'ter', 'sta', 'fel'};
%
% e.g: forallrats('runAdaptFilter', 'all', 'adapt')


if nargin < 3; error('not enough arguments given'); end
if nargin < 4;
    prefix={'kyl', 'ter', 'sta', 'fel'};
end

global fmaux % file manager auxillary variable

% set search paths
setRoot
if strfind(dirname, 'data');
    subdir=  '/';
else
    subdir=  '/results/';
end

olddir= pwd;
nRats= length(prefix);
err= 0;
for r=1:nRats
    % change directory, set local options
    name= [root '/' prefix{r} subdir dirname];
    if exist(name)==7
        cd(name);
        setLocalOptions;

        % select file
        fmaux.selectid= selectid;
        fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];
        
        % run command
        eval(fctname);
    else
        disp(['directory ' name ' does not exist.']);
        err= 1;
    end
end
cd(olddir)

