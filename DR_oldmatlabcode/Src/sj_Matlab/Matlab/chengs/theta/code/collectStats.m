function [mAll,mIndiv,sAll,sIndiv]= collectStats(prefix, varname, indices, minspikes)
%function [mAll,mIndiv,sAll,sIndiv]= collectStats(prefix, varname, indices, minspikes)
%  
%  Collect summary of statistic for several animals.
%  Each statistic is usually computed for each cell at several time points and
%  saved per cell, i.e. stat{d}{e}{t}{c}.xyz= [1 2 3].
%  This function collects the statistics across all cells at particular time
%  indices. 
%  
% Output:
%  sAll{i}   stat across all animals.
%  sIndividual{rat}{i}   stat for particular animal.
%  mAll      summary (mean, std) across all animals.
%  mIndiv    summary (mean, std) for particular animal.
%
% Input:
%  prefix    animals for which to collect statistics
%  varname   name of statistics
%  indices   time points at which 

debugging= 1;

nindices= length(indices);
nRats= length(prefix);
olddir= pwd;
mAll= [];
mIndiv= cell(nRats,1);
sAll= cell(nindices,1);
sIndiv= cell(nRats,1);
for r=1:nRats
    % directories
    cd(['/home/chengs/theta/' prefix{r} '/adapt']);
    setLocalOptions;

    fmaux.selectid= 'all-placefields';
    fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];

    sIndiv{r}= cell(nindices,1);
    if nargin < 4
        stmp= extractStats(varname,indices);
    else
        stmp= extractStats(varname,indices,minspikes);
    end
    for p= indices
        if ~isfield(stmp{p},'val'); continue; end
        sIndiv{r}{p}= [sIndiv{r}{p},stmp{p}.val];
        sAll{p}= [sAll{p},stmp{p}.val];
    end
    mIndiv{r}= auxAverage(sIndiv{r});
end
mAll= auxAverage(sAll);
cd(olddir)

if debugging
    data= [indices; mAll.mean(indices)];
    fstring= 'index %2d: all= %.2f';
    for r=1:nRats
        fstring= [fstring ', ' prefix{r} '= %.2f'];
        data= [data; mIndiv{r}.mean(indices)];
    end
    fstring= [fstring, '\n'];
    fprintf(1, fstring, data);
end

function s= auxAverage(data)
s=[];
for p= 1:length(data)
    if isempty(data{p}); 
        s.mean(p)= nan;
        s.std(p)= nan;
        s.n(p)= nan;
        continue; 
    end
    s.mean(p)= mean(data{p});
    s.std(p)= std(data{p});
    s.n(p)= length(data{p});
end
