function [s, analist]= auxLoadStats(selectid, adaptid)
%function [s, analist]= auxLoadStats(selectid, adaptid)

s= {};
%% load stats
load(['/home/chengs/theta/data/analist-' selectid])
nC= length(analist.rat);
oldrat= '';
for iC= 1:nC
    % set up data
    rat= analist.rat{iC};
%    if ~strcmp(rat, 'kyl'); return; end %@@
    num= analist.cellnum(iC,:); d=num(1); e=num(2); t=num(3); c=num(4);
    % load real stats
    if ~strcmp(rat, oldrat)
        load(['/home/chengs/theta/' rat '/' adaptid '/stats-' analist.selectid '.mat']);
        oldrat= rat;
    end
    s{iC,1}= stats{d}{e}{t}{c};
end  % for iC
