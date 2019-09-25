function showSelect(selectId)
%function showSelect(selectId)

global fmaux
if nargin < 1;
    selectId= fmaux.selectid
end
load([fmaux.data2dir '/select-' selectId]);
nC= size(select.cellnum,1);
nS= 0;

for iC=1:nC
    fprintf(1, '%2.d: [%d %d %d %d], ', iC, select.cellnum(iC,:));
    nA= length(select.a{iC});
    fprintf(1, '%d selections.\n', nA);
    for iA=1:nA
        if isfield(select.a{iC}{iA}, 'linpos')
            fprintf(1, '\t traj= %d, x=[%.2f %.2f].\n', select.a{iC}{iA}.traj, ...
            select.a{iC}{iA}.linpos);
        else
            fprintf(1, '\t traj= %d .\n', select.a{iC}{iA}.traj);
        end
        nS= nS+1;
    end
    fprintf(1, '\n');
end
fprintf(1, '---\ntotal selections= %d\n', nS);
