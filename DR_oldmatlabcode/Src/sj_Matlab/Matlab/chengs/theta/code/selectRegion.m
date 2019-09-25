function selectRegion(selectIn, selectOut, minx, maxx)
%function selectRegion(selectIn, selectOut, minx, maxx)
%
%  [20, 70] home arm
%  [90,140] outside arm

% parameters
cleanup=       1;



global fmaux linbehav linmeanvel behavdata nbinx select task tetloc spikes
fmaux.selectid= selectIn;
fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];


loadFile(fmaux.select);
loadVar(fmaux.datadir, 'tetloc', 0, fmaux.prefix, 1);
[d,e,t,c]= startCellList;
if ~isfield(select, 'x');
    select.x= {};
end

while ~isempty(d)
    n= fmaux.currentCell;
    loadVar(fmaux.datadir, 'spikes', d, fmaux.prefix, 1);

    if  n>length(select.a) | isempty(select.a{n}) | ...
        isempty(tetloc{d}) | isempty(tetloc{d}{t}) | ...
        (e~= 2) | ...                              % only familiar arm
        ~strcmp(tetloc{d}{t}.region, 'CA1') | ...  % only CA1 cells count
        ~strcmp(spikes{d}{e}{t}{c}.celltype, 'PE') % only excitatory cells
        % reject traj
        select.a{n}= {};
        select.x{n}= {};
    else 
        % accept traj
        nana= length(select.a{n});
        for ia=1:nana
            select.a{n}{ia}.linpos= [minx; maxx];
        end
%        keyboard
    end

    [d,e,t,c]= getNextCell;
end

if(cleanup)
    disp('**** clean up ****')
    select= cleanSelect(select);
end

save([fmaux.data2dir '/select-' selectIn '-' selectOut], 'select');
