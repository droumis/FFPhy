
function [task] = sj_addthetamodtag(animno,animdirect,fileprefix)

% Similar to sj_addripmodtag: Addsa tag based on significance of theta phase locking.
% sj_addthetamodtag(1,'/data25/sjadhav/HPExpt/HPa_direct/','HPa');
% sj_addthetamodtag(2,'/data25/sjadhav/HPExpt/HPb_direct/','HPb');
% sj_addthetamodtag(3,'/data25/sjadhav/HPExpt/HPc_direct/','HPc');
% Shantanu
% Add tag to identify theta modulated and non-modulated cells. Whether a cell is modulated or not is obtained from *thetamod_gather files.
% Only for run epochs. Currently, the significance is calculated by combining across run epochs in a day

% Update - Propagate thetamodtag to all epochs in the day. Call it runthetamodtag.


% Start with a clean slate. Set tag for all to Undefined = "u"
load([animdirect, fileprefix,'cellinfo']);
o = cellfetch(cellinfo,'numspikes');
targetcells = o.index;
for i = 1:size(targetcells,1)
    % thetamodtag will mark significance for current epoch: for runs only
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.thetamodtag = 'u';
    % The following tag will be propagated across epochs: both PFC abd CA1 cells
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.runthetamodtag = 'u';

end


procdatadir = '/data25/sjadhav/HPExpt/ProcessedData/';


%1) Do PFC cells first

%1a) RUN DATA FOR PFC
load([procdatadir, 'HP_thetamod_PFC_gather']);
animidxs = find(allthetamod_idx(:,1)==animno);
daytetcell_list = allthetamod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=1:length(cellinfo)
    
    if d==1,
        runepochs = [4 6];
        allepochs = [1:7]; % To propagate runthetamodtag
    else
        runepochs = [2 4];
        allepochs = [1:5]; % To propagate runthetamodtag
    end
    
    for e=1:length(runepochs)
        ep = runepochs(e);
        for tet=1:length(cellinfo{d}{ep})
            if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
                for c=1:length(cellinfo{d}{ep}{tet})
                    
                    currdaytetcell=[d tet c];
                    match = rowfind(currdaytetcell, daytetcell_list);
                    
                    % There should be only 1 match. If no match, its undefined for that cell.
                    if match~=0
                        
                        getdataidx = animidxs(match);
                        
                        if allthetamod(getdataidx).prayl < 0.05
                            cellinfo{d}{ep}{tet}{c}.thetamodtag='y';
                            for tmpe=allepochs
                                cellinfo{d}{tmpe}{tet}{c}.runthetamodtag='y';
                            end
                        else
                            cellinfo{d}{ep}{tet}{c}.thetamodtag='n';
                            for tmpe=allepochs
                                cellinfo{d}{tmpe}{tet}{c}.runthetamodtag='n';
                            end
                        end
                    end
                    
                end % end c
            end % end c
        end % end tet
    end % end ep
end %end d




%1) Do CA1 and iCA1 cells next

%1a) RUN DATA FOR Hipp
load([procdatadir, 'HP_thetamod_CA1_gather']);
animidxs = find(allthetamod_idx(:,1)==animno);
daytetcell_list = allthetamod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=1:length(cellinfo)
    
    if d==1,
        runepochs = [4 6];
        allepochs = [1:7]; % To propagate runthetamodtag
    else
        runepochs = [2 4];
        allepochs = [1:5];
    end
    
    for e=1:length(runepochs)
        ep = runepochs(e);
        for tet=1:length(cellinfo{d}{ep})
            if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
                for c=1:length(cellinfo{d}{ep}{tet})
                    
                    currdaytetcell=[d tet c];
                    match = rowfind(currdaytetcell, daytetcell_list);
                    
                    % There should be only 1 match. If no match, its undefined for that cell.
                    if match~=0
                        
                        getdataidx = animidxs(match);
                        
                        if allthetamod(getdataidx).prayl < 0.05
                            cellinfo{d}{ep}{tet}{c}.thetamodtag='y';
                            for tmpe=allepochs
                                cellinfo{d}{tmpe}{tet}{c}.runthetamodtag='y';
                            end
                        else
                            cellinfo{d}{ep}{tet}{c}.thetamodtag='n';
                            for tmpe=allepochs
                                cellinfo{d}{tmpe}{tet}{c}.runthetamodtag='n';
                            end
                        end
                    end
                    
                end % end c
            end % end c
        end % end tet
    end % end ep
end %end d





i=1;
% Save updated cellinfo file
save([animdirect, fileprefix,'cellinfo'], 'cellinfo');



