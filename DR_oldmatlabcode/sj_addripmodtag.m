
function [task] = sj_addripmodtag(animno,animdirect,fileprefix)

% sj_addripmodtag(1,'/mnt/data25/sjadhav/HPExpt/HPa_direct/','HPa');
% sj_addripmodtag(2,'/mnt/data25/sjadhav/HPExpt/HPb_direct/','HPb');
% sj_addripmodtag(3,'/mnt/data25/sjadhav/HPExpt/HPc_direct/','HPc');
% Shantanu
% Add tag to identify ripple modulated and non-modulated cells. Whether a cell is modulated or not is obtained from *ripplemod_gather files.
% For run epochs, currently, the significance is calculated by combining across run epochs in a day
% For sleep epochs, pre-sleep is kept separate from post-sleep epochs.

% Update - I need to propogate the sleepripmodtag to other epochs as well for theta cov analysis. 
% So make separate fields "postsleepripmodtag", "presleepripmodtag" to transfer to all epochs.
% Also do the same for "runripmodtag". Might be useful later.
% Note: Only doing this for PFC cells right now. "ripmodtag" for CA1 not that useful currently.

% Note that you are not using the "ripplemod" structure files, which contains the ripple-aligned responses for each epoch separately, 
% but you are not using statistics from that file.

% Start with a clean slate. Set tag for all to Undefined = "u"
load([animdirect, fileprefix,'cellinfo']);
o = cellfetch(cellinfo,'numspikes');
targetcells = o.index;
for i = 1:size(targetcells,1)
    % ripmodtag will mark signficance for current epoch
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.ripmodtag = 'u';
    % The following tags will be propagated across epochs - Doing this only for PFC cells?
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.runripmodtag = 'u';
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.postsleepripmodtag = 'u';
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.presleepripmodtag = 'u';
end


procdatadir = '/mnt/data25/sjadhav/HPExpt/ProcessedData/';

%1) Do PFC cells first

%1a) RUN DATA FOR PFC
load([procdatadir, 'HP_ripplemod_PFC_alldata_gather']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=1:length(cellinfo)
    
    if d==1,
        runepochs = [4 6]; % only Wtracks for now.
        allepochs = [1:7]; % To propagate runripmodtag
    else
        runepochs = [2 4];
        allepochs = [1:5]; % To propagate runripmodtag
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
                        
                        if allripplemod(getdataidx).sig_shuf==1
                            cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
                            for tmpe=allepochs
                                tmpe;
                                cellinfo{d}{tmpe}{tet}{c}.runripmodtag='y';
                            end
                            
                        else
                            cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
                            for tmpe=allepochs
                                tmpe;
                                cellinfo{d}{tmpe}{tet}{c}.runripmodtag='n';
                            end
                            
                        end
                    end
                    
                end % end c
            end % end c
        end % end tet
    end % end ep
end %end d



%1b) SLEEP DATA FOR PFC
% CHange for sleep - presleep (epoch 1) is now kept separate, so check separately
load([procdatadir, 'HP_ripplemodsleep_PFC_gather']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=1:length(cellinfo)
    
    if d==1,
        slepochs = [1 7];
        allepochs = [1:7]; % To propogate post/pre/sleepripmodtag
    else
        slepochs = [1 5];
        allepochs = [1:5];
    end
    
    for e=1:length(slepochs)
        ep = slepochs(e);
        for tet=1:length(cellinfo{d}{ep})
            if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
                for c=1:length(cellinfo{d}{ep}{tet})
                    
                    currdaytetcell=[d tet c];
                    match = rowfind(currdaytetcell, daytetcell_list);
                    
                    % There should be only 1 match. If no match, its undefined for that cell.
                    if match~=0
                        
                        getdataidx = animidxs(match);
                        
                        if ep~=1 % not pre-sleep
                            
                            if allripplemod(getdataidx).sig_shuf==1
                                cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
                                for tmpe=allepochs
                                    cellinfo{d}{tmpe}{tet}{c}.postsleepripmodtag='y';
                                end
                                    
                            else
                                cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
                                for tmpe=allepochs
                                    cellinfo{d}{tmpe}{tet}{c}.postsleepripmodtag='n';
                                end
                            end
                            
                        else % pre-sleep
                            
                            if allripplemod(getdataidx).sig_shuf_pre==1
                                cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
                                for tmpe=allepochs
                                    cellinfo{d}{tmpe}{tet}{c}.presleepripmodtag='y';
                                end
                            else
                                cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
                                for tmpe=allepochs
                                    cellinfo{d}{tmpe}{tet}{c}.presleepripmodtag='n';
                                end
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
load([procdatadir, 'HP_ripplemod_CA1_gather']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=1:length(cellinfo)
    
    if d==1,
        runepochs = [4 6];
    else
        runepochs = [2 4];
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
                        
                        if allripplemod(getdataidx).sig_shuf==1
                            cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
                        else
                            cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
                        end
                    end
                    
                end % end c
            end % end c
        end % end tet
    end % end ep
end %end d



%1b) SLEEP DATA FOR Hipp
% CHange for sleep - presleep (epoch 1) is now kept separate, so check separately

load([procdatadir, 'HP_ripplemodsleep_CA1_gather']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=1:length(cellinfo)
    
    if d==1,
        slepochs = [1 7];
    else
        slepochs = [1 5];
    end
    
    for e=1:length(slepochs)
        ep = slepochs(e);
        for tet=1:length(cellinfo{d}{ep})
            if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
                for c=1:length(cellinfo{d}{ep}{tet})
                    
                    currdaytetcell=[d tet c];
                    match = rowfind(currdaytetcell, daytetcell_list);
                    
                    % There should be only 1 match. If no match, its undefined for that cell.
                    if match~=0
                        
                        getdataidx = animidxs(match);
                        
                        if ep~=1 % not pre-sleep
                            
                            if allripplemod(getdataidx).sig_shuf==1
                                cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
                            else
                                cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
                            end
                            
                        else % pre-sleep
                            
                            if allripplemod(getdataidx).sig_shuf_pre==1
                                cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
                            else
                                cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
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



