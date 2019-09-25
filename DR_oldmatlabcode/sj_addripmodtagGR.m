
function [task] = sj_addripmodtagGR(animno,animdirect,fileprefix)
% gideon:  [task] = sj_addripmodtagGR(1,'/data15/gideon/Ndl/','Ndl')
startDay=1;
if strcmp('Ndl',fileprefix), startDay=8;end

% sj_addripmodtag(1,'/data25/sjadhav/HPExpt/HPa_direct/','HPa');
% sj_addripmodtag(2,'/data25/sjadhav/HPExpt/HPb_direct/','HPb');
% sj_addripmodtag(3,'/data25/sjadhav/HPExpt/HPc_direct/','HPc');
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
    
    % Rip-modulation according to the variance test, calle ripmodtag2
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.ripmodtag2 = 'u';
    % The following tags will be propagated across epochs - Doing this only for PFC cells?
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.runripmodtag2 = 'u';
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.postsleepripmodtag2 = 'u';
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.presleepripmodtag2 = 'u';

end


procdatadir = '/mnt/data25/sjadhav/HPExpt/ProcessedData/';
datadirprefix=[animdirect '/' fileprefix]; 

%1) Do PFC cells first

%1a) RUN DATA FOR PFC
load([procdatadir, 'HP_ripplemod_PFC_alldata_gather']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches
% Gideon: changed a bit here. Going only from day 8, and identifying
% run/sleep epochs from the task struct rather from hardcoded numbers
for d=startDay:length(cellinfo)
    if d<10
     load([datadirprefix 'task0' num2str(d) '.mat'])
    else
     load([datadirprefix 'task' num2str(d) '.mat'])
    end
     numDailyEpochs= length(task{d});
    allepochs=1:numDailyEpochs;
        runepochs=[];
 for i=1:numDailyEpochs
     i
     if strcmp(task{d}{i}.type,'run')
         runepochs=[runepochs i];  
     end
 end

    % Gideon: should we really be iterating across epochs? I think the data
    % in allripplemod is already grouped across run epochs, right?? I think
    % it doesn't matter, just does it too many times.
    for e=1:length(runepochs)
        ep = runepochs(e);
        for tet=1:length(cellinfo{d}{ep})
            if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
                for c=1:length(cellinfo{d}{ep}{tet})
                    
                    currdaytetcell=[d tet c];
                    match = rowfind(currdaytetcell, daytetcell_list);
                    if d==15&&tet==17
                        keyboard
                    end
                    % There should be only 1 match. If no match, its undefined for that cell.
                    if match~=0
                        
                        getdataidx = animidxs(match);
                        % the regular test
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
                                % ERROR???
                                % was:
                                %cellinfo{d}{tmpe}{tet}{c}.runripmodtag='y';
                                % but I think it should be:
                                cellinfo{d}{tmpe}{tet}{c}.runripmodtag='n';
                            end
                            
                        end
                        % the new var test
                        if allripplemod(getdataidx).rasterShufP<0.05
                            cellinfo{d}{ep}{tet}{c}.ripmodtag2='y';
                            for tmpe=allepochs
                                tmpe;
                                cellinfo{d}{tmpe}{tet}{c}.runripmodtag2='y';
                            end
                            
                        else
                            cellinfo{d}{ep}{tet}{c}.ripmodtag2='n';
                            for tmpe=allepochs
                                tmpe;
                                % ERROR???
                                % was:
                                %cellinfo{d}{tmpe}{tet}{c}.runripmodtag='y';
                                % but I think it should be:
                                cellinfo{d}{tmpe}{tet}{c}.runripmodtag2='n';
                            end
                            
                        end
                        
                    end
                    
                end % end c
            end % end c
        end % end tet
    end % end ep
end %end d
save([animdirect, fileprefix,'cellinfo'], 'cellinfo');

'NOTE: NOT DOING sleep NOW!!'
return;


%1b) SLEEP DATA FOR PFC
% CHange for sleep - presleep (epoch 1) is now kept separate, so check separately
load([procdatadir, 'NdlBrg_ripplemodsleep_PFC_gather']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=startDay:length(cellinfo)
    if d<10
        load([datadirprefix 'task0' num2str(d) '.mat'])
    else
        load([datadirprefix 'task' num2str(d) '.mat'])
    end
    numDailyEpochs= length(task{d});
    allepochs=1:numDailyEpochs;
    slepochs=[];
    for i=1:numDailyEpochs
        if strcmp(task{d}{i}.type,'sleep')&strcmp(task{d}{i}.audprot,'0')
            slepochs=[slepochs i];
        end;
        
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
                            % regular test
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
                            
                            % var test
                             if allripplemod(getdataidx).rasterShufP<0.05
                                cellinfo{d}{ep}{tet}{c}.ripmodtag2='y';
                                for tmpe=allepochs
                                    cellinfo{d}{tmpe}{tet}{c}.postsleepripmodtag2='y';
                                end
                                    
                            else
                                cellinfo{d}{ep}{tet}{c}.ripmodtag2='n';
                                for tmpe=allepochs
                                    cellinfo{d}{tmpe}{tet}{c}.postsleepripmodtag2='n';
                                end
                            end
                            
                        
                            % Gideon: the pre-sleep should be empty for me
                        else % pre-sleep
                            %regular test
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
                            %var test
                            if allripplemod(getdataidx).rasterShufP<0.05
                                cellinfo{d}{ep}{tet}{c}.ripmodtag2='y';
                                for tmpe=allepochs
                                    cellinfo{d}{tmpe}{tet}{c}.presleepripmodtag2='y';
                                end
                            else
                                cellinfo{d}{ep}{tet}{c}.ripmodtag2='n';
                                for tmpe=allepochs
                                    cellinfo{d}{tmpe}{tet}{c}.presleepripmodtag2='n';
                                end
                            end
                            
                            
                        end
                        
                        
                        
                        
                        
                    end
                    
                end % end c
            end % end c
        end % end tet
    end % end ep
end %end d



'NOTE: NOT DOING HC CELLS NOW!!'
return;





%1) Do CA1 and iCA1 cells next

%1a) RUN DATA FOR Hipp
load([procdatadir, 'Ndl_ripplemod_CA1_gather']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=startDay:length(cellinfo)
    if d<10
     load([datadirprefix 'task0' num2str(d) '.mat'])
    else
     load([datadirprefix 'task' num2str(d) '.mat'])
    end
     numDailyEpochs= length(task{d});
    allepochs=1:numDailyEpochs;
        runepochs=[];
 for i=1:numDailyEpochs
     
     if strcmp(task{d}{i}.type,'run')
         runepochs=[runepochs i];  
     end
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

load([procdatadir, 'Ndl_ripplemodsleep_CA1_gather']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

for d=startDay:length(cellinfo)
    
   if d<10
        load([datadirprefix 'task0' num2str(d) '.mat'])
    else
        load([datadirprefix 'task' num2str(d) '.mat'])
    end
    numDailyEpochs= length(task{d});
    allepochs=1:numDailyEpochs;
    slepochs=[];
    for i=1:numDailyEpochs
        if strcmp(task{d}{i}.type,'sleep')&strcmp(task{d}{i}.audprot,'0')
            slepochs=[slepochs i];
        end;
        
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



