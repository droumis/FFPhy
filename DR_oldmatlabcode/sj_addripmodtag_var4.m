
function sj_addripmodtag_var4(animno,animdirect,fileprefix)
% Add the var2 variable as ripplemodtag2: RE_INITALIZE ONLY var, or else this will erase everything previous
% ONLY Runs, and only PFC implemented right now. Comment Sleeps out for now

%sj_addripmodtag_var4(1,'/data25/sjadhav/HPExpt/HPa_direct/','HPa')
%sj_addripmodtag_var4(2,'/data25/sjadhav/HPExpt/HPb_direct/','HPb')
%sj_addripmodtag_var4(3,'/data25/sjadhav/HPExpt/HPc_direct/','HPc')
%sj_addripmodtag_var4(4,'/data25/sjadhav/HPExpt/Ndl_direct/','Ndl')

%sj_addripmodtag_var4(1,'/mnt/data25new/sjadhav/HPExpt/HPa_direct/','HPa')
%sj_addripmodtag_var4(2,'/mnt/data25new/sjadhav/HPExpt/HPb_direct/','HPb')
%sj_addripmodtag_var4(3,'/mnt/data25new/sjadhav/HPExpt/HPc_direct/','HPc')
%sj_addripmodtag_var4(4,'/data15/gideon/Ndl/','Ndl')


% Add tag to identify ripple modulated and non-modulated cells. Whether a cell is modulated or not is obtained from *ripplemod_gather files.
% For run epochs, currently, the significance is calculated by combining across run epochs in a day
% For sleep epochs, pre-sleep is kept separate from post-sleep epochs.

% Update - I need to propogate the sleepripmodtag to other epochs as well for theta cov analysis. 
% So make separate fields "postsleepripmodtag", "presleepripmodtag" to transfer to all epochs.
% Also do the same for "runripmodtag". Might be useful later.
% Note: Only doing this for PFC cells right now. "ripmodtag" for CA1 not that useful currently.

% Note that you are not using the "ripplemod" structure files, which contains the ripple-aligned responses for each epoch separately, 
% but you are not using statistics from that file.

% Start with a clean slate. Set ripplemodtag2 for all to Undefined = "u"
load([animdirect, fileprefix,'cellinfo']);
o = cellfetch(cellinfo,'numspikes');
targetcells = o.index;
for i = 1:size(targetcells,1)
    % ripmodtag will mark signficance for current epoch
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.ripmodtag2 = 'u';
    % The following tags will be propagated across epochs - Doing this only for PFC cells?
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.runripmodtag2 = 'u';
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.postsleepripmodtag2 = 'u';
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.presleepripmodtag2 = 'u';
end

procdatadir = '/data25/sjadhav/HPExpt/ProcessedData/';
%procdatadir = '/data15/gideon/ProcessedData/';

%1) Do PFC cells first

%1a) RUN DATA FOR PFC
load([procdatadir, 'HP_ripplemod_PFC_alldata_std3_speed4_ntet2_Nspk50_gather_2-12-2014']);
%load([procdatadir, 'HP_ripplemod_PFC_alldata_gather_var2']);
animidxs = find(allripplemod_idx(:,1)==animno);
daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);

% Instead of looping through daytetcell_list, loop through cellinfo and find matches

mindays = 1; maxdays = length(cellinfo);
if animno==4, 
    mindays = 8; maxdays = 17;
end

for d=mindays:maxdays
    
    if animno~=4 % If not Nadal
        if d==1,
            runepochs = [4 6]; % only Wtracks for now.
            allepochs = [1:7]; % To propagate runripmodtag
        else
            runepochs = [2 4];
            allepochs = [1:5]; % To propagate runripmodtag
        end
    else
        d % day starts from 8
        taskfile = sprintf('%s/%stask%02d.mat', animdirect, fileprefix, d);
        load(taskfile);
        allepochs = 1:length(task{d});
        tmpepochs = length(task{d}); 
        runepochs = [];
        for tmpep = 1:tmpepochs
            if strcmp(task{d}{tmpep}.type,'run') % If a run epoch
                runepochs = [runepochs; tmpep];
            end
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
   
                        if allripplemod(getdataidx).rasterShufP2 < 0.05
                            cellinfo{d}{ep}{tet}{c}.ripmodtag2='y';
                            for tmpe=allepochs
                                tmpe;
                                cellinfo{d}{tmpe}{tet}{c}.runripmodtag2='y';
                            end
                            
                        else
                            cellinfo{d}{ep}{tet}{c}.ripmodtag2='n';
                            for tmpe=allepochs
                                tmpe;
                                cellinfo{d}{tmpe}{tet}{c}.runripmodtag2='n';
                            end
                            
                        end
                    end
                    
                end % end c
            end % end c
        end % end tet
    end % end ep
end %end d


% 
% %1b) SLEEP DATA FOR PFC
% % CHange for sleep - presleep (epoch 1) is now kept separate, so check separately
% load([procdatadir, 'HP_ripplemodsleep_PFC_gather']);
% animidxs = find(allripplemod_idx(:,1)==animno);
% daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);
% 
% % Instead of looping through daytetcell_list, loop through cellinfo and find matches
% 
% for d=mindays:maxdays
%     
%   if animno~=4 % If not Nadal
%        if d==1,
%             slepochs = [1 7];
%           allepochs = [1:7]; % To propogate post/pre/sleepripmodtag
%        else
%             slepochs = [1 5];
%             allepochs = [1:5];
%        end
%   else
%       d % day starts from 7
%       taskfile = sprintf('%s/%stask%02d.mat', animdirect, fileprefix, d);
%       load(taskfile);
%       allepochs = 1:length(task{d});
%       tmpepochs = length(task{d});
%       runepochs = [];
%       for tmpep = 1:tmpepochs
%           if strcmp(task{d}{tmpep}.type,'sleep') % If a run epoch
%               runepochs = [runepochs; tmpep];
%           end
%       end
%   end

%     
%     for e=1:length(slepochs)
%         ep = slepochs(e);
%         for tet=1:length(cellinfo{d}{ep})
%             if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
%                 for c=1:length(cellinfo{d}{ep}{tet})
%                     
%                     currdaytetcell=[d tet c];
%                     match = rowfind(currdaytetcell, daytetcell_list);
%                     
%                     % There should be only 1 match. If no match, its undefined for that cell.
%                     if match~=0
%                         
%                         getdataidx = animidxs(match);
%                         
%                         if ep~=1 % not pre-sleep
%                             
%                             if allripplemod(getdataidx).sig_shuf==1
%                                 cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
%                                 for tmpe=allepochs
%                                     cellinfo{d}{tmpe}{tet}{c}.postsleepripmodtag='y';
%                                 end
%                                     
%                             else
%                                 cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
%                                 for tmpe=allepochs
%                                     cellinfo{d}{tmpe}{tet}{c}.postsleepripmodtag='n';
%                                 end
%                             end
%                             
%                         else % pre-sleep
%                             
%                             if allripplemod(getdataidx).sig_shuf_pre==1
%                                 cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
%                                 for tmpe=allepochs
%                                     cellinfo{d}{tmpe}{tet}{c}.presleepripmodtag='y';
%                                 end
%                             else
%                                 cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
%                                 for tmpe=allepochs
%                                     cellinfo{d}{tmpe}{tet}{c}.presleepripmodtag='n';
%                                 end
%                             end
%                             
%                         end
%                         
%                     end
%                     
%                 end % end c
%             end % end c
%         end % end tet
%     end % end ep
% end %end d






% %1) Do CA1 and iCA1 cells next
% 
% %1a) RUN DATA FOR Hipp
% load([procdatadir, 'HP_ripplemod_CA1_gather']);
% animidxs = find(allripplemod_idx(:,1)==animno);
% daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);
% 
% % Instead of looping through daytetcell_list, loop through cellinfo and find matches
% 
% for d=mindays:maxdays
%     
%     if animno~=4 % If not Nadal
%         if d==1,
%             runepochs = [4 6]; % only Wtracks for now.
%             allepochs = [1:7]; % To propagate runripmodtag
%         else
%             runepochs = [2 4];
%             allepochs = [1:5]; % To propagate runripmodtag
%         end
%     else
%         d % day starts from 7
%         taskfile = sprintf('%s/%stask%02d.mat', animdirect, fileprefix, d);
%         load(taskfile);
%         allepochs = 1:length(task{d});
%         tmpepochs = length(task{d}); 
%         runepochs = [];
%         for tmpep = 1:tmpepochs
%             if strcmp(task{d}{tmpep}.type,'run') % If a run epoch
%                 runepochs = [runepochs; tmpep];
%             end
%         end
%    end
%     
%     for e=1:length(runepochs)
%         ep = runepochs(e);
%         for tet=1:length(cellinfo{d}{ep})
%             if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
%                 for c=1:length(cellinfo{d}{ep}{tet})
%                     
%                     currdaytetcell=[d tet c];
%                     match = rowfind(currdaytetcell, daytetcell_list);
%                     
%                     % There should be only 1 match. If no match, its undefined for that cell.
%                     if match~=0
%                         
%                         getdataidx = animidxs(match);
%                         
%                         if allripplemod(getdataidx).rasterShufP2 < 0.05
%                             cellinfo{d}{ep}{tet}{c}.ripmodtag2='y';
%                         else
%                             cellinfo{d}{ep}{tet}{c}.ripmodtag2='n';
%                         end
%                     end
%                     
%                 end % end c
%             end % end c
%         end % end tet
%     end % end ep
% end %end d


% 
% %1b) SLEEP DATA FOR Hipp
% % CHange for sleep - presleep (epoch 1) is now kept separate, so check separately
% 
% load([procdatadir, 'HP_ripplemodsleep_CA1_gather']);
% animidxs = find(allripplemod_idx(:,1)==animno);
% daytetcell_list = allripplemod_idx(animidxs,[2 3 4]);
% 
% % Instead of looping through daytetcell_list, loop through cellinfo and find matches
% 
% for d=mindays:maxdays
%     
%   if animno~=4 % If not Nadal
%        if d==1,
%             slepochs = [1 7];
%           allepochs = [1:7]; % To propogate post/pre/sleepripmodtag
%        else
%             slepochs = [1 5];
%             allepochs = [1:5];
%        end
%   else
%       d % day starts from 7
%       taskfile = sprintf('%s/%stask%02d.mat', animdirect, fileprefix, d);
%       load(taskfile);
%       allepochs = 1:length(task{d});
%       tmpepochs = length(task{d});
%       runepochs = [];
%       for tmpep = 1:tmpepochs
%           if strcmp(task{d}{tmpep}.type,'sleep') % If a run epoch
%               runepochs = [runepochs; tmpep];
%           end
%       end
%   end
%     
%     for e=1:length(slepochs)
%         ep = slepochs(e);
%         for tet=1:length(cellinfo{d}{ep})
%             if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
%                 for c=1:length(cellinfo{d}{ep}{tet})
%                     
%                     currdaytetcell=[d tet c];
%                     match = rowfind(currdaytetcell, daytetcell_list);
%                     
%                     % There should be only 1 match. If no match, its undefined for that cell.
%                     if match~=0
%                         
%                         getdataidx = animidxs(match);
%                         
%                         if ep~=1 % not pre-sleep
%                             
%                             if allripplemod(getdataidx).sig_shuf==1
%                                 cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
%                             else
%                                 cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
%                             end
%                             
%                         else % pre-sleep
%                             
%                             if allripplemod(getdataidx).sig_shuf_pre==1
%                                 cellinfo{d}{ep}{tet}{c}.ripmodtag='y';
%                             else
%                                 cellinfo{d}{ep}{tet}{c}.ripmodtag='n';
%                             end
%                             
%                         end
%                     end
%                     
%                 end % end c
%             end % end c
%         end % end tet
%     end % end ep
% end %end d


i=1;
% Save updated cellinfo file
save([animdirect, fileprefix,'cellinfo'], 'cellinfo');



