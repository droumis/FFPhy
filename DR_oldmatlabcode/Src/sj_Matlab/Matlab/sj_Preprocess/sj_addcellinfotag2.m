
function [task] = sj_addcellinfotag2(animdirect,fileprefix)

% Shantanu
% Add tag2 to cell info based on tag given during matclust, or using meanfiringrate, etc

load([animdirect, fileprefix,'cellinfo']);

% Dont want to use cellfetch like following, since you want to examine tag and will have to write a loop
% o = cellfetch(cellinfo,'tag');
% targetcells = o.index(find(ismember(o.index(:,1),days)),:);

for d=1:length(cellinfo)
    for ep=1:length(cellinfo{d})
        for tet=1:length(cellinfo{d}{ep})
            if ~isempty(cellinfo{d}{ep}{tet}) % if cells for the tet
                for c=1:length(cellinfo{d}{ep}{tet})
                    if ~isfield(cellinfo{d}{ep}{tet}{c},'tag2') % Only proceed if tag2 is not in there
                        if isfield(cellinfo{d}{ep}{tet}{c},'tag') % check if tag exists

                            currtag = cellinfo{d}{ep}{tet}{c}.tag;
                            if iscell(currtag)
                                currtag = currtag{1};
                            end                      
                            if length(currtag)<3, % if old MU tag is still in, or for empty but not deleted clusters, call it Undef
                               d,ep,tet,c, keyboard; 
                               cellinfo{d}{ep}{tet}{c}.tag2 = 'Undef';
                               cellinfo{d}{ep}{tet}{c}.tag = 'Undef';
                               currtag = 'Undef';
                            end
                            
                            %PFC
                            if strcmp(currtag(1:3),'PFC')
                                cellinfo{d}{ep}{tet}{c}.tag2 = 'PFC'; %Call all PFC, regardless pf PFCSil, PFCRun, PFCpr, PFCps, PFCp
                            end
                            % CA1
                            if strcmp(currtag(1:3),'CA1')
                                if strcmp(currtag,'CA1Int')
                                    cellinfo{d}{ep}{tet}{c}.tag2 = 'CA1Int'; % If tag is Int, tag2 is Int
                                else
                                    cellinfo{d}{ep}{tet}{c}.tag2 = 'CA1Pyr'; %all else is CA1Pyr, regardless of CA1Sil, CA1Run, etc
                                end
                            end
                            % iCA1
                            if strcmp(currtag(1:3),'iCA')
                                if strcmp(currtag,'iCA1Int')
                                    cellinfo{d}{ep}{tet}{c}.tag2 = 'iCA1Int'; % If tag is Int, tag2 is Int
                                else
                                    cellinfo{d}{ep}{tet}{c}.tag2 = 'iCA1Pyr'; %all else is CA1Pyr, regardless of CA1Sil, CA1Run, etc
                                end
                            end
                            
                        else % If cell exists but no tag field
                            if isfield(cellinfo{d}{ep}{tet}{c},'tag') % Check if area field exists
                                area = cellinfo{d}{ep}{tet}{c}.area;
                                %PFC
                                if strcmp(area,'PFC')
                                    cellinfo{d}{ep}{tet}{c}.tag2 = 'PFC';
                                end
                                %CA1
                                if strcmp(area,'CA1')
                                    if cellinfo{d}{ep}{tet}{c}.meanrate <= 7 % meanrate less than 7Hz for Hippocampus
                                        cellinfo{d}{ep}{tet}{c}.tag2 = 'CA1Pyr';
                                    else
                                        cellinfo{d}{ep}{tet}{c}.tag2 = 'CA1Int';
                                    end
                                end
                                %iCA1
                                if strcmp(area,'iCA1')
                                    if cellinfo{d}{ep}{tet}{c}.meanrate <= 7 % meanrate less than 7Hz for Hippocampus
                                        cellinfo{d}{ep}{tet}{c}.tag2 = 'iCA1Pyr';
                                    else
                                        cellinfo{d}{ep}{tet}{c}.tag2 = 'iCA1Int';
                                    end
                                end
                            else % if tag and area field does not exist, call it Undefined. Go to Keyboard.
                                cellinfo{d}{ep}{tet}{c}.tag2 = 'Undef';
                                %keyboard;
                            end
                            
                        end % end isfield tag
                    end % end Check for existence of tag2
                end % end cell
            end % end if tet has cells
        end % end tet
    end % end ep
end % end day


% Save updated cellinfo file
save([animdirect, fileprefix,'cellinfo'], 'cellinfo')